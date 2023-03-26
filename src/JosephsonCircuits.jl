__precompile__(true)

module JosephsonCircuits

import Graphs
import FFTW
import KLU
import UUIDs
import Symbolics: Symbolic, Sym, Num, @variables, @syms, @register_symbolic
import Symbolics
import AxisKeys
import SnoopPrecompile
import OrderedCollections
import Statistics

using LinearAlgebra
using SparseArrays
using Printf

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
include("fftutils.jl")
include("qesparams.jl")
include("keyedarrayutils.jl")

# old single tone harmonic balance solver code
include("hbsolveold.jl")

# new mult-tone harmonic balance solver code
include("hbsolve.jl")

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
    push!(circuit,("R1","1","0",:Rleft))
    push!(circuit,("C1","1","2",:Cc)) 
    push!(circuit,("Lj1","2","0",:Lj)) 
    push!(circuit,("C2","2","0",:Cj))

    circuitdefs = Dict{Symbol,Complex{Float64}}(
        :Lj =>1000.0e-12,
        :Cc => 100.0e-15,
        :Cj => 1000.0e-15,
        :Rleft => 50.0,
    )

    return hbsolve(2*pi*(4.5:0.5:5.0)*1e9,2*pi*4.75001*1e9,0.00565e-6,2,2,circuit,circuitdefs,pumpports=[1]);

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

function warmupsymsnew()

    @variables R Cc Lj Cj
    circuit = [
        ("P1","1","0",1),
        ("R1","1","0",R),
        ("C1","1","2",Cc),
        ("Lj1","2","0",Lj),
        ("C2","2","0",Cj)]

    circuitdefs = Dict(
        Lj =>1000.0e-12,
        Cc => 100.0e-15,
        Cj => 1000.0e-15,
        R => 50.0)

    ws = 2*pi*(4.5:0.5:5.0)*1e9
    wp = (2*pi*4.75001*1e9,)
    sources = [(mode=(1,),port=1,current=0.00565e-6)]
    Nmodulationharmonics = (2,)
    Npumpharmonics = (4,)

    return hbsolve(ws, wp, sources, Nmodulationharmonics,
        Npumpharmonics, circuit, circuitdefs)
end



function warmupsymsold()

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

    return hbsolveold(2*pi*(4.5:0.5:5.0)*1e9, 2*pi*4.75001*1e9, 0.00565e-6, 2, 2, circuit, circuitdefs, pumpports=[1]);

end

function warmupparse()

    @variables Rleft Cc Lj Cj w L1
    # define the circuit components
    circuit = Array{Tuple{String,String,String,Num},1}(undef,0)

    # port on the left side
    push!(circuit,("P1","1","0",1))
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
    )

    return parsecircuit(circuit)
end

function warmupparsesort()

    @variables Rleft Cc Lj Cj w L1
    # define the circuit components
    circuit = Array{Tuple{String,String,String,Num},1}(undef,0)

    # port on the left side
    push!(circuit,("P1","1","0",1))
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
    )

    return parsesortcircuit(circuit)
end

function warmupnumericmatrices()

    @variables Rleft Cc Lj Cj w L1
    # define the circuit components
    circuit = Array{Tuple{String,String,String,Num},1}(undef,0)

    # port on the left side
    push!(circuit,("P1","1","0",1))
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
    )

    return numericmatrices(circuit,circuitdefs)
end

function warmuphblinsolve()

    @variables Rleft Cc Lj Cj w L1
    # define the circuit components
    circuit = Array{Tuple{String,String,String,Num},1}(undef,0)

    # port on the left side
    push!(circuit,("P1","1","0",1))
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
    )

    return hblinsolve(2*pi*(4.5:0.1:5.0)*1e9,circuit,circuitdefs)
end

function warmupvvn()

    @variables Rleft Cc Lj Cj w L1
    # define the circuit components
    circuit = Array{Tuple{String,String,String,Num},1}(undef,0)

    # port on the left side
    push!(circuit,("P1","1","0",1))
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
    )

    psc = parsesortcircuit(circuit,sorting=:number)

    return componentvaluestonumber(psc.componentvalues,circuitdefs)
end


export @syms, hbsolve, hbnlsolve, hblinsolve, hbsolveold, hbnlsolveold,
    hblinsolveold, parsecircuit, parsesortcircuit, calccircuitgraph,
    symbolicmatrices, numericmatrices, LjtoIc, IctoLj, @variables,
    @register_symbolic, Num, Symbolics


# the below precompile directives are to help the compiler perform type inference
# during the precompilation stage (when the package is installed) instead of
# when it is loaded or the functions are run. this is a helpful guide:
# https://timholy.github.io/SnoopCompile.jl/stable/snoopi_deep_analysis/#inferrability
# and the basic commands to look at the inference triggers
# julia> using SnoopCompile, JosephsonCircuits
# julia> tinf = @snoopi_deep JosephsonCircuits.warmupsyms();
# julia> itrigs = inference_triggers(tinf)

SnoopPrecompile.@precompile_all_calls begin
    warmup()
    warmupsyms()
    warmupsymsold()
    warmupsymsnew()
end

#end module
end
