__precompile__(true)

module JosephsonCircuits

import Graphs
import FFTW
import KLU
import UUIDs
import Symbolics: Symbolic, Sym, Num, @variables, @syms, @register_symbolic, @wrapped
import Symbolics
import SymbolicUtils
import AxisKeys
import PrecompileTools
import OrderedCollections
import StaticArrays
import Statistics

using LinearAlgebra
using SparseArrays
using Printf

# define the zero for symbolic numbers so that we can view the sparse arrays
# Base.zero(::Type{Symbolic{Number}}) = 0
# Base.zero(::Type{Any}) = 0
# Base.zero(::Type{Nothing}) = 0


"""
    const phi0

A constant for phi0, the reduced magnetic flux quantum in Weber, H*A:
phi0 = hbar/(2*charge of electron).
"""
const phi0 = 3.29105976e-16

"""
    const Phi0

A constant for Phi0, the magnetic flux quantum in Weber, H*A:
Phi0 = h/(2*charge of electron).
"""
const Phi0 = 2.067833848e-15

"""
    const speed_of_light

A constant for the speed of light which is 2.99792458e8 m/s
"""
const speed_of_light = 2.99792458e8

# The file structure below mimics the typical analysis flow. An input file
# is first parsed, then the incidence matrix is calculated, then the 
# capacitance, inverse inductance, and other matrices are calculated, then
# then circuit is solved using the harmonic balance method, and postprocessed
# to determine the scattering parameters and quantum efficiency. 
include("parseinput.jl")
include("graphproc.jl")
include("capindmat.jl")
include("matutils.jl")

include("networkparamconversion.jl")
include("networks.jl")
include("networkconnection.jl")

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
    # return hbsolve(2*pi*(4.5:0.5:5.0)*1e9,(2*pi*4.75001*1e9,),[(mode=(1,),port=1,current=0.00565e-6)],(2,),(2,),circuit,circuitdefs;keyedarrays=Val(false));

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
    # return hbsolve(2*pi*(4.5:0.5:5.0)*1e9,(2*pi*4.75001*1e9,),[(mode=(1,),port=1,current=0.00565e-6)],(2,),(2,),circuit,circuitdefs);
end

function warmupsymsnew()

    @variables R Cc Lj Cj
    circuit = Tuple{String,String,String,Num}[
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

function warmupnetwork()
    # StoZ, StoY, StoA, StoB, StoABCD
    # the different functions we want to test
    for f in [
            (JosephsonCircuits.ZtoS,JosephsonCircuits.StoZ),
            (JosephsonCircuits.YtoS,JosephsonCircuits.StoY),
            (JosephsonCircuits.AtoS,JosephsonCircuits.StoA),
            (JosephsonCircuits.BtoS,JosephsonCircuits.StoB),
            (JosephsonCircuits.ABCDtoS,JosephsonCircuits.StoABCD),
        ]
        # single matrix input
        for portimpedances in [
                rand(Complex{Float64}), rand(Complex{Float64},2),
            ]
            for arg1 in [rand(Complex{Float64},2,2), (StaticArrays.@MMatrix rand(Complex{Float64},2,2))]
                f[1](arg1,portimpedances=portimpedances)
                f[2](arg1,portimpedances=portimpedances)
                f[1](arg1)
                f[2](arg1)
            end
        end
        # array input
        for portimpedances in [rand(Complex{Float64}), rand(Complex{Float64},2,10)]
            for arg1 in [rand(Complex{Float64},2,2,10)]
                f[1](arg1,portimpedances=portimpedances)
                f[2](arg1,portimpedances=portimpedances)
                f[1](arg1)
                f[2](arg1)
            end
        end
        # vector of matrices
        for portimpedances in [rand(Complex{Float64}), rand(Complex{Float64},2) ]
            for arg1 in [
                    [rand(Complex{Float64},2,2) for i in 1:10],
                ]
                [f[1](arg1[i],portimpedances=portimpedances) for i in 1:10]
                [f[2](arg1[i],portimpedances=portimpedances) for i in 1:10]
                [f[1](arg1[i]) for i in 1:10]
                [f[2](arg1[i]) for i in 1:10]
            end
        end
    end

    # StoT, AtoB, ZtoA, YtoA, YtoB, ZtoB, ZtoY
    # the different functions we want to test
    for f in [
            (JosephsonCircuits.StoT,JosephsonCircuits.TtoS),
            (JosephsonCircuits.AtoB,JosephsonCircuits.BtoA),
            (JosephsonCircuits.ZtoA,JosephsonCircuits.AtoZ),
            (JosephsonCircuits.YtoA,JosephsonCircuits.AtoY),
            (JosephsonCircuits.YtoB,JosephsonCircuits.BtoY),
            (JosephsonCircuits.ZtoB,JosephsonCircuits.BtoZ),
            (JosephsonCircuits.ZtoY,JosephsonCircuits.YtoZ),
        ]
        # single matrix input
        for arg1 in [rand(Complex{Float64},2,2), (StaticArrays.@MMatrix rand(Complex{Float64},2,2))]
            f[1](arg1)
            f[2](arg1)
        end
        # array input
        for arg1 in [rand(Complex{Float64},2,2,10)]
            f[1](arg1)
            f[2](arg1)
        end
        # vector of matrices
        for arg1 in [
                [rand(Complex{Float64},2,2) for i in 1:10],
            ]
            [f[1](arg1[i]) for i in 1:10]
            [f[2](arg1[i]) for i in 1:10]
        end
    end


    #  network devices, consistency check
    x1 = rand(Complex{Float64})
    x2 = rand(Complex{Float64})
    x3 = rand(Complex{Float64})
    x4 = rand(Complex{Float64})
    JosephsonCircuits.ABCD_seriesZ(x1)
    JosephsonCircuits.YtoA(JosephsonCircuits.Y_seriesY(1/x1))

    JosephsonCircuits.ABCD_shuntY(1/x1)
    JosephsonCircuits.ZtoA(JosephsonCircuits.Z_shuntZ(x1))

    JosephsonCircuits.ABCD_tline(x1,x2)
    JosephsonCircuits.ZtoA(JosephsonCircuits.Z_tline(x1,x2))

    JosephsonCircuits.ABCD_PiY(x1,x2,x3)
    JosephsonCircuits.YtoA(JosephsonCircuits.Y_PiY(x1,x2,x3))

    JosephsonCircuits.ABCD_TZ(x1,x2,x3)
    JosephsonCircuits.ZtoA(JosephsonCircuits.Z_TZ(x1,x2,x3))

    JosephsonCircuits.ABCD_coupled_tline(x1,x2,x3,x4)
    JosephsonCircuits.ZtoA(JosephsonCircuits.Z_coupled_tline(x1,x2,x3,x4))


    return true
end

function warmupconnect()
    # define an open
    Sopen = ones(Complex{Float64},1,1)

    # and a short
    Sshort = -ones(Complex{Float64},1,1)

    # and a match
    Smatch = zeros(Complex{Float64},1,1)

    # a splitter
    Ssplitter = Complex{Float64}[-1/3 2/3 2/3;2/3 -1/3 2/3;2/3 2/3 -1/3]

    S1 = rand(Complex{Float64},3,3)
    S2 = rand(Complex{Float64},2,2)

    # with symbols
    networks = [(:S1,S1),(:S2,S2),(:S3,Ssplitter),(:S4,Sopen)]
    connections = [(:S1,:S1,1,2),(:S1,:S2,3,1),(:S3,:S2,2,2),(:S3,:S4,3,1)]
    networkdata, ports = JosephsonCircuits.connectS(networks,connections)
    Sout1 = networkdata[1]

    S = JosephsonCircuits.connectS(Ssplitter,Sopen,3,1)
    S = JosephsonCircuits.connectS(S,S2,2,2)
    S = JosephsonCircuits.connectS(S1,S,3,2)
    S = JosephsonCircuits.connectS(S,1,2)

    # many frequencies
    N = 100

    # define an open
    Sopen = ones(Complex{Float64},1,1,N)

    # and a short
    Sshort = -ones(Complex{Float64},1,1,N)

    # and a match
    Smatch = zeros(Complex{Float64},1,1,N)

    # a splitter
    Ssplitter = zeros(Complex{Float64},3,3,N)
    for i in 1:N
        Ssplitter[:,:,i] .= Complex{Float64}[-1/3 2/3 2/3;2/3 -1/3 2/3;2/3 2/3 -1/3]
    end

    S1 = rand(Complex{Float64},3,3,N)
    S2 = rand(Complex{Float64},2,2,N)

    networks = [("S1",S1),("S2",S2),("S3",Ssplitter),("S4",Sopen)]
    connections = [("S1","S1",1,2),("S1","S2",3,1),("S3","S2",2,2),("S3","S4",3,1)]
    JosephsonCircuits.connectS(networks,connections)
    JosephsonCircuits.solveS(networks,connections)

    networks = [("S1",S1),("S2",S2),("S4",Sopen)]
    connections = [[("S1",1),("S1",2)],[("S1",3),("S2",2),("S4",1)]]
    JosephsonCircuits.connectS(networks,connections)
    JosephsonCircuits.solveS(networks,connections)

    S = JosephsonCircuits.connectS(Ssplitter,Sopen,3,1)
    S = JosephsonCircuits.connectS(S,S2,2,2)
    S = JosephsonCircuits.connectS(S1,S,3,2)
    S = JosephsonCircuits.connectS(S,1,2)


    return true
end

export @syms, hbsolve, hbnlsolve, hblinsolve, parsecircuit, parsesortcircuit,
    calccircuitgraph, symbolicmatrices, numericmatrices, LjtoIc, IctoLj,
    @variables, @register_symbolic, Num, Symbolics, connectS, solveS

# Export the new functions
export NonlinearElement

export identify_nonlinear_elements

export debug_log, get_debug_log, clear_debug_log

export warning_log, get_warning_log, clear_warning_log

# the below precompile directives are to help the compiler perform type inference
# during the precompilation stage (when the package is installed) instead of
# when it is loaded or the functions are run. this is a helpful guide:
# https://timholy.github.io/SnoopCompile.jl/stable/tutorials/invalidations/#Tutorial-on-@snoop_invalidations
# https://timholy.github.io/SnoopCompile.jl/stable/tutorials/snoop_inference/#Tutorial-on-@snoop_inference
# and the basic commands to look at the inference triggers
# julia> using SnoopCompileCore, JosephsonCircuits
# julia> tinf = @snoop_inference JosephsonCircuits.warmupconnect();
# julia> using SnoopCompile, AbstractTrees
# julia> print_tree(tinf, maxdepth=100)

#=
PrecompileTools.@compile_workload begin
    warmup()
    # warmupsyms()
    # warmupsymsold()
    # warmupsymsnew()
    warmupnetwork()
    warmupconnect()
end
=#

#end module
end