__precompile__(true)

module JosephsonCircuits

import Graphs
import FFTW
import Statistics
import NLsolve
import KLU
import UUIDs
import SymbolicUtils: @syms, substitute, Symbolic
import OrderedCollections: OrderedDict

using LinearAlgebra
using SparseArrays
using Printf
using Requires
using BandedMatrices

# define the zero for symbolic numbers so that we can view the sparse arrays
Base.zero(::Type{Symbolic{Number}}) = 0
Base.zero(::Type{Any}) = 0
Base.zero(::Type{Nothing}) = 0


# Symbolics is needed for calculating the inverse inductance matrix if there
# are mutual inductors. For large systems with many mutual inductors this may
# result in very large equations so is not recommended. 
function __init__()
    @require Symbolics="0c5d862f-8b57-4792-8d23-62f2024744c7" begin
        function calcsymbolicinvLn(L,Lb,Rbn)
            s =  sparse(transpose(Rbn[Lb.nzind,:])*(Symbolics.sym_lu(L)\ Matrix(Rbn[Lb.nzind,:])))
            return SparseMatrixCSC(s.m, s.n, s.colptr, s.rowval,Symbolics.value.(s.nzval))
        end

        """
            checkissymbolic(a::Num)

        Check if a is a symbolic variable, when the Symbolics package is 
        loaded. The default version for SymbolicUtils is defined in matutils.jl.

        # Examples
        ```jldoctest
        julia> import Symbolics;Symbolics.@variables w;JosephsonCircuits.checkissymbolic(w)
        true

        julia> JosephsonCircuits.checkissymbolic(1.0)
        false
        ```
        """
        function checkissymbolic(a::Symbolics.Num)
            return !(Symbolics.value(a) isa Number)
        end

    end
end


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
include("hbsolve.jl")
include("fftutils.jl")
include("qesparams.jl")

# Experimental multi-tone harmonic balance code.
# include("hbsolve2.jl")

# These are for exporting SPICE netlists and running simulations in
# WRSPICE or Xyce. 
include("exportnetlist.jl")
include("spiceutils.jl")
include("spicewrapper.jl")
include("spiceraw.jl")

# This is for reading and writing files in the touchstone format. 
include("touchstone.jl")

# Precompile directives to improve time to first use. 
# include("precompile.jl")


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
    #reduced flux quantum in Weber, H*A
    # hbar/(2*charge of electron)
    # phi0 = 3.29105976e-16
    return phi0./Ic
end

function warmup()

    # define the circuit components
    circuit = Array{Tuple{String,String,String,Union{Complex{Float64}, Symbol,Int64}},1}(undef,0)

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

    return hbsolve(2*pi*(4.5:0.1:5.0)*1e9,2*pi*4.75001*1e9,0.00565e-6,8,8,circuit,circuitdefs,pumpports=[1],solver=:klu);

end

function warmup2()

    # define the circuit components
    circuit = Array{Tuple{String,String,String,Union{Complex{Float64}, Symbol,Int64}},1}(undef,0)

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

    return hbsolve2(2*pi*(4.5:0.1:5.0)*1e9,2*pi*4.75001*1e9,0.00565e-6,8,8,circuit,circuitdefs,pumpports=[1],solver=:klu);

end


function warmupsyms()

    @syms Ipump Rleft Cc Lj Cj

    # define the circuit components
    circuit = Array{Tuple{String,String,String,Any},1}(undef,0)

    # port on the left side
    push!(circuit,("P1","1","0",1))
    push!(circuit,("I1","1","0",Ipump))
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

    return hbsolve(2*pi*(4.5:0.1:5.0)*1e9,2*pi*4.75001*1e9,0.00565e-6,8,8,circuit,circuitdefs,pumpports=[1],solver=:klu);

end

export @syms, hbnlsolve, hblinsolve, hbsolve, hbsolve2, hbnlsolve2, parsecircuit,
    parsesortcircuit, calccircuitgraph, symbolicmatrices, numericmatrices,
    LjtoIc, IctoLj


#end module
end
