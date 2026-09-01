__precompile__(true)

module JosephsonCircuits

import Graphs
import FFTW
import KLU
import KernelAbstractions
import KernelAbstractions: @kernel, @index, @Const, CPU, Backend
import Atomix
import UUIDs
import RuntimeGeneratedFunctions
import AxisKeys
import PrecompileTools
import OrderedCollections
import StaticArrays
import Statistics
import FastInterpolations

using LinearAlgebra
using SparseArrays
using Printf
using Touchstone

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

A constant for the speed of light which is 2.99792458e8 m/s.
"""
const speed_of_light = 2.99792458e8


"""
    const planck_constant

A constant for the Planck constant which is 6.62607015e-34 J*s.
"""
const planck_constant = 6.62607015e-34

"""
    const reduced_planck_constant

A constant for the reduced Planck constant which is planck_constant/(2*pi).
"""
const reduced_planck_constant = planck_constant/(2*pi)

"""
    const boltzmann_constant

A constant for the Boltzmann constant which is 1.380649×10−23 J/K.
"""
const boltzmann_constant = 1.380649e-23


# The file structure below mimics the typical analysis flow. An input file
# is first parsed, then the incidence matrix is calculated, then the 
# capacitance, inverse inductance, and other matrices are calculated, then
# then circuit is solved using the harmonic balance method, and postprocessed
# to determine the scattering parameters and quantum efficiency. 
# How a component value is written and how it becomes a number.
# --- the circuit -------------------------------------------------------
# How a component value is written and how it becomes a number.
include("circuitvalue.jl")
# What each component is: the models, their noise, and the matrix providers
# a scattering block reads its parameters from.
include("circuitmodel.jl")
# The input path: the Circuit a user writes and how one level of it is
# parsed, then the hierarchy flattened and lowered to the compiled tables.
include("parseinput.jl")
include("circuitcompile.jl")
# Scattering parameter stamps for the solvers, and the legacy tuple netlist
# adapted into a Circuit.
include("scatteringstamp.jl")
include("legacyadapter.jl")

# --- from a compiled circuit to matrices --------------------------------
include("graphproc.jl")
include("capindmat.jl")
include("matutils.jl")
include("circuitbind.jl")
include("mna.jl")
include("dcconductance.jl")

# --- the pieces a harmonic balance system is assembled from -------------
include("realcomplexconv.jl")
include("complexjacobian.jl")
include("devicepattern.jl")
include("structureassembly.jl")
include("nonlinearterm.jl")
include("nonlineartermtranspose.jl")
include("fftutils.jl")

# --- the system itself --------------------------------------------------
include("hbsystem.jl")

# --- the linear algebra it is handed to ---------------------------------
include("nlsolve.jl")
include("krylov.jl")
include("batchedblocks.jl")
include("cudss.jl")
include("modepreconditioner.jl")

# --- the harmonic balance solves ----------------------------------------
# the canonical state layout the direct current block is carried in; it
# supplies a preconditioner, so it follows the abstraction it implements
include("compositelayout.jl")
# the entry point, the result types both solves return, and the shared
# pieces of their docstrings, which are interpolated at definition time
include("hbsolve.jl")
include("hbnlsolve.jl")
include("hblinsolve.jl")
# the device sweep dispatches on the linearized solve's own array types, so
# it follows it
include("devicelinsolve.jl")
include("stagedsolve.jl")
include("hbcache.jl")
include("hbnonlinearproblem.jl")

# --- sensitivities ------------------------------------------------------
include("sensitivities.jl")
include("designsensitivities.jl")

# --- turning a solution into what was asked for -------------------------
include("networkparamconversion.jl")
include("networks.jl")
include("networkconnection.jl")
include("qesparams.jl")
include("keyedarrayutils.jl")
include("quantumoptics.jl")
# From DSP.jl https://github.com/JuliaDSP/DSP.jl/blob/master/src/unwrap.jl
include("unwrap.jl")

# --- exporting SPICE netlists and running them --------------------------
include("exportnetlist.jl")
include("spiceutils.jl")
include("spicewrapper.jl")
include("spiceraw.jl")

# These are deprecated functions
include("deprecated.jl")

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

# The circuit the warmups share: a Josephson parametric amplifier whose port
# carries its own reference impedance. Written in the typed format because
# that is the input path worth precompiling; the legacy tuple netlist is
# adapted into a `Circuit` anyway and is not expected to see much use.
function warmupcircuit(Rleft, Cc, Lj, Cj)
    return Circuit(
        ["P1" => Port(1; Z0 = Rleft), "C1" => Capacitor(Cc),
         "Lj1" => JosephsonJunction(Lj), "C2" => Capacitor(Cj)],
        [Net("1", [("P1",1), ("C1",1)]),
         Net("2", [("C1",2), ("Lj1",1), ("C2",1)]),
         Net("0", [("P1",2), ("Lj1",2), ("C2",2), Ground])])
end

# the values the warmups solve at, as a symbol keyed dictionary
warmupdefs(Rleft, Cc, Lj, Cj) = Dict(
    Lj => 1000.0e-12,
    Cc => 100.0e-15,
    Cj => 1000.0e-15,
    Rleft => 50.0,
)

function warmup()

    circuit = warmupcircuit(:Rleft, :Cc, :Lj, :Cj)
    circuitdefs = Dict{Symbol,Complex{Float64}}(
        :Lj => 1000.0e-12,
        :Cc => 100.0e-15,
        :Cj => 1000.0e-15,
        :Rleft => 50.0,
    )

    ws = 2*pi*(4.5:0.5:5.0)*1e9
    wp = (2*pi*4.75001*1e9,)
    sources = [(mode=(1,),port=1,current=0.00565e-6)]
    Nmodulationharmonics = (2,)
    Npumpharmonics = (4,)

    return hbsolve(ws, wp, sources, Nmodulationharmonics,
        Npumpharmonics, circuit, circuitdefs;ftol=1e-12)
end


function warmupsyms()

    @params R Cc Lj Cj
    circuit = warmupcircuit(R, Cc, Lj, Cj)
    circuitdefs = warmupdefs(R, Cc, Lj, Cj)

    ws = 2*pi*(4.5:0.5:5.0)*1e9
    wp = (2*pi*4.75001*1e9,)
    sources = [(mode=(1,),port=1,current=0.00565e-6)]
    Nmodulationharmonics = (2,)
    Npumpharmonics = (4,)

    return hbsolve(ws, wp, sources, Nmodulationharmonics,
        Npumpharmonics, circuit, circuitdefs;ftol=1e-12)
end

function warmupcompile()

    @params Rleft Cc Lj Cj
    return compile(warmupcircuit(Rleft, Cc, Lj, Cj))
end

function warmupnumericmatrices()

    @params Rleft Cc Lj Cj
    circuit = warmupcircuit(Rleft, Cc, Lj, Cj)
    return numericmatrices(circuit, warmupdefs(Rleft, Cc, Lj, Cj))
end

function warmuphblinsolve()

    @params Rleft Cc Lj Cj
    circuit = warmupcircuit(Rleft, Cc, Lj, Cj)
    return hblinsolve(2*pi*(4.5:0.1:5.0)*1e9, circuit,
        warmupdefs(Rleft, Cc, Lj, Cj))
end

function warmupvvn()

    @params Rleft Cc Lj Cj
    psc = compile(warmupcircuit(Rleft, Cc, Lj, Cj); sorting = :number)

    return componentvaluestonumber(psc.componentvalues,
        warmupdefs(Rleft, Cc, Lj, Cj))
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

    S = JosephsonCircuits.interconnectS(Ssplitter,Sopen,3,1)
    S = JosephsonCircuits.interconnectS(S,S2,2,2)
    S = JosephsonCircuits.interconnectS(S1,S,3,2)
    S = JosephsonCircuits.intraconnectS(S,1,2)

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

    S = JosephsonCircuits.interconnectS(Ssplitter,Sopen,3,1)
    S = JosephsonCircuits.interconnectS(S,S2,2,2)
    S = JosephsonCircuits.interconnectS(S1,S,3,2)
    S = JosephsonCircuits.intraconnectS(S,1,2)

    return true
end

export hbsolve, hbnlsolve, hblinsolve, compile,
    calccircuitgraph, symbolicmatrices, numericmatrices, LjtoIc, IctoLj,
    connectS, solveS

# The CircuitValues expression type is INTERNAL: it is the lowering target
# of the Symbolics extension and the representation of parameterized
# netlist-file expressions. Users parameterize circuits with ordinary
# Julia functions (builders) and numbers; there is deliberately no public
# symbolic-looking parameter type, whose closed operator set would be a
# confusing boundary. `@params` is used by internal warmup code only.
import .CircuitValues: @params
export FrequencyDependent, designsensitivities, designjacobian,
    hbcache, hbsolve!,
    hbnonlinearproblem, JacobianOperator, preconditioner, hbresidual!,
    hbjvp!, hbvjp!, hbjacobian!, hbd2F!, hbd3F!, hbdFdp!, jacobianprototype,
    setdrive!, drivenresidual!, NewtonKrylov, Newton, QuasiNewton,
    ExternalSolver, InternalGMRES, KrylovJL

# the typed circuit representation
export Circuit, Interface, Instance, Ground, Net, PortRef, PinRef,
    Inductor, Capacitor, Resistor, CurrentSource, VoltageSource, Port,
    MutualInductor, JosephsonJunction, NonlinearInductor, PolynomialCPR,
    ScatteringParameters, GaussianChannel, TransmissionLine, Passive, Lossless,
    ScatteringLimit, OpenDC, ShortDC, ThroughDC, ScatteringDC,
    ThermalEquilibrium, NoiseCovariance, ConjugateSymmetry, Native,
    elaborate, ElaboratedCircuit, quadraturetransform,
    ComponentNotSupportedError


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

PrecompileTools.@compile_workload begin
    warmup()
    warmupsyms()
    # warmupnetwork() is deliberately NOT part of the workload: it costs
    # about 20 seconds of precompile time (a third of the total) compiling
    # every network-parameter conversion for every input shape, while a
    # cold first call of one conversion costs a quarter of a second. The
    # function itself is kept as a manual probe and for the test suite.
    warmupconnect()
end

#end module
end
