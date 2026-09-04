__precompile__(true)

"""
    JosephsonCircuits

A frequency domain simulator for superconducting circuits containing
Josephson junctions, capacitors, inductors, mutual inductors, resistors,
and multiport scattering parameter blocks.

Circuits are solved by harmonic balance in a modified nodal analysis
formulation in the node flux basis. A strong periodic drive (the pump) is
solved with [`hbnlsolve`](@ref), the circuit is linearized about that
operating point and swept over weak signal frequencies with
[`hblinsolve`](@ref), and [`hbsolve`](@ref) runs the two in sequence. From
the linearized solution the package computes scattering parameters, noise
scattering parameters, quantum efficiency, commutation relations, and
adjoint-method sensitivities with respect to component values.

A circuit is written as a [`Circuit`](@ref) of typed component models, or
as a legacy netlist of `(name, node1, node2, value)` tuples. The stages a
circuit passes through, and the files that implement them, are listed
next to the `include` statements below.
"""
module JosephsonCircuits

import Graphs
import FFTW
import KLU
# the orderings of `kluordered` come from the CHOLMOD library that ships
# with SparseArrays
import SparseArrays.CHOLMOD
import SparseArrays.LibSuiteSparse
import KernelAbstractions
import KernelAbstractions: @kernel, @index, @Const, CPU, Backend
import Atomix
import UUIDs
import AxisKeys
import PrecompileTools
import OrderedCollections
import StaticArrays
import Statistics
import FastInterpolations

using LinearAlgebra
using SparseArrays
using Touchstone

# === physical constants ===

"""
    const phi0

The reduced magnetic flux quantum `hbar/(2e)` in Weber (equivalently
H*A). This is the flux scale that relates a Josephson junction's inductance
to its critical current, `Ic = phi0/Lj`, and the unit in which the branch
phases of the harmonic balance solution are measured.
"""
const phi0 = 3.29105976e-16

"""
    const Phi0

The magnetic flux quantum `h/(2e)` in Weber (equivalently H*A), equal to
`2*pi*phi0`.
"""
const Phi0 = 2.067833848e-15

"""
    const speed_of_light

The speed of light in vacuum, 2.99792458e8 m/s, the default phase velocity
of a [`TransmissionLine`](@ref).
"""
const speed_of_light = 2.99792458e8


"""
    const planck_constant

The Planck constant `h`, 6.62607015e-34 J*s.
"""
const planck_constant = 6.62607015e-34

"""
    const reduced_planck_constant

The reduced Planck constant `hbar = h/(2*pi)` in J*s.
"""
const reduced_planck_constant = planck_constant/(2*pi)

"""
    const boltzmann_constant

The Boltzmann constant `k_B`, 1.380649e-23 J/K, used with
[`reduced_planck_constant`](@ref) to convert a physical temperature into
the thermal occupation of a noise channel.
"""
const boltzmann_constant = 1.380649e-23


# === source files, in the order a circuit passes through them ===
#
# Files are included roughly in the order of the analysis flow: a circuit is
# written, parsed and compiled to flat tables; the tables are turned into
# the incidence, capacitance and inverse inductance matrices; a harmonic
# balance system is assembled from those and solved, first for the pump
# and then linearized for the signals; and the solution is post-processed
# into scattering parameters, noise, quantum efficiency, sensitivities and
# exported netlists. Where a file is placed for a reason other than flow
# order, the reason is noted next to it.

# --- writing a circuit --------------------------------------------------
# How a component value is written (a number, a symbol, a parameterized
# expression, or a callable of frequency) and how it becomes a number.
include("circuitvalue.jl")
# The component models: lumped elements, ports, nonlinear inductors, and
# multiport scattering and Gaussian channel blocks with the matrix
# providers their frequency dependent data comes from.
include("circuitmodel.jl")
# The typed `Circuit` the user writes, the parse of one hierarchy level
# and the node naming and sorting helpers `compile` uses.
include("parseinput.jl")
# Flattening the hierarchy (`elaborate`) and lowering it to the integer
# indexed tables the matrix builders read (`compile`).
include("circuitcompile.jl")
# Stamps of multiport scattering blocks into the harmonic balance system.
include("scatteringstamp.jl")
# The legacy tuple netlist, adapted into a `Circuit`.
include("legacyadapter.jl")

# --- from a compiled circuit to matrices --------------------------------
include("graphproc.jl")        # incidence matrix, spanning tree, loops
include("capindmat.jl")        # capacitance and inverse inductance matrices
include("matutils.jl")         # sparse matrix helpers shared by the solvers
include("solveroptions.jl")    # the methods, preconditioners and factorizations a caller composes
include("circuitbind.jl")      # binding values and pattern-fixed assembly
include("mna.jl")              # the modified nodal analysis augmentation
include("dcconductance.jl")    # the explicit direct current block

# --- pieces a harmonic balance system is assembled from -----------------
include("realcomplexconv.jl")  # the equivalent real representation
include("complexjacobian.jl")  # the holomorphic (complex) Jacobian
include("devicepattern.jl")    # sparsity patterns of device stamps
include("structureassembly.jl")# assembly of the real Jacobian structure
include("nonlinearterm.jl")    # the Josephson nonlinearity, forward map
include("nonlineartermtranspose.jl") # ...and its transpose for adjoints
include("fftutils.jl")         # frequency grids and transform plans

# --- the system itself --------------------------------------------------
include("hbsystem.jl")

# --- the linear algebra it is handed to ---------------------------------
include("nlsolve.jl")          # Newton and quasi-Newton with factorizations
include("krylov.jl")           # Newton-Krylov and its preconditioner types
include("floquetdeflation.jl") # residual-image A-DEF1 with physical candidates
include("batchedblocks.jl")    # batched small block factorizations
include("cudss.jl")            # the cuDSS factorization type (host stubs)
include("modepreconditioner.jl")
include("blockfactorization.jl") # dense node blocks over the circuit graph

# --- the harmonic balance solves ----------------------------------------
# The canonical state layout the direct current block is carried in. It
# supplies a preconditioner, so it must follow the abstraction it implements
# in krylov.jl.
include("compositelayout.jl")
# The entry point `hbsolve`, the result types both solves return, and the
# docstring fragments shared between the three solver docstrings.
include("hbsolve.jl")
include("hbnlsolve.jl")
include("hblinsolve.jl")
# The device sweep dispatches on the linearized solve's own array types, so
# it follows hblinsolve.jl.
include("devicelinsolve.jl")
include("stagedsolve.jl")      # source continuation on a growing harmonic grid
include("hbcache.jl")          # reusable workspace for repeated solves
include("hbnonlinearproblem.jl") # the system exposed to external solvers

# --- sensitivities ------------------------------------------------------
include("sensitivities.jl")
include("designsensitivities.jl")

# --- turning a solution into what was asked for -------------------------
include("networkparamconversion.jl") # S, Z, Y, ABCD, ... conversions
include("networks.jl")         # closed form networks (lines, couplers, ...)
include("networkconnection.jl") # connecting scattering parameter networks
include("qesparams.jl")        # scattering parameters, noise and quantum efficiency
include("keyedarrayutils.jl")  # keyed array output helpers
include("quantumoptics.jl")    # symplectic and Bogoliubov utilities
# Phase unwrapping, copied from DSP.jl (see the license header in the file).
include("unwrap.jl")

# --- exporting SPICE netlists and running them --------------------------
include("exportnetlist.jl")
include("spiceutils.jl")
include("spicewrapper.jl")
include("spiceraw.jl")

# Deprecated entry points, kept so that older scripts keep running with a
# warning.
include("deprecated.jl")

# Helpers the test suite uses to print and compare solver output.
include("testutils.jl")


"""
    LjtoIc(Lj)

The critical current `Ic = phi0/Lj` in Amperes of a Josephson junction
with junction inductance `Lj` in Henries.

# Examples
```jldoctest
julia> LjtoIc(100e-12)
3.29105976e-6
```
"""
function LjtoIc(Lj)
    return phi0./Lj
end

"""
    IctoLj(Ic)

The junction inductance `Lj = phi0/Ic` in Henries of a Josephson junction
with critical current `Ic` in Amperes.

# Examples
```jldoctest
julia> IctoLj(3.29105976e-6)
1.0e-10
```
"""
function IctoLj(Ic)
    # Lj = phi0/Ic has the same form as Ic = phi0/Lj
    return LjtoIc(Ic)
end

# === precompilation workloads ===
#
# The `warmup*` functions below exercise the main code paths so that
# PrecompileTools can compile them when the package is installed rather
# than on first use. The test suite also calls them and compares their
# output against stored reference values, so changing what they compute
# requires updating test/JosephsonCircuits.jl.

# The circuit every warmup shares: a single junction parametric amplifier,
# capacitively coupled to a port which owns its own matched termination.
# It is written in the typed format because that is the input path worth
# precompiling; a legacy tuple netlist is adapted into a `Circuit` first
# and then takes the same path.
function warmupcircuit(Rleft, Cc, Lj, Cj)
    return Circuit(
        ["P1" => Port(1; Z0 = Rleft), "C1" => Capacitor(Cc),
         "Lj1" => JosephsonJunction(Lj), "C2" => Capacitor(Cj)],
        [Net("1", [("P1",1), ("C1",1)]),
         Net("2", [("C1",2), ("Lj1",1), ("C2",1)]),
         Net("0", [("P1",2), ("Lj1",2), ("C2",2), Ground])])
end

# The component values the warmups solve at, keyed by whatever parameter
# objects `warmupcircuit` was given (symbols or `@params` parameters).
warmupdefs(Rleft, Cc, Lj, Cj) = Dict(
    Lj => 1000.0e-12,
    Cc => 100.0e-15,
    Cj => 1000.0e-15,
    Rleft => 50.0,
)

# The full pump plus signal solve with symbol valued components.
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


# The same solve with `CircuitValue` parameters, which is the form a
# Symbolics `Num` is lowered to.
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

# Elaboration and compilation of the typed circuit alone.
function warmupcompile()

    @params Rleft Cc Lj Cj
    return compile(warmupcircuit(Rleft, Cc, Lj, Cj))
end

# The capacitance and inverse inductance matrices at numeric values.
function warmupnumericmatrices()

    @params Rleft Cc Lj Cj
    circuit = warmupcircuit(Rleft, Cc, Lj, Cj)
    return numericmatrices(circuit, warmupdefs(Rleft, Cc, Lj, Cj))
end

# A linear (no pump) frequency sweep.
function warmuphblinsolve()

    @params Rleft Cc Lj Cj
    circuit = warmupcircuit(Rleft, Cc, Lj, Cj)
    return hblinsolve(2*pi*(4.5:0.1:5.0)*1e9, circuit,
        warmupdefs(Rleft, Cc, Lj, Cj))
end

# Resolving the compiled component values to numbers.
function warmupvvn()

    @params Rleft Cc Lj Cj
    psc = compile(warmupcircuit(Rleft, Cc, Lj, Cj); sorting = :number)

    return componentvaluestonumber(psc.componentvalues,
        warmupdefs(Rleft, Cc, Lj, Cj))
end

# Every network parameter conversion, for every input shape it accepts.
# This is not part of the precompile workload (see the note there) but is
# kept as a manual probe and for the test suite.
function warmupnetwork()
    # conversions to and from scattering parameters, which take a port
    # impedance keyword
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

    # conversions between the other representations, which take no port
    # impedance
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


    # closed form two port networks in their ABCD, Z and Y forms
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

# Connecting scattering parameter networks, with symbol and string names,
# with single matrices and with frequency indexed arrays, through both the
# graph based `connectS` and the linear system based `solveS`.
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

# The `CircuitValues` expression type is internal: it is what a Symbolics
# `Num` and a parameterized netlist file expression are lowered to. It is
# deliberately not exported as a user facing symbolic type, because its
# closed operator set (see circuitvalue.jl) would make a confusing public
# boundary; users parameterize circuits with symbols, numbers, and ordinary
# Julia functions. `@params` is imported for the warmups above only.
import .CircuitValues: @params
export FrequencyDependent, designsensitivities, designjacobian,
    hbcache, hbsolve!,
    hbnonlinearproblem, JacobianOperator, preconditioner, hbresidual!,
    hbjvp!, hbvjp!, hbjacobian!, hbd2F!, hbd3F!, hbdFdp!, jacobianprototype,
    setdrive!, drivenresidual!, NewtonKrylov, Newton, QuasiNewton,
    ExternalSolver, GMRES, KrylovJL, Staged,
    BlockDiagonal, FullJacobian, HarmonicBand, MeasuredBand, Clusters,
    CoupledModes, CouplingMask, Automatic,
    Recycling, Floquet, Always, Probe, Never, KLUfactorization, LUfactorization,
    QRfactorization, CUDSSFactorization, BlockFactorization

# the typed circuit representation and its component models
export Circuit, Interface, Instance, Ground, Net, PortRef, PinRef,
    Inductor, Capacitor, Resistor, CurrentSource, VoltageSource, Port,
    MutualInductor, JosephsonJunction, NonlinearInductor, PolynomialCPR,
    ScatteringParameters, GaussianChannel, TransmissionLine, Passive, Lossless,
    ScatteringLimit, OpenDC, ShortDC, ThroughDC, ScatteringDC,
    ThermalEquilibrium, NoiseCovariance, ConjugateSymmetry, Native,
    elaborate, ElaboratedCircuit, quadraturetransform,
    ComponentNotSupportedError


# The precompile workload runs the warmups when the package is installed so
# that type inference and compilation happen then rather than at load time
# or on first call. To see what still gets inferred at run time:
#
#   julia> using SnoopCompileCore, JosephsonCircuits
#   julia> tinf = @snoop_inference JosephsonCircuits.warmupconnect();
#   julia> using SnoopCompile, AbstractTrees
#   julia> print_tree(tinf, maxdepth=100)
#
# See the SnoopCompile.jl tutorials on inference and invalidations.

PrecompileTools.@compile_workload begin
    warmup()
    warmupsyms()
    # `warmupnetwork()` is deliberately not part of the workload. It
    # compiles every network parameter conversion for every input shape,
    # which is a large fraction of the total precompile time, while a cold
    # first call of any one conversion is cheap.
    warmupconnect()
end

end # module JosephsonCircuits
