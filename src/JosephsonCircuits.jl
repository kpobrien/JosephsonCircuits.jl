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
adjoint-method sensitivities with respect to component values or, through
[`designsensitivities`](@ref), to the design parameters of a circuit
builder.

The same compiled circuit, on the same node flux unknowns, can also be
integrated directly in physical time with [`transientsolve`](@ref), for
pulsed drives and for drives with more tones than a harmonic grid can
hold, with the exact tangent and adjoint of the recorded time steps.

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
import FunctionWrappers: FunctionWrapper

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
# The directories are the chapters of the analysis: a circuit is written,
# parsed and compiled to flat tables and matrices (circuit/); a harmonic
# balance system is assembled from them (harmonics/); it is solved for the
# pump (solvers/), then linearized and swept for the signals and their
# noise, quantum efficiency and sensitivities (linearized/); and the
# scattering parameters meet the network library (networks/) and the SPICE
# tools (spice/). Within a chapter the files are in reading order where
# the dependencies allow it; where a file is placed for a reason other
# than flow order, the reason is noted next to it.

# --- circuit/: writing a circuit and turning it into matrices -----------
# How a component value is written (a number, a symbol, a parameterized
# expression, or a callable of frequency) and how it becomes a number.
include("circuit/values.jl")
# The component models: lumped elements, ports, nonlinear inductors, and
# multiport scattering and Gaussian channel blocks with the matrix
# providers their frequency dependent data comes from.
include("circuit/components.jl")
include("circuit/vectorfit.jl")
# The typed `Circuit` the user writes, the parse of one hierarchy level
# and the node naming and sorting helpers `compile` uses.
include("circuit/parse.jl")
# Flattening the hierarchy (`elaborate`) and lowering it to the integer
# indexed tables the matrix builders read (`compile`).
include("circuit/compile.jl")
# Stamps of multiport scattering blocks into the harmonic balance system:
# a linearized/ concern, included here because `compile` and the legacy
# adapter need its block types.
include("linearized/scatteringblocks.jl")
# The legacy tuple netlist, adapted into a `Circuit`.
include("circuit/legacy.jl")
include("circuit/graph.jl")      # incidence matrix, spanning tree, loops
include("circuit/matrices.jl")   # capacitance and inverse inductance matrices
include("harmonics/sparse.jl")   # sparse matrix helpers shared by the solvers
# The methods, preconditioners and factorizations a caller composes: a
# solvers/ concern, included here because binding stores a method.
include("solvers/options.jl")
include("circuit/bind.jl")       # binding values and pattern-fixed assembly
include("circuit/mna.jl")        # the modified nodal analysis augmentation

# --- harmonics/: the pieces a harmonic balance system is assembled from -
include("harmonics/layout.jl")   # the equivalent real representation and the canonical state
include("harmonics/directcurrent.jl") # everything about direct current, gauge to block
include("harmonics/complexjacobian.jl") # the holomorphic (complex) Jacobian
include("harmonics/pattern.jl")  # sparsity patterns of device stamps
include("harmonics/assembly.jl") # assembly of the real Jacobian structure
include("harmonics/nonlinearterm.jl") # the Josephson nonlinearity, forward map
include("harmonics/nonlineartermtranspose.jl") # ...and its transpose for adjoints
include("harmonics/frequencies.jl") # frequency grids and transform plans
# The system itself: residual, products, assembled Jacobians.
include("harmonics/system.jl")
# The linearized system: the operating point assembled into one matrix per
# signal frequency.
include("linearized/system.jl")

# --- solvers/: the linear algebra the system is handed to ---------------
include("solvers/solverinfo.jl") # the per stage records and stall diagnostics
include("solvers/factorizations.jl") # sparse factorizations, their cache and solves
include("solvers/linesearch.jl") # the backtracking line search both loops share
include("solvers/newton.jl")     # Newton and quasi-Newton with Anderson acceleration
include("solvers/preconditioners.jl") # the preconditioner interface and the solve record
include("solvers/gmres.jl")      # GMRES and the linear solver objects
include("solvers/newtonkrylov.jl") # the Newton-Krylov driver
include("solvers/floquetdeflation.jl") # residual-image A-DEF1 with physical candidates
include("solvers/cudss.jl")      # the cuDSS factorization type (host stubs)
include("solvers/modecoupling.jl") # the mode coupling preconditioner family
include("solvers/blockclusters.jl") # dense node blocks over the circuit graph, per cluster
include("linearized/blockfactorization.jl") # the same blocks batched over a sweep
# The canonical operators and the preconditioner wrapper of the direct
# current block; it must follow the abstraction it implements in
# solvers/preconditioners.jl.
include("harmonics/canonical.jl")

# --- the harmonic balance solves ----------------------------------------
# The entry point `hbsolve`, the result types both solves return, and the
# docstring fragments shared between the three solver docstrings.
include("linearized/outputs.jl") # scattering parameters, noise and quantum efficiency, what the sweep computes at each frequency
include("linearized/hbsolve.jl")
include("solvers/hbnlsolve.jl")
include("linearized/hblinsolve.jl")
# The device sweep dispatches on the linearized solve's own array types, so
# it follows hblinsolve.jl; its scattering block evaluation and noise
# reductions come first because the sweep drives them.
include("linearized/devicescattering.jl")
include("linearized/devicenoise.jl")
include("linearized/devicesweep.jl")
include("solvers/staged.jl")     # source continuation on a growing harmonic grid
include("solvers/cache.jl")      # reusable workspace for repeated solves
include("solvers/problem.jl")    # the system exposed to external solvers

# --- linearized/: sensitivities and outputs -----------------------------
include("linearized/operatingpoint.jl") # the operating point and its implicit differentiation
include("linearized/sensitivities.jl") # the fixed point stamps and the contraction
include("linearized/designsensitivities.jl")
include("linearized/keyed.jl")   # keyed array output helpers

# --- transient/: the circuit integrated in time -------------------------
# The same compiled circuit, on the same node flux state and matrices as
# harmonic balance at one mode, stepped in physical time for pulsed and
# many tone drives; the Jacobian of a step is the real Jacobian of
# harmonic balance with the step's linear term folded in.
include("transient/system.jl")
include("transient/solve.jl")
include("transient/gauss.jl")
include("transient/batch.jl")
include("transient/sensitivity.jl")
include("transient/iq.jl")       # windowed I/Q of a port trace and its transpose
include("transient/quantum.jl")  # photon normalized temporal modes of a port trace
include("transient/noise.jl")    # the physical baths, and the noise

# --- networks/: the network library -------------------------------------
include("networks/parameters.jl") # S, Z, Y, ABCD, ... conversions
include("networks/networks.jl")  # closed form networks (lines, couplers, ...)
include("networks/connections.jl") # connecting scattering parameter networks
include("networks/quantumoptics.jl") # symplectic and Bogoliubov utilities
# Phase unwrapping, copied from DSP.jl (see the license header in the file).
include("networks/unwrap.jl")

# --- spice/: exporting SPICE netlists and running them ------------------
include("spice/export.jl")
include("spice/utils.jl")
include("spice/wrapper.jl")
include("spice/raw.jl")
include("spice/transient.jl") # the transient run through WRspice

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
# requires updating test/JosephsonCircuits.jl, which holds the fixtures
# of the same circuit that are not part of the workload.

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
# closed operator set (see circuit/values.jl) would make a confusing public
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
    Floquet, Always, Probe, Never, KLUfactorization, LUfactorization,
    QRfactorization, CUDSSFactorization, BlockFactorization

# the typed circuit representation and its component models
export Circuit, Interface, Instance, Ground, Net, PortRef, PinRef,
    Inductor, Capacitor, Resistor, CurrentSource, VoltageSource, Port,
    MutualInductor, JosephsonJunction, NonlinearInductor, PolynomialCPR,
    ScatteringParameters, GaussianChannel, TransmissionLine, RationalScattering, Passive, Lossless,
    ScatteringLimit, OpenDC, ShortDC, ThroughDC, ScatteringDC,
    ThermalEquilibrium, NoiseCovariance, ConjugateSymmetry, Native,
    elaborate, ElaboratedCircuit, quadraturetransform,
    ComponentNotSupportedError

# the circuit integrated in time
export TransientSource, transientproblem, transientstate, transientsolve,
    transientdemodulate, transienttangent, transientadjoint, transientinjection,
    Trapezoidal, GaussLegendre, BackwardEuler, WRspice, TransientReuse, TransientBatchSolution,
    transientiqplan, transientiq!, transientiq, transientiqvjp!,
    transientquantumplan, transientquantum, transientquantum!, transientquantumvjp!,
    transientnoisebaths, transientnoise, transientgain, transientquantumdiagnostics,
    transientquantumefficiency


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

# The circuit in time on the same amplifier: the Gauss-Legendre and the
# trapezoidal steps, the records the responses read, the tangent, the
# adjoint and the noise of a short record, a batch of two conditions,
# which are the paths a first transient pays for otherwise.
function warmuptransient()
    circuit = warmupcircuit(50.0, 100.0e-15, 1000.0e-12, 1000.0e-15)
    fp, ip = 4.75e9, 0.00565e-6
    pump(t) = t <= 0 ? 0.0 : 2ip*cospi(2fp*t)
    base = transientproblem(circuit; sources = [TransientSource(1, t -> 0.0)])
    problem = transientproblem(base; sources = [TransientSource(1, pump)])
    n, T = 64, 2e-9
    solution = transientsolve(problem, (0.0, T*(n - 1)/n); dt = T/n, method = GaussLegendre(), record = :phases)
    transientsolve(problem, (0.0, T*(n - 1)/n); dt = T/n)
    checkpointed = transientsolve(problem, (0.0, T*(n - 1)/n); dt = T/n, method = GaussLegendre(), record = :checkpoints)
    currents = [1e-8*cospi(2*4.7e9*t) for _ in 1:1, t in solution.times]
    weights = [cospi(2*4.7e9*t) for _ in 1:1, t in solution.times]
    transienttangent(solution, currents)
    transientadjoint(solution, weights)
    transienttangent(checkpointed, currents)
    transientadjoint(checkpointed, weights)
    # the measured tone on a bin of the record, and the bath on two bins
    plan = transientquantumplan(solution.times, [9/T])
    transientnoise(solution, plan; frequencies = [8/T, 9/T], weights = [1/T, 1/T], inputs = plan)
    batch = transientsolve([problem, transientproblem(base; sources = [TransientSource(1, t -> pump(t)/2)])],
        (0.0, T*(n - 1)/n); dt = T/n, record = :phases)
    transienttangent(batch, currents)
    transientadjoint(batch, weights)
    # the batch's noise on both methods, of a member and of a range of
    # members, which are views of the batch, and of a checkpointed batch;
    # the pulsed gain of the solution and of the batch
    transientnoise(batch, plan; frequencies = [8/T, 9/T], weights = [1/T, 1/T], inputs = plan)
    transientnoise(batch, plan; frequencies = [8/T, 9/T], weights = [1/T, 1/T], inputs = plan, method = :forward)
    transientnoise(batch[1], plan; frequencies = [8/T, 9/T], weights = [1/T, 1/T], inputs = plan)
    transientnoise(batch[1:2], plan; frequencies = [8/T, 9/T], weights = [1/T, 1/T], inputs = plan)
    transientnoise(transientsolve([problem, problem], (0.0, T*(n - 1)/n); dt = T/n, record = :checkpoints),
        plan; frequencies = [8/T, 9/T], weights = [1/T, 1/T], inputs = plan)
    transientgain(solution, plan, plan)
    transientgain(batch, plan, plan)
    transienttangent(batch[1], currents)
    # several directions at once, and directions given at the stages
    transienttangent(solution, cat(currents, currents; dims = 3))
    transienttangent(solution, [currents[1, i] for _ in 1:1, _ in 1:3, i in 1:n, _ in 1:1])
    # a windowed measurement with an envelope
    half = solution.times[1:n÷2]
    windowed = transientquantumplan(half, [9/T]; envelopes = reshape(sinpi.((half .- half[1]) ./ (T/2)) .^ 2, :, 1))
    transientnoise(solution, windowed; frequencies = [8/T, 9/T], weights = [1/T, 1/T])
    transientgain(solution, windowed, windowed)
    # a line and a lossy rational block ahead of the junction: the line
    # history, the block states and their noise, on both noise methods
    a = 2pi*4e9
    block = RationalScattering(-a .* Matrix(1.0I, 2, 2), a .* Matrix(1.0I, 2, 2), 0.8 .* [0.0 1.0; 1.0 0.0],
        zeros(2, 2); zref = 50.0)
    front = Circuit([("p1", "1", "0", Port(1; Z0 = 50.0)), ("line", "1", "2", TransmissionLine(50.0, 0.02)),
        ("b", "2", "3", block), ("cc", "3", "4", Capacitor(100.0e-15)),
        ("jj", "4", "0", JosephsonJunction(1000.0e-12)), ("cj", "4", "0", Capacitor(1000.0e-15))])
    fronted = transientsolve(transientproblem(front; sources = [TransientSource(1, pump)]),
        (0.0, T*(n - 1)/n); dt = T/n, method = GaussLegendre(), record = :phases)
    transienttangent(fronted, currents)
    transientadjoint(fronted, weights)
    transientnoise(fronted, plan; frequencies = [8/T, 9/T], weights = [1/T, 1/T], inputs = plan)
    transientnoise(fronted, plan; frequencies = [8/T, 9/T], weights = [1/T, 1/T], inputs = plan, method = :forward)
    return nothing
end

PrecompileTools.@compile_workload begin
    warmup()
    warmupsyms()
    warmuptransient()
    # `warmupnetwork()` is deliberately not part of the workload. It
    # compiles every network parameter conversion for every input shape,
    # which is a large fraction of the total precompile time, while a cold
    # first call of any one conversion is cheap.
    warmupconnect()
end

end # module JosephsonCircuits
