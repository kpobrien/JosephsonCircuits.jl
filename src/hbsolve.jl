
"""
    NonlinearHB(nodeflux, Rbnm, Ljb, Lb, Ljbm, Nmodes, Nbranches, S)

A simple structure to hold the nonlinear harmonic balance solutions.

# Fields
- `w`: a tuple containing the the angular frequency of the pump in radians/s.
- `frequencies`:
- `nodeflux`: the node fluxes resulting from inputs at each frequency and
    port.
- `Rbnm`: incidence matrix to convert between the node and branch basis.
- `Ljb`: sparse vector of Josephson junction inductances.
- `Lb`: sparse vector of linear inductances.
- `Ljbm`: sparse vector of linear inductances with each element duplicated
    Nmodes times.
- `Nmodes`: the number of signal and idler frequencies.
- `Nbranches`: the number of branches in the circuit.
- `nodenames`: the vector of unique node name strings.
- `componentnames`: the vector of component name strings
- `portnumbers`: vector of port numbers.
- `portindices`: 
- `modes`: tuple of the pump mode indices where (1,) is the pump in the single
    pump case.
- `S`: the scattering matrix relating inputs and outputs for each combination
    of port and frequency. Its zero frequency entries are identically zero
    and carry no information: the waves are in units of
    sqrt(photons/second), whose normalization `1/sqrt(|w|)` has no limit at
    zero (see [`portwavescale`](@ref)), so there is no direct current wave
    to report. The direct current operating point is in `dcnodevoltage`,
    which is a voltage and not a wave.
- `solverinfo`: diagnostics describing the nonlinear solution process.
    See [`SolverInfo`](@ref).
- `operatingpoint`: the converged operating point and the exact real
    Jacobian there, when the solver is called with
    `returnoperatingpoint = true`, otherwise `nothing`. See
    [`HBOperatingPoint`](@ref).
"""
struct NonlinearHB
    w
    frequencies
    nodeflux
    Rbnm
    Ljb
    Lb
    Ljbm
    Nmodes
    Nbranches
    nodes
    ports
    modes
    S
    solverinfo
    operatingpoint
    # The average node voltage in volts from the explicit direct current
    # block, or `nothing` when the circuit has no such block. Ground is
    # excluded and the values are keyed by node name, as `nodeflux` is, so
    # that node indexed code reads the two the same way.
    #
    # Distinct from `nodeflux`, whose zero mode is the static periodic flux
    # which sets inductor currents and junction phases. A vector of zeros is
    # an answer -- a node shorted to ground sits at zero volts -- and is not
    # the same as `nothing`.
    dcnodevoltage
end

function NonlinearHB(w, frequencies, nodeflux, Rbnm, Ljb, Lb, Ljbm, Nmodes,
    Nbranches, nodes, ports, modes, S, solverinfo, operatingpoint)
    return NonlinearHB(w, frequencies, nodeflux, Rbnm, Ljb, Lb, Ljbm,
        Nmodes, Nbranches, nodes, ports, modes, S, solverinfo,
        operatingpoint, nothing)
end

# backwards compatible constructors without the solver diagnostics and
# without the operating point
function NonlinearHB(w, frequencies, nodeflux, Rbnm, Ljb, Lb, Ljbm, Nmodes,
    Nbranches, nodes, ports, modes, S)
    return NonlinearHB(w, frequencies, nodeflux, Rbnm, Ljb, Lb, Ljbm,
        Nmodes, Nbranches, nodes, ports, modes, S,
        SolverInfo(IterationInfo[], NaN, NaN, true, NaN), nothing)
end

function NonlinearHB(w, frequencies, nodeflux, Rbnm, Ljb, Lb, Ljbm, Nmodes,
    Nbranches, nodes, ports, modes, S, solverinfo)
    return NonlinearHB(w, frequencies, nodeflux, Rbnm, Ljb, Lb, Ljbm,
        Nmodes, Nbranches, nodes, ports, modes, S, solverinfo, nothing)
end

"""
    LinearizedHB(w, modes, S, Snoise, Ssensitivity, QE, QEideal, CM,
        nodeflux, nodefluxadjoint, voltage, voltageadjoint, Nmodes, Nnodes,
        Nbranches, Nports, signalindex)

A simple structure to hold the linearized harmonic balance solutions. An
output which was not requested is an empty array.

# Fields
- `w`: the signal frequencies.
- `modes`: tuple of the signal mode indices where (0,) is the signal.
- `S`: the scattering matrix relating inputs and outputs for each combination
    of port and frequency.
- `Snoise`: the scattering matrix relating inputs at the noise channels of the
    dissipative elements and outputs at the physical ports, for each
    combination of port and frequency. The channels are the dissipative
    lumped components first, then one per port of each dissipative
    [`ScatteringParameters`](@ref). Being a scattering matrix it describes a
    transformation and does not depend on temperature.
- `Cnoise`: the added noise covariance at the output ports,
    `sum_c occupation[c] Snoise[c,i] conj(Snoise[c,j])`, returned only when
    `returnCnoise = true`. With `S` this is the Gaussian channel the circuit
    implements: the map takes an input covariance to `S sigma S' + Cnoise`.
    This is where the temperature of each channel enters, `Snoise` supplying
    the transformation and the occupation the state. At zero temperature its
    diagonal is the noise term in the denominator of the quantum efficiency.
- `Ssensitivity`: the derivative of the scattering matrix with respect to a
    relative (logarithmic) perturbation of each component value in
    `sensitivitynames`, at the fixed pump operating point, or the total
    derivative including the shift of the operating point when
    `sensitivityoperatingpoint = true`.
- `QE`: the quantum efficiency for each combination of port and frequency.
- `QEideal`: the quantum efficiency for an ideal amplifier with the same level
    of gain, for each combination of port and frequency.
- `CM`: the commutation relations (equal to ±1), for each combination of port
    and frequency.
- `nodeflux`: the node fluxes resulting from inputs at each frequency and port.
- `nodefluxadjoint`: the node fluxes resulting from inputs at each frequency
    and port with a time reversed modulation.
- `voltage`: the node voltages resulting from inputs at each frequency and port.
- `voltageadjoint`: the node fluxes resulting from inputs at each frequency
    and port with a time reversed modulation.
- `nodenames`: the vector of unique node strings.
- `nodeindices`:
- `componentnames`:
- `componenttypes`:
- `componentnamedict`:
- `mutualinductorbranchnames`:
- `portnumbers`: vector of port numbers.
- `portindices`:
- `portimpedances`: the reference impedance of each port, ordered by port
    number, which is what the scattering parameters are normalized to
- `noiseportimpedanceindices`:
- `sensitivitynames`:
- `sensitivityindices`:
- `Nmodes`: the number of signal and idler frequencies.
- `Nnodes`: the number of nodes in the circuit (including the ground node).
- `Nbranches`: the number of branches in the circuit.
- `Nports`: the number of ports.
- `signalindex`: the index of the signal mode.
"""
struct LinearizedHB
    w
    modes
    S
    Snoise
    Cnoise
    Ssensitivity
    QE
    QEideal
    CM
    nodeflux
    nodefluxadjoint
    voltage
    voltageadjoint
    nodenames
    nodeindices
    componentnames
    componenttypes
    componentnamedict
    mutualinductorbranchnames
    portnumbers
    portindices
    portimpedances
    noiseportimpedanceindices
    sensitivitynames
    sensitivityindices
    Nmodes
    Nnodes
    Nbranches
    Nports
    signalindex
end

# a solution which carries no added noise covariance, which is every one
# taken before `returnCnoise` existed and every one which does not ask for it
LinearizedHB(w, modes, S, Snoise, Ssensitivity, QE, QEideal, CM, nodeflux, nodefluxadjoint, voltage, voltageadjoint, nodenames, nodeindices, componentnames, componenttypes, componentnamedict, mutualinductorbranchnames, portnumbers, portindices, portimpedances, noiseportimpedanceindices, sensitivitynames, sensitivityindices, Nmodes, Nnodes, Nbranches, Nports, signalindex) =
    LinearizedHB(w, modes, S, Snoise, Array{Complex{Float64},3}(undef, 0, 0, 0), Ssensitivity, QE, QEideal, CM, nodeflux, nodefluxadjoint, voltage, voltageadjoint, nodenames, nodeindices, componentnames, componenttypes, componentnamedict, mutualinductorbranchnames, portnumbers, portindices, portimpedances, noiseportimpedanceindices, sensitivitynames, sensitivityindices, Nmodes, Nnodes, Nbranches, Nports, signalindex)

"""
    HB(nonlinear, linearized)

A simple structure to hold the nonlinear and linearized harmonic balance
solutions.

# Fields
- `nonlinear`: nonlinear harmonic balance solution for pump and pump
    harmonics. See [`NonlinearHB`](@ref).
- `linearized`: linearized harmonic balance solution.
    See [`LinearizedHB`](@ref).
"""
struct HB
    nonlinear
    linearized
end

# ---------------------------------------------------------------------------
# Shared docstring fragments, interpolated into the docstrings of `hbsolve`,
# `hbnlsolve` and `hblinsolve` so the two descriptions of the same keyword
# cannot drift.
# ---------------------------------------------------------------------------
const _DOC_FTOL = """
- `ftol = 1e-8`: the residual tolerance `norm(F) <= ftol` at which the
    nonlinear solution is considered converged. `F` is scaled by `Z0/w0`;
    see [`calcsolverscale`](@ref). A tolerance below the rounding error of
    the scaled source, which a circuit whose interior is far from its port
    impedances can reach, is raised to it: no iteration resolves a residual
    below the rounding of the terms which make it."""

const _DOC_METHOD = """
- `method = :newtonkrylov`: the nonlinear solver, a symbol or a solver
    value carrying its own options ([`NewtonKrylov`](@ref),
    [`Newton`](@ref), [`QuasiNewton`](@ref), [`ExternalSolver`](@ref)).
    `:quasinewton` uses the complex holomorphic Jacobian `Jx` only, an
    approximation to the full Jacobian. `:newton` solves the equivalent
    real system with the full Jacobian. `:newtonkrylov` uses the
    matrix-free real Jacobian.
    `:staged` runs [`stagedhbnlsolve`](@ref): source continuation on an
    adaptively grown harmonic grid, spending the many Newton iterations a
    near critical drive needs on small cheap truncations and warm starting
    each larger grid from the last, with `:newtonkrylov` solving every
    stage. It is the strategy for operating points the direct methods fail
    outright, and the one that distinguishes a hard operating point from a
    nonexistent one: a stall on the finest grid with the minimum drive step
    brackets a fold -- the self oscillation threshold -- and errors with
    the bracketed drive fraction rather than returning garbage."""

const _DOC_STAGEDKWARGS = """
- `stagedkwargs::NamedTuple = (;)`: options for `method = :staged`; see
    [`stagedhbnlsolve`](@ref) (`grids`, `s0`, `smin`, `interiorftol`,
    `interioriterations`, `innermethod`, `verbose`)."""

const _DOC_ANDERSON = """
- `andersondepth::Integer = method == :quasinewton ? 5 : 0`: the depth of the
    Anderson acceleration of the Newton fixed point iteration, the maximum
    number of previous iterates used for the extrapolation. Values less than
    one disable the acceleration."""

const _DOC_SORTING = """
- `sorting = :number`: sort the nodes by:
    `:name`: Sort the vector of strings. This always works but leads
    to results like "101" comes before "11".
    `:number`: Convert the node strings to integer and sort by these
    (this errors if the nodes names cannot be converted to integers).
    `:none`: Don't perform any sorting except to place the ground node
    first. In other words, order the nodes in the order they are found in
    `circuit`."""

"""
    SolverInfo(stages, initialresidual, finalresidual, converged,
        sourcefold)

Diagnostics describing the nonlinear solution process of
[`hbnlsolve`](@ref).

# Fields
- `stages`: a vector of per-stage records (subtypes of
    `AbstractStageInfo`), one for each invocation of the nonlinear solver,
    in the order they ran. The direct and Krylov solvers push
    [`IterationInfo`](@ref); `method = :staged` pushes one
    [`StagedStageInfo`](@ref) per attempted continuation stage, each
    carrying its inner solver records. Every record has `label`,
    `converged` and `iterations` fields; the rest is method specific.
- `initialresidual`: the norm of the residual at the initial value.
- `finalresidual`: the norm of the residual at the returned solution.
- `converged`: whether the solver reported convergence.
- `sourcefold`: reserved for source stepping continuation, which is not
    currently implemented, and always `NaN`. It is retained so the field
    layout of this struct does not change; do not depend on it.
"""
struct SolverInfo
    stages
    initialresidual
    finalresidual
    converged
    sourcefold
end


"""
    hbsolve(ws, wp::NTuple{N,Number}, sources::Vector,
        Nmodulationharmonics::NTuple{M,Int}, Npumpharmonics::NTuple{N,Int},
        circuit, circuitdefs;dc = false, threewavemixing = false,
        fourwavemixing = true, maxintermodorder=Inf, iterations = 1000,
        ftol = 1e-8, symfreqvar = nothing, nbatches = Base.Threads.nthreads(),
        sorting = :number, returnS = true, returnSnoise = false, returnQE = true,
        returnCM = true, returnnodeflux = false, returnvoltage = false,
        returnnodefluxadjoint = false, returnvoltageadjoint = false,
        keyedarrays = true, sensitivitynames::Vector{String} = String[],
        sensitivityoperatingpoint = true, sensitivitymode = :auto,
        returnSsensitivity = false,
        factorization = nothing, backend = CPU()) where {N,M}

Calls the harmonic balance solvers, [`hbnlsolve`](@ref) and
[`hblinsolve`](@ref), which work for an arbitrary number of modes and ports,
and for both three and four wave mixing processes. See also [`hbnlsolve`](@ref)
and [`hblinsolve`](@ref).

The system is solved in a modified nodal analysis formulation in the node
flux basis: resistors with constant real values and mutually coupled
inductor branches are assigned auxiliary branch current variables (keeping
the system matrix bounded as coupling coefficients approach one), and one
gauge fixing equation per floating inductive/Josephson subnetwork and
zero-frequency mode makes circuits with no inductive path to ground
exactly solvable in the nonlinear solve. The nonlinear system is
nondimensionalized by the scale `Z0/w0` (see [`calcsolverscale`](@ref)),
making `ftol` independent of the unit system. The
linearized solve throws an informative `ArgumentError` when any signal plus
pump mode frequency total is (numerically) zero; estimate DC limits from a
sequence of decreasing nonzero frequencies. All returned quantities contain
only the node coordinates. See `src/mna.jl`.

# Arguments
- `ws`: the angular frequency or frequencies of the signal in Hz such as
    2\\*pi\\*5.0e9 or 2\\*pi\\*(4.5:0.001:5.0)\\*1e9.
- `wp::NTuple{N,Number}`: a tuple containing the angular frequencies of the
    strong tones (or pumps) such as (2\\*pi\\*5.0e9,) for a single pump at 5 GHz
    (2\\*pi\\*5.0e9,2\\*pi\\*6.0e9) for a pump at 5 GHz and a pump at 6 GHz. The
    frequencies should be non-commensurate. For commensurate pumps, the lowest
    pump frequency should be provided here, and the other pumps added to
    `sources` with a mode index equal to the ratio.
- `sources::Vector`: a vector of named tuples specifying the mode index,
    port, and current for each source. The named tuple(s) have names
    mode, port, and current. mode is a tuple specifying the mode or harmonic
    indices of the pumps, port is an integer specifying the port, and current
    is a number specifying the current. Note that the current is a complex
    number 
    For example:
    [(mode=(1,0),port=1,current=Ip1),(mode=(0,1),port=1,current=Ip2)]
    specifies two pumps where the frequency of the first pump would be
    1\\*wp1 + 0\\*wp2 and the second 0\\*wp1+1\\*wp2 where wp1 is the first
    pump frequency and wp2 is the second pump frequency. Both of the pumps are
    applied to port 1 with currents Ip1 and Ip2, respectively. 
- `Nmodulationharmonics::NTuple{M,Int}`: a tuple of integers describing how
    many signal and idler modes.
- `Npumpharmonics::NTuple{N,Int}`: a tuple of integers describing how many
    harmonics to simulate for each of the pumps. The length of the tuple must
    equal the number of non-commensurate pumps.
- `circuit`: vector of tuples each of which contain the component name, the
    first node, the second node, and the component value. The first three must
    be strings.
- `circuitdefs`: a dictionary where the keys are symbols or symbolic
    variables for component values and the values are the numerical values
    for the components.

# Keywords
- `dc = false`: include 0 frequency terms in the harmonic balance analysis.
- `threewavemixing = false`: simulate three wave mixing processes. 
- `fourwavemixing = true`: simulate four wave mixing processes.
- `maxintermodorder=Inf`: the maximum intermod order as defined by the sum of
    the absolute values of the integers multiplying each of the frequencies
    being less than or equal to `maxintermodorder`. This performs a diamond
    truncation of the discrete Fourier space.
- `iterations = 1000`: the number of iterations before the nonlinear solver
    returns an error.
$(_DOC_FTOL)
$(_DOC_METHOD)
$(_DOC_STAGEDKWARGS)
- `krylovcouplingmodes = :none`: forwarded to `hbnlsolve`; see there.
- `krylovkwargs::NamedTuple = (;)`: forwarded to `hbnlsolve`; see there.
$(_DOC_ANDERSON)
- `symfreqvar = nothing`: the symbolic frequency variable, eg `w`.
- `nbatches = Base.Threads.nthreads()`: the number of batches to split the
    signal frequencies into for multi-threading. Set to 1 for singled threaded
    evaluation.
$(_DOC_SORTING)
- `returnS = true`: return the scattering parameters from the linearized
    simulations.
- `returnSnoise = false`: return the noise scattering parameters from the
    linearized simulations.
- `returnCnoise = false`: return the added noise covariance at the output
    ports, `sum_c occupation[c] Snoise[c,i] conj(Snoise[c,j])`. With the
    scattering matrix this is the Gaussian channel the circuit implements,
    taking an input covariance to `S sigma S' + Cnoise`, and it is where the
    temperature of each channel enters.
- `returnQE = true`: return the quantum efficiency from the linearized
    simulations.
- `returnCM = true`: return the commutation relations from the linearized
    simulations.
- `returnnodeflux = false`: return the node fluxes from the linearized
    simulations.
- `returnvoltage = false`: return the node voltages from the linearized
    simulations.
- `returnnodefluxadjoint = false`: return the node fluxes from the linearized
    adjoint simulations.
- `returnvoltageadjoint = false`: return the node voltages from the linearized
    adjoint simulations.
- `keyedarrays = true`: when true return the output matrices
    and vectors as keyed arrays for more intuitive indexing. When false
    return normal matrices and vectors.
- `temperature = 0.0`: the physical temperature in Kelvin of every
    dissipative element, and so of the noise it adds. A channel at
    temperature `T` carries `coth(hbar*w/(2*k*T))` times its vacuum noise,
    which at the default of zero temperature is the vacuum noise itself.
    Raising it lowers the quantum efficiency and leaves `Snoise` alone:
    `Snoise` is a scattering matrix, describing a transformation rather than
    the state of anything, and the occupation is applied where the noise
    power is asked for. The commutation relations are likewise a statement
    about the transformation, so they stay at one at every temperature and
    remain a check on the adjoint solve, the covariance factorization and the
    port normalization.

    The explicitly defined ports are vacuum by definition. A non-vacuum input
    at a port is had by tracing that port out and assuming one there.
    A component may state its own temperature instead, and then takes that:
    a lumped one as `Resistor(R; temperature = T)`, a
    [`ScatteringParameters`](@ref) as `noise = ThermalEquilibrium(T)`. Only the
    typed circuit format carries those; a netlist of tuples states none, and
    everything in it takes this default. A block with a
    [`NoiseCovariance`](@ref) takes no temperature at all: its covariance is
    given rather than derived, and scaling it would modify an answer the
    caller already has.
- `sensitivitynames::Vector{String} = String[]`: the component names for which
    to calculate sensitivities. Supported component types are C, L, R and Lj
    with numeric values. A circuit may contain a [`ScatteringParameters`](@ref)
    while the sensitivity is taken with respect to its lumped components; a
    block itself has no scalar value to perturb, so naming one is an error.
- `sensitivitymode = :auto`: the order in which the contribution of the pump
    operating point shift is contracted. Both orders run with the pump solve
    and the signal sweep on the host, on a backend, or split between them. `:forward` costs one product against
    the sparsity structure of the linearized system per component and per
    signal frequency; `:reverse` pushes the output functionals through the
    transposed pump Jacobian once per output port mode pair, leaving a sparse
    inner product per component, so its cost does not grow with the number of
    components. `:auto` selects `:reverse` when there are more components
    than output port mode pairs, which is an operation count rather than a
    measured crossover. Both orders support any number of pump tones.
- `sensitivityoperatingpoint = true`: include the shift of the pump
    operating point in the sensitivities, making the result the total
    derivative rather than the derivative at a fixed operating point. Near
    the gain peak of a strongly pumped amplifier the operating point
    contribution is comparable to or larger than the fixed operating point
    term, so this is the relevant quantity for optimization and tolerance
    analysis. Requires the exact real Jacobian of the nonlinear solution,
    which is assembled and factorized once.
- `returnSsensitivity = false`: return `dS/dr`, the derivative of the
    scattering matrix with respect to a relative (logarithmic) perturbation
    `r` of each component value in `sensitivitynames` (`p -> r*p`, evaluated
    at `r = 1`), calculated with the adjoint method. The pump operating point
    is held fixed unless `sensitivityoperatingpoint = true`.
- `factorization = KLUfactorization()`: the type of factorization to use for
    the nonlinear and the linearized simulations.
- `backend = CPU()`: the backend both solves run on. The nonlinear solve
    assembles, factorizes and iterates there; the linearized sweep solves its
    signal frequencies there in batches which share a sparsity pattern, and
    falls back to the host for the requests it cannot serve there (see
    [`hblinsolve`](@ref)).

# Returns
- `HB`: A simple structure to hold the harmonic balance solutions. See
    [`HB`](@ref).

"""
function hbsolve(ws, wp::NTuple{N,Number}, sources::Vector,
    Nmodulationharmonics::NTuple{M,Int}, Npumpharmonics::NTuple{N,Int},
    circuit, circuitdefs;dc::Bool = false, threewavemixing::Bool = false,
    fourwavemixing::Bool = true, maxpumpintermodorder=Inf,
    maxmodulationintermodorder=Inf,
    maxpumpharmonics::NTuple{N,Number} = Npumpharmonics,
    maxmodulationharmonics::NTuple{M,Number} = Nmodulationharmonics,
    iterations = 1000, ftol = 1e-8, switchofflinesearchtol = nothing,
    alphamin = nothing, method = :newtonkrylov,
    andersondepth::Integer = method == :quasinewton ? 5 : 0,
    x0 = nothing, symfreqvar = nothing, nbatches = Base.Threads.nthreads(),
    sorting = :number, returnS::Bool = true, returnSnoise::Bool = false,
    returnQE::Bool = true, returnCM::Bool = true, returnnodeflux::Bool = false,
    returnvoltage::Bool = false, returnnodefluxadjoint::Bool = false,
    returnvoltageadjoint::Bool = false, keyedarrays::Bool = true,
    temperature = 0.0, returnCnoise::Bool = false,
    sensitivitynames::Vector{String} = String[],
    sensitivitypairs::AbstractVector =
        Tuple{String,Int,Complex{Float64}}[],
    sensitivityblockpairs::AbstractVector = Tuple{String,Int,Any}[],
    nsensitivityparameters::Integer = 0,
    sensitivitylabels::Union{Nothing,Vector{String}} = nothing,
    sensitivityoperatingpoint::Bool = true,
    sensitivitymode::Symbol = :auto,
    returnSsensitivity::Bool = false, returnZ = nothing,
    returnZadjoint = nothing, returnZsensitivity = nothing,
    returnZsensitivityadjoint = nothing,
    krylovcouplingmodes = :none, krylovkwargs::NamedTuple = (;),
    stagedkwargs::NamedTuple = (;),
    factorization = nothing, backend = CPU(),
    precision::Type{<:AbstractFloat} = Float64) where {N,M}

    # calculate the Frequencies struct
    freq = removeconjfreqs(
        truncfreqs(
            calcfreqsrdft(Npumpharmonics); dc = dc, odd = fourwavemixing,
            even = threewavemixing,
            maxintermodorder = maxpumpintermodorder,
            maxharmonics = maxpumpharmonics,
        )
    )

    indices = fourierindices(freq)

    Nmodes = length(freq.modes)

    # parse and sort the circuit
    # parse, graph, and assemble; a typed circuit takes the compiled path
    psc, cg, nm = preparecircuit(circuit, circuitdefs;
        sorting = sorting, Nmodes = Nmodes)


    # solve the nonlinear problem. `:staged` goes through the source
    # continuation driver, which re-derives its truncations from the same
    # arguments the pump truncation above was built from and manages its own
    # warm starts; every other method solves the assembled system directly.
    nonlinear = if method === :staged
        stagedhbnlsolve(wp, Npumpharmonics, sources, circuit, circuitdefs;
            iterations = iterations, ftol = ftol,
            maxharmonics = maxpumpharmonics,
            maxintermodorder = maxpumpintermodorder,
            dc = dc, odd = fourwavemixing, even = threewavemixing,
            symfreqvar = symfreqvar, sorting = sorting,
            keyedarrays = keyedarrays, sensitivitynames = sensitivitynames,
            returnoperatingpoint = sensitivityoperatingpoint &&
                returnSsensitivity && !isempty(nm.Ljb.nzind),
            krylovcouplingmodes = krylovcouplingmodes,
            krylovkwargs = krylovkwargs, factorization = factorization,
            backend = backend, precision = precision, stagedkwargs...)
    else
        hbnlsolve(wp, sources, freq, indices, psc, cg, nm;
            iterations = iterations, x0 = x0, ftol = ftol,
            switchofflinesearchtol = switchofflinesearchtol,
            alphamin = alphamin,
            method = method, andersondepth = andersondepth,
            krylovcouplingmodes = krylovcouplingmodes,
            krylovkwargs = krylovkwargs,
                symfreqvar = symfreqvar, keyedarrays = keyedarrays,
            sensitivitynames = sensitivitynames,
            returnoperatingpoint = sensitivityoperatingpoint &&
                returnSsensitivity && !isempty(nm.Ljb.nzind),
            factorization = factorization, backend = backend,
            precision = precision)
    end

    # the derivative of the harmonic balance residual with respect to each
    # component value, for total (operating point inclusive) sensitivities.
    # The frequencies of the nonlinear solution are the truncated `freq`, so
    # the numeric matrices `nm` built above have the mode count
    # `calcresidualsensitivity` expects and everything is reused. The
    # operating point derivatives dx/dr themselves are computed inside
    # `hblinsolve` only if the forward contraction order is selected there,
    # so the reverse order does not pay for the per-component solves.
    # Without Josephson junctions the linearized system matrix has no
    # dependence on the solved operating point (its only state dependence
    # is the junction modulation term), so the operating point contribution
    # is identically zero and the fixed component stamps are already the
    # total derivative: skip the residual derivatives entirely.
    sensitivityresidual = if sensitivityoperatingpoint && returnSsensitivity &&
            (!isempty(sensitivitynames) || !isempty(sensitivitypairs) ||
             !isempty(sensitivityblockpairs)) &&
            !isempty(nm.Ljb.nzind)
        cols = if isempty(sensitivitypairs) && isempty(sensitivityblockpairs)
            calcresidualsensitivity(nonlinear.operatingpoint, psc, cg, nm,
                [psc.componentnamedict[name] for name in sensitivitynames])
        elseif isempty(sensitivitypairs)
            nothing
        else
            # per-pair columns with the design parameter direction folded
            # in; hblinsolve merges them per parameter with the stamps
            calcresidualsensitivity(nonlinear.operatingpoint, psc, cg, nm,
                [psc.componentnamedict[String(t[1])]
                 for t in sensitivitypairs],
                [Complex{Float64}(t[3]) for t in sensitivitypairs])
        end
        # the scattering block pairs contribute their own columns, after
        # the lumped ones, matching the pair order of hblinsolve
        if !isempty(sensitivityblockpairs)
            blockcols = calcblockresidualsensitivity(
                nonlinear.operatingpoint, psc, sensitivityblockpairs)
            cols = isnothing(cols) ? blockcols : hcat(cols, blockcols)
        end
        cols
    else
        nothing
    end

    # generate the signal modes
    signalfreq =truncfreqs(
        calcfreqsdft(Nmodulationharmonics); dc = true, odd = threewavemixing,
        even = fourwavemixing, maxintermodorder = maxmodulationintermodorder,
        maxharmonics = maxmodulationharmonics,
    )

    # solve the linearized problem. the component values were already
    # resolved for the nonlinear solve, so the signal side numeric matrices
    # are built from them directly instead of redoing the symbolic value
    # resolution at the signal mode count.
    linearized = hblinsolve(ws, psc, cg, nm.vvn, signalfreq;
        nonlinear = nonlinear, symfreqvar = symfreqvar, nbatches = nbatches,
        returnS = returnS, returnSnoise = returnSnoise, returnQE = returnQE,
        returnCM = returnCM, returnnodeflux = returnnodeflux,
        returnnodefluxadjoint = returnnodefluxadjoint,
        returnvoltage = returnvoltage,
        returnvoltageadjoint = returnvoltageadjoint, 
        keyedarrays = keyedarrays, temperature = temperature,
        returnCnoise = returnCnoise, sensitivitynames = sensitivitynames,
        sensitivitypairs = sensitivitypairs,
        sensitivityblockpairs = sensitivityblockpairs,
        nsensitivityparameters = nsensitivityparameters,
        sensitivitylabels = sensitivitylabels,
        sensitivityresidual = sensitivityresidual,
        sensitivitymode = sensitivitymode,
        returnSsensitivity = returnSsensitivity, returnZ = returnZ,
        returnZadjoint = returnZadjoint,
        returnZsensitivity = returnZsensitivity,
        returnZsensitivityadjoint = returnZsensitivityadjoint,
        # the linearized sweep solves its frequencies on the backend when it
        # can (see hblinsolve), and falls back to the host otherwise, so it
        # takes the backend as well. The host factorization stays the default
        # because that fallback is a host solve of complex matrices; an
        # explicit one is still honored for both solves.
        factorization = isnothing(factorization) ? KLUfactorization() :
            factorization,
        backend = backend)

    return HB(nonlinear, linearized)
end

# A fully numeric circuit (such as the output of a circuit builder) needs no
# component definitions, so `circuitdefs` is optional.
function hbsolve(ws, wp::NTuple{N,Number}, sources::Vector,
    Nmodulationharmonics::NTuple{M,Int}, Npumpharmonics::NTuple{N,Int},
    circuit; kwargs...) where {N,M}
    return hbsolve(ws, wp, sources, Nmodulationharmonics, Npumpharmonics,
        circuit, Dict{Any,Any}(); kwargs...)
end

"""
    hbsolve(ws, wp, sources, Nmodulationharmonics, Npumpharmonics,
        circuit::Circuit, circuitdefs = Dict{Symbol,Number}();
        sorting = :name, keyword arguments...)

Harmonic balance solution of a typed [`Circuit`](@ref). The circuit is
elaborated and lowered with [`compile`](@ref) and solved with the
existing solver; all keyword arguments of the legacy method are supported.
`circuitdefs` is only needed when component values are symbolic. The
default `sorting` is `:name` because hierarchical net names are not
integers.
"""
function hbsolve(ws, wp::NTuple{N,Number}, sources::Vector,
        Nmodulationharmonics::NTuple{M,Int}, Npumpharmonics::NTuple{N,Int},
        circuit::Circuit, circuitdefs::AbstractDict = Dict{Symbol,Number}();
        sorting::Symbol = :name, kwargs...) where {N,M}
    return hbsolve(ws, wp, sources, Nmodulationharmonics, Npumpharmonics,
        elaborate(circuit), circuitdefs; sorting = sorting, kwargs...)
end
