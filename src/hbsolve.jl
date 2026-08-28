
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
    of port and frequency.
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
    [`ScatteringBlock`](@ref). Being a scattering matrix it describes a
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
- `portimpedanceindices`:
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
    portimpedanceindices
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
LinearizedHB(w, modes, S, Snoise, Ssensitivity, QE, QEideal, CM, nodeflux, nodefluxadjoint, voltage, voltageadjoint, nodenames, nodeindices, componentnames, componenttypes, componentnamedict, mutualinductorbranchnames, portnumbers, portindices, portimpedanceindices, noiseportimpedanceindices, sensitivitynames, sensitivityindices, Nmodes, Nnodes, Nbranches, Nports, signalindex) =
    LinearizedHB(w, modes, S, Snoise, Array{Complex{Float64},3}(undef, 0, 0, 0), Ssensitivity, QE, QEideal, CM, nodeflux, nodefluxadjoint, voltage, voltageadjoint, nodenames, nodeindices, componentnames, componenttypes, componentnamedict, mutualinductorbranchnames, portnumbers, portindices, portimpedanceindices, noiseportimpedanceindices, sensitivitynames, sensitivityindices, Nmodes, Nnodes, Nbranches, Nports, signalindex)

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
- `ftol = 1e-8`: the residual tolerance `norm(F) <= ftol` at which the
    nonlinear solution is considered converged. `F` is  is scaled by `Z0/w0`.
    Aee [`calcsolverscale`](@ref)).
- `method = :newtonkrylov`: the nonlinear solver. `:quasinewton`
    uses the complex holomorphic Jacobian `Jx` only, an approximation to the
    full Jacobian. `:newton` solves the equivalent real system with the full
    Jacobian. `:newtonkrylov` uses the matrix-free real Jacobian.
    `:staged` solves the nonlinear stage by [`stagedhbnlsolve`](@ref):
    source continuation on an adaptively grown harmonic grid, for strongly
    driven operating points the direct methods cannot reach from a cold
    start; a fold (self oscillation threshold) below the requested drive is
    reported as an error with the bracketed drive fraction.
- `stagedkwargs::NamedTuple = (;)`: options for `method = :staged`; see
    [`stagedhbnlsolve`](@ref).
- `krylovcouplingmodes = :none`: forwarded to `hbnlsolve`; see there.
- `krylovkwargs::NamedTuple = (;)`: forwarded to `hbnlsolve`; see there.
- `andersondepth::Integer = method == :quasinewton ? 5 : 0`: the depth of the
    Anderson acceleration of the Newton fixed point iteration, the maximum
    number of previous iterates used for the extrapolation. Values less than
    one disable the acceleration.
- `symfreqvar = nothing`: the symbolic frequency variable, eg `w`.
- `nbatches = Base.Threads.nthreads()`: the number of batches to split the
    signal frequencies into for multi-threading. Set to 1 for singled threaded
    evaluation.
- `sorting = :number`: sort the nodes by:
    `:name`: Sort the vector of strings. This always works but leads
    to results like "101" comes before "11".
    `:number`: Convert the node strings to integer and sort by these
    (this errors if the nodes names cannot be converted to integers).
    `:none`: Don't perform any sorting except to place the ground node
    first. In other words, order the nodes in the order they are found in
    `circuit`.
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
    [`ScatteringBlock`](@ref) as `noise = ThermalEquilibrium(T)`. Only the
    typed circuit format carries those; a netlist of tuples states none, and
    everything in it takes this default. A block with a
    [`NoiseCovariance`](@ref) takes no temperature at all: its covariance is
    given rather than derived, and scaling it would modify an answer the
    caller already has.
- `sensitivitynames::Vector{String} = String[]`: the component names for which
    to calculate sensitivities. Supported component types are C, L, R and Lj
    with numeric values. A circuit may contain a [`ScatteringBlock`](@ref)
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
    psc = parsesortcircuit(circuit, sorting = sorting)

    # calculate the circuit graph
    cg = calccircuitgraph(psc)

    # calculate the numeric matrices
    nm=numericmatrices(psc, cg, circuitdefs, Nmodes = Nmodes)


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
            !isempty(sensitivitynames) && !isempty(nm.Ljb.nzind)
        calcresidualsensitivity(nonlinear.operatingpoint, psc, cg, nm,
            [psc.componentnamedict[name] for name in sensitivitynames])
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

"""
    hblinsolve(w, circuit, circuitdefs; Nmodulationharmonics = (0,),
        nonlinear = nothing, symfreqvar = nothing, threewavemixing = false,
        fourwavemixing = true, maxharmonics = Nmodulationharmonics,
        maxintermodorder = Inf,
        nbatches::Integer = Base.Threads.nthreads(), sorting = :number,
        returnS = true, returnSnoise = false, returnQE = true,
        returnCM = true, returnnodeflux = false,
        returnnodefluxadjoint = false, returnvoltage = false,
        returnvoltageadjoint = false, keyedarrays = true,
        sensitivitynames::Vector{String} = String[],
        sensitivitynodeflux = nothing, sensitivityresidual = nothing,
        sensitivitymode = :auto, returnSsensitivity = false,
        factorization = KLUfactorization(), backend = CPU())

Harmonic balance solver supporting an arbitrary number of small signals (weak
tones) linearized around `nonlinear`, the solution of the nonlinear system
consisting of an arbitrary number of large signals (strong tones) from
`hbnlsolve`.

# Arguments
- `w`: the small signal frequency or frequencies.
- `circuit`: vector of tuples each of which contain the component name, the
    first node, the second node, and the component value. The first three must
    be strings.
- `circuitdefs`: a dictionary where the keys are symbols or symbolic
    variables for component values and the values are the numerical values
    for the components.

# Keywords
- `Nmodulationharmonics::NTuple{M,Int}`: a tuple of integers describing how
    many signal and idler modes.
- `nonlinear=nothing`: solution to the nonlinear system from `hbnlsolve`.
- `symfreqvar = nothing`: the symbolic frequency variable, eg `w`.
- `threewavemixing = false`: simulate three wave mixing processes. 
- `fourwavemixing = true`: simulate four wave mixing processes.
- `maxintermodorder=Inf`: the maximum intermod order as defined by the sum of
    the absolute values of the integers multiplying each of the frequencies
    being less than or equal to `maxintermodorder`. This performs a diamond
    truncation of the discrete Fourier space.
- `nbatches = Base.Threads.nthreads()`: the number of batches to split the
    signal frequencies into for multi-threading. Set to 1 for singled threaded
    evaluation.
- `sorting = :number`: sort the nodes by:
    `:name`: Sort the vector of strings. This always works but leads
    to results like "101" comes before "11".
    `:number`: Convert the node strings to integer and sort by these
    (this errors if the nodes names cannot be converted to integers).
    `:none`: Don't perform any sorting except to place the ground node
    first. In other words, order the nodes in the order they are found in
    `circuit`.
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
    [`ScatteringBlock`](@ref) as `noise = ThermalEquilibrium(T)`. Only the
    typed circuit format carries those; a netlist of tuples states none, and
    everything in it takes this default. A block with a
    [`NoiseCovariance`](@ref) takes no temperature at all: its covariance is
    given rather than derived, and scaling it would modify an answer the
    caller already has.
- `sensitivitynames::Vector{String} = String[]`: the component names for which
    to calculate sensitivities. Supported component types are C, L, R and Lj
    with numeric values. A circuit may contain a [`ScatteringBlock`](@ref)
    while the sensitivity is taken with respect to its lumped components; a
    block itself has no scalar value to perturb, so naming one is an error.
- `sensitivitynodeflux = nothing`: the derivatives of the pump operating
    point with respect to each component value, the columns of
    [`calcnodefluxsensitivity`](@ref), to include the operating point shift
    in the sensitivities. `hbsolve` supplies the residual derivatives
    instead, from which these are computed only if the forward contraction
    order is selected.
- `sensitivityresidual = nothing`: the derivatives of the harmonic balance
    residual with respect to each component value, the columns of
    [`calcresidualsensitivity`](@ref). Required by the reverse contraction
    order; sufficient by itself for either order.
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
- `backend = CPU()`: the backend the frequency sweep is solved on. On a device
    backend the system matrices of a batch of signal frequencies, which share
    one sparsity pattern, are assembled by one kernel and factorized and solved
    as a uniform batch. Falls back to the host when the component values depend
    on the symbolic frequency variable, or when an output requires the adjoint
    (transposed) solve: the noise scattering parameters, the scattering
    parameter sensitivities, or the adjoint node outputs.

# Returns
- `LinearizedHB`: A simple structure to hold the harmonic balance solutions.
    See [`LinearizedHB`](@ref).

# Examples
```jldoctest
circuit = Tuple{String,String,String,Union{Complex{Float64},Symbol,Int64}}[]
push!(circuit,("P1","1","0",1))
push!(circuit,("R1","1","0",:Rleft))
push!(circuit,("L1","1","0",:Lm)) 
push!(circuit,("K1","L1","L2",:K1))
push!(circuit,("C1","1","2",:Cc)) 
push!(circuit,("L2","2","3",:Lm)) 
push!(circuit,("Lj3","3","0",:Lj)) 
push!(circuit,("Lj4","2","0",:Lj)) 
push!(circuit,("C2","2","0",:Cj))
circuitdefs = Dict{Symbol,Complex{Float64}}(
    :Lj =>2000e-12,
    :Lm =>10e-12,
    :Cc => 200.0e-15,
    :Cj => 900e-15,
    :Rleft => 50.0,
    :Rright => 50.0,
    :K1 => 0.9,
)

Idc = 1e-6*0
Ip=5.0e-6
wp=2*pi*5e9
ws=2*pi*5.2e9
symfreqvar = nothing

# modulation settings
Npumpharmonics = (16,)
Nmodulationharmonics = (2,)
threewavemixing=false
fourwavemixing=true

nonlinear=hbnlsolve(
    (wp,),
    Npumpharmonics,
    [
        (mode=(0,),port=1,current=Idc),
        (mode=(1,),port=1,current=Ip),
    ],
    circuit,circuitdefs;dc=true,odd=fourwavemixing,even=threewavemixing)

linearized = JosephsonCircuits.hblinsolve(ws,
    circuit, circuitdefs; Nmodulationharmonics = Nmodulationharmonics,
    nonlinear = nonlinear, symfreqvar=nothing, threewavemixing=false,
    fourwavemixing=true, returnnodeflux=true, keyedarrays = false)
isapprox(linearized.nodeflux,
    ComplexF64[9.901008591291e-12 - 6.40587007644028e-14im 2.164688307719963e-14 - 2.90852607344097e-16im 6.671563044645655e-14 - 8.585524364135119e-16im; 2.1633104519765224e-14 - 8.251861334047893e-16im 1.0099063486905209e-11 - 1.948847859339803e-13im -8.532003011745068e-15 + 3.234788465760295e-16im; 6.671648606599472e-14 + 7.892709980649199e-16im -8.53757633177974e-15 - 9.748395563374129e-17im 9.856580758892428e-12 + 5.859984004390703e-14im; 1.5888896262186103e-11 - 1.0303480614499543e-13im -2.557126237504446e-12 + 1.759201163407723e-14im -8.475819811683215e-12 + 5.3531443609574795e-14im; -2.5781681021577177e-13 + 4.757590640631487e-15im 2.36818731889176e-12 - 4.569646499606389e-14im 1.116372367616482e-13 - 2.039935997276492e-15im; -1.0210743447568219e-11 - 5.905490368441375e-14im 1.3377918536056493e-12 + 7.190105205618706e-15im 2.5392856657302323e-11 + 1.5143842454586225e-13im; 2.4781693042536835e-11 - 1.6057018472176702e-13im -2.5342360504077476e-12 + 1.7306764301173096e-14im -8.40554044664581e-12 + 5.269404591748149e-14im; -2.348528974341763e-13 + 3.949450668269274e-15im 1.1449271118157543e-11 - 2.2093702114766968e-13im 1.0261871618968225e-13 - 1.7240213938923877e-15im; -1.0140560031409567e-11 - 5.828587508192886e-14im 1.3288225860409326e-12 + 7.0954601524623594e-15im 3.423954321087654e-11 + 2.0403371894291513e-13im],
    atol = 1e-6)

# output
true
```
"""
function hblinsolve(w, circuit,circuitdefs; Nmodulationharmonics = (0,),
    nonlinear = nothing, symfreqvar = nothing, threewavemixing::Bool = false,
    fourwavemixing::Bool = true,
    maxharmonics = Nmodulationharmonics, maxintermodorder = Inf,
    nbatches::Integer = Base.Threads.nthreads(), sorting = :number,
    returnS::Bool = true, returnSnoise::Bool = false, returnQE::Bool = true,
    returnCM::Bool = true, returnnodeflux::Bool = false,
    returnnodefluxadjoint::Bool = false, returnvoltage::Bool = false,
    returnvoltageadjoint::Bool = false, keyedarrays::Bool = true,
    temperature = 0.0, returnCnoise::Bool = false,
    sensitivitynames::Vector{String} = String[],
    sensitivitynodeflux = nothing, sensitivityresidual = nothing,
    sensitivitymode::Symbol = :auto,
    returnSsensitivity::Bool = false, returnZ = nothing,
    returnZadjoint = nothing, returnZsensitivity = nothing,
    returnZsensitivityadjoint = nothing,
    factorization = KLUfactorization(), backend = CPU())

    # parse and sort the circuit
    psc = parsesortcircuit(circuit, sorting = sorting)

    # calculate the circuit graph
    cg = calccircuitgraph(psc)

    # generate the signal modes
    signalfreq =truncfreqs(
        calcfreqsdft(Nmodulationharmonics); dc = true, odd = threewavemixing,
        even = fourwavemixing, maxintermodorder = maxintermodorder,
        maxharmonics = maxharmonics,
    )

return hblinsolve(w, psc, cg, circuitdefs, signalfreq; nonlinear = nonlinear,
        symfreqvar = symfreqvar, nbatches = nbatches,
        returnS = returnS, returnSnoise = returnSnoise, returnQE = returnQE,
        returnCM = returnCM, returnnodeflux = returnnodeflux,
        returnnodefluxadjoint = returnnodefluxadjoint,
        returnvoltage = returnvoltage,
        returnvoltageadjoint = returnvoltageadjoint,
        keyedarrays = keyedarrays, temperature = temperature,
        returnCnoise = returnCnoise, sensitivitynames = sensitivitynames,
        sensitivitynodeflux = sensitivitynodeflux,
        sensitivityresidual = sensitivityresidual,
        sensitivitymode = sensitivitymode,
        returnSsensitivity = returnSsensitivity,returnZ = returnZ,
        returnZadjoint = returnZadjoint,
        returnZsensitivity = returnZsensitivity,
        returnZsensitivityadjoint = returnZsensitivityadjoint,
        factorization = factorization, backend = backend)
end

"""
    hblinsolve(w, psc::ParsedSortedCircuit,
        cg::CircuitGraph, circuitdefs, signalfreq::Frequencies{N};
        nonlinear=nothing, symfreqvar=nothing,
        nbatches::Integer = Base.Threads.nthreads(), sorting = :number,
        returnS = true, returnSnoise = false, returnQE = true, returnCM = true,
        returnnodeflux = false, returnnodefluxadjoint = false,
        returnvoltage = false,
        )

Harmonic balance solver supporting an arbitrary number of small signals (weak
tones) linearized around `pump`, the solution of the nonlinear system consisting
of an arbitrary number of large signals (strong tones).

# Examples
```jldoctest
circuit = Tuple{String,String,String,Union{Complex{Float64},Symbol,Int64}}[]
push!(circuit,("P1","1","0",1))
push!(circuit,("R1","1","0",:Rleft))
push!(circuit,("L1","1","0",:Lm)) 
push!(circuit,("K1","L1","L2",:K1))
push!(circuit,("C1","1","2",:Cc)) 
push!(circuit,("L2","2","3",:Lm)) 
push!(circuit,("Lj3","3","0",:Lj)) 
push!(circuit,("Lj4","2","0",:Lj)) 
push!(circuit,("C2","2","0",:Cj))
circuitdefs = Dict{Symbol,Complex{Float64}}(
    :Lj =>2000e-12,
    :Lm =>10e-12,
    :Cc => 200.0e-15,
    :Cj => 900e-15,
    :Rleft => 50.0,
    :Rright => 50.0,
    :K1 => 0.9,
)

Idc = 1e-6*0
Ip = 5.0e-6
wp = 2*pi*5e9
ws = 2*pi*5.2e9
Npumpharmonics = (2,)
Nmodulationharmonics = (2,)
threewavemixing = false
fourwavemixing = true

frequencies = JosephsonCircuits.removeconjfreqs(
    JosephsonCircuits.truncfreqs(
        JosephsonCircuits.calcfreqsrdft(Npumpharmonics),
        dc = true, odd = true, even = false, maxintermodorder = Inf,
    )
)
fi = JosephsonCircuits.fourierindices(frequencies)
Nmodes = length(frequencies.modes)
psc = JosephsonCircuits.parsesortcircuit(circuit)
cg = JosephsonCircuits.calccircuitgraph(psc)
nm = JosephsonCircuits.numericmatrices(psc, cg, circuitdefs, Nmodes = Nmodes)
nonlinear = hbnlsolve(
    (wp,),
    [
        (mode=(0,),port=1,current=Idc),
        (mode=(1,),port=1,current=Ip),
    ],
    frequencies, fi, psc, cg, nm)
signalfreq =JosephsonCircuits.truncfreqs(
    JosephsonCircuits.calcfreqsdft(Nmodulationharmonics),
    dc = true, odd = threewavemixing, even = fourwavemixing,
    maxintermodorder = Inf,
)
linearized = JosephsonCircuits.hblinsolve(ws, psc, cg, circuitdefs,
    signalfreq;nonlinear = nonlinear, returnnodeflux=true, keyedarrays = false)
isapprox(linearized.nodeflux,
    ComplexF64[9.901008591291e-12 - 6.40587007644028e-14im 2.164688307719963e-14 - 2.90852607344097e-16im 6.671563044645655e-14 - 8.585524364135119e-16im; 2.1633104519765224e-14 - 8.251861334047893e-16im 1.0099063486905209e-11 - 1.948847859339803e-13im -8.532003011745068e-15 + 3.234788465760295e-16im; 6.671648606599472e-14 + 7.892709980649199e-16im -8.53757633177974e-15 - 9.748395563374129e-17im 9.856580758892428e-12 + 5.859984004390703e-14im; 1.5888896262186103e-11 - 1.0303480614499543e-13im -2.557126237504446e-12 + 1.759201163407723e-14im -8.475819811683215e-12 + 5.3531443609574795e-14im; -2.5781681021577177e-13 + 4.757590640631487e-15im 2.36818731889176e-12 - 4.569646499606389e-14im 1.116372367616482e-13 - 2.039935997276492e-15im; -1.0210743447568219e-11 - 5.905490368441375e-14im 1.3377918536056493e-12 + 7.190105205618706e-15im 2.5392856657302323e-11 + 1.5143842454586225e-13im; 2.4781693042536835e-11 - 1.6057018472176702e-13im -2.5342360504077476e-12 + 1.7306764301173096e-14im -8.40554044664581e-12 + 5.269404591748149e-14im; -2.348528974341763e-13 + 3.949450668269274e-15im 1.1449271118157543e-11 - 2.2093702114766968e-13im 1.0261871618968225e-13 - 1.7240213938923877e-15im; -1.0140560031409567e-11 - 5.828587508192886e-14im 1.3288225860409326e-12 + 7.0954601524623594e-15im 3.423954321087654e-11 + 2.0403371894291513e-13im],
    atol = 1e-6)

# output
true
```
"""
function hblinsolve(w, psc::ParsedSortedCircuit,
    cg::CircuitGraph, circuitdefs, signalfreq::Frequencies; nonlinear = nothing,
    symfreqvar = nothing, nbatches::Integer = Base.Threads.nthreads(),
    returnS::Bool = true, returnSnoise::Bool = false, returnQE::Bool = true,
    returnCM::Bool = true, returnnodeflux::Bool = false,
    returnnodefluxadjoint::Bool = false, returnvoltage::Bool = false,
    returnvoltageadjoint::Bool = false, keyedarrays::Bool = true,
    temperature = 0.0, returnCnoise::Bool = false,
    sensitivitynames::Vector{String} = String[],
    sensitivitynodeflux = nothing, sensitivityresidual = nothing,
    sensitivitymode::Symbol = :auto,
    returnSsensitivity::Bool = false, returnZ = nothing,
    returnZadjoint = nothing, returnZsensitivity = nothing,
    returnZsensitivityadjoint = nothing,
    factorization = KLUfactorization(), backend = CPU(), debuglsys = false)


    # deprecation warnings for  `returnZ`, `returnZadjoint`,
    # `returnZsensitivity`, and `returnZsensitivityadjoint`
    if !isnothing(returnZ) || !isnothing(returnZadjoint) || 
        !isnothing(returnZsensitivity) || !isnothing(returnZsensitivityadjoint)
        Base.depwarn(lazy"The `returnZ`, `returnZadjoint`, `returnZsensitivity`, and `returnZsensitivityadjoint` kwargs have been removed. Please compute them from scattering parameters matrices.", :hblinsolve; force=true)
    end

    Nsignalmodes = length(signalfreq.modes)
    # calculate the numeric matrices
    signalnm = numericmatrices(psc, cg, circuitdefs, Nmodes = Nsignalmodes)

    if isnothing(nonlinear)

        allpumpfreq = calcfreqsrdft((0,))
        Amatrixmodes, Amatrixindices = hbmatind(allpumpfreq, signalfreq)
        Nwtuple = NTuple{length(allpumpfreq.Nw)+1,Int}((allpumpfreq.Nw..., length(signalnm.Ljb.nzval)))
        phimatrix = ones(Complex{Float64}, Nwtuple)
        # wpumpmodes = calcmodefreqs((0.0,),signalfreq.modes)

    else

        pumpfreq = nonlinear.frequencies

        # pumpfreq = JosephsonCircuits.calcfreqsrdft((2*Npumpmodes,))
        allpumpfreq = calcfreqsrdft(pumpfreq.Nharmonics)
        pumpindices = fourierindices(pumpfreq)
        Npumpmodes = length(pumpfreq.modes)

        Amatrixmodes, Amatrixindices = hbmatind(allpumpfreq, signalfreq)

        # calculate the dimensions of the array which holds the frequency
        # domain information for the fourier transform
        Nwtuple = NTuple{length(pumpfreq.Nw)+1,Int}((pumpfreq.Nw..., length(nonlinear.Ljb.nzval)))

        # create an array to hold the frequency domain data for the
        # fourier transform
        phimatrix = zeros(Complex{Float64}, Nwtuple)

        # create an array to hold the time domain data for the RFFT. also generate
        # the plans.
        phimatrixtd, irfftplan, rfftplan = plan_applynl(phimatrix, CPU())

        # convert the branch flux vector to a matrix with the terms arranged
        # in the correct way for the inverse rfft including the appropriate
        # complex conjugates.
        branchflux = nonlinear.Rbnm*nonlinear.nodeflux[:]
        phivectortomatrix!(
            branchflux[nonlinear.Ljbm.nzind], phimatrix,
            pumpindices.vectomatmap,
            pumpindices.conjsourceindices,
            pumpindices.conjtargetindices,
            length(nonlinear.Ljb.nzval)
        )

        # apply the sinusoidal nonlinearity when evaluaing the function
        applynl!(
            phimatrix,
            phimatrixtd,
            cos,
            irfftplan,
            rfftplan,
        )

        # wpumpmodes = calcmodefreqs(nonlinear.w,signalfreq.modes)
    end

    # to avoid wpumpmodes being boxed
    wpumpmodes::Vector{Float64} = if isnothing(nonlinear)
        calcmodefreqs((0.0,),signalfreq.modes)
    else
        calcmodefreqs(nonlinear.w,signalfreq.modes)
    end

    if !all(isfinite, w)
        throw(ArgumentError("All signal frequencies must be finite."))
    end
    # make sure the signal frequency vector or iterator isn't empty.
    if isempty(w)
        throw(ArgumentError("At least one signal frequency is required."))
    end
    # make sure there is at least one batch.
    if nbatches < 1
        throw(ArgumentError(lazy"`nbatches` = $(nbatches) must be at least 1."))
    end
    # the nonlinear solution's fields are untyped, so convert the pump
    # frequencies to a concrete tuple here: everything below is per
    # (frequency, mode) and would box every value otherwise.
    let wpumptuple = isnothing(nonlinear) ? (0.0,) :
            map(x -> Float64(real(x)), nonlinear.w)
        # the contributing term magnitudes of each mode do not depend on
        # the signal frequency, so their sum and the term count are
        # precomputed per mode; only abs(wi) joins per frequency. this
        # keeps the O(Nfrequencies*Nmodes) validation allocation free.
        all(isfinite, wpumptuple) || throw(ArgumentError("Every pump frequency must be finite."))
        modescales = Float64[sum(j -> abs(mode[j]*wpumptuple[j]),
            eachindex(wpumptuple)) for mode in signalfreq.modes]
        nterms = 1 + length(wpumptuple)
        for wi in w
            for (mi, wm) in enumerate(wpumpmodes)
                mode = signalfreq.modes[mi]
                if isnumericallyzero(wi + wm,
                        abs(float(real(wi))) + modescales[mi], nterms)
                    throw(ArgumentError("hblinsolve cannot evaluate a mode at (numerically) zero total frequency (signal frequency plus pump mode frequency, here signal $(wi/(2*pi)) Hz with mode $(mode)) because the node flux basis represents voltages as v = im*w*phi. Zero-frequency small-signal analysis is not supported; to estimate a DC limit, evaluate a sequence of decreasing nonzero frequencies and verify that the requested network parameters converge. For frequency independent resistive networks the result at any nonzero frequency equals the DC limit."))
                end
            end
        end
    end

    # this is the first signal frequency. we will use it for various setup tasks
    wmodes = w[1] .+ wpumpmodes

    # extract the elements we need
    Nnodes = psc.Nnodes
    nodenames = psc.nodenames
    nodeindices = psc.nodeindices
    componentnames = psc.componentnames
    componentnamedict = psc.componentnamedict
    componenttypes = psc.componenttypes
    mutualinductorbranchnames = psc.mutualinductorbranchnames
    Nbranches = cg.Nbranches
    edge2indexdict = cg.edge2indexdict
    Ljb = signalnm.Ljb
    Rbnm = signalnm.Rbnm
    Cnm = signalnm.Cnm
    Gnm = signalnm.Gnm
    invLnm = signalnm.invLnm

    #error if any component value contains symbolic variables which were not
    # assigned numerical values in circuitdefs (values depending only on the
    # symbolic frequency variable are ok).
    checkcomponentvaluesdefined(psc.componentnames, signalnm.vvn,
        symfreqvar)

    # set up the modified nodal analysis (MNA), which assigns auxiliary branch
    # current variables to the port resistors with constant real values and to the
    # mutually coupled inductor branches (whose inverse inductance entries
    # would otherwise diverge as the coupling coefficient approaches one).
    # eliminating the auxiliary variables recovers the nodal equations. no
    # gauge fixing equations are needed because modes at (numerically) zero
    # total frequency are not permitted.
    checkstaticstiffnessvalues(psc.componenttypes, signalnm.vvn)
    mnaindices = mnaportresistorindices(psc.componenttypes, psc.nodeindices,
        psc.mutualinductorbranchnames, signalnm.vvn)
    coupledbranches = mnacoupledbranches(signalnm.Mb)
    Nauxmnar = length(mnaindices)*Nsignalmodes
    Nauxscattering = countscatteringports(psc)*Nsignalmodes
    Nauxmna = Nauxmnar + length(coupledbranches)*Nsignalmodes +
        Nauxscattering
    Nnodalmna = (psc.Nnodes-1)*Nsignalmodes
    Amna0, AmnaG = calcAmnasplit(mnaindices, psc.nodeindices, signalnm.vvn,
        Nsignalmodes, psc.Nnodes)
    Amna0 = mnapad(Amna0, length(coupledbranches)*Nsignalmodes +
        Nauxscattering)
    AmnaG = mnapad(AmnaG, length(coupledbranches)*Nsignalmodes +
        Nauxscattering)
    if !isempty(coupledbranches)
        # the coupled inductor branches: numericmatrices already excludes
        # them from the inverse inductance matrix; add the branch flux
        # constitutive equations and Kirchhoff current law couplings, with
        # unscaled branch currents as the auxiliary variables (Lscale = 1),
        # matching the unscaled matrices of the linearized solver.
        AmnaL = calcAmnaind(coupledbranches, signalnm.Lb, signalnm.Mb,
            cg.Rbn, Nsignalmodes, Nnodalmna + Nauxmnar,
            Nnodalmna + Nauxmna, 1)
        Amna0 = spaddkeepzeros(Amna0, AmnaL)
    end
    if !isempty(mnaindices)
        Gnmp = calcGn(psc.componenttypes[mnaindices],
            psc.nodeindices[:, mnaindices], signalnm.vvn[mnaindices],
            Nsignalmodes, psc.Nnodes)
        Gnm = mnasubtractpromoted(Gnm, Gnmp)
    end
    if Nauxmna > 0
        Cnm = mnapad(Cnm, Nauxmna)
        Gnm = mnapad(Gnm, Nauxmna)
        invLnm = mnapad(invLnm, Nauxmna)
    end
    # the incidence matrix used for the sparsity structure and the pump
    # modulation contribution gains empty columns for the auxiliary
    # variables; the nodal Rbnm is kept for the source assembly.
    Rbnmmna = hcat(Rbnm, spzeros(eltype(Rbnm), size(Rbnm,1), Nauxmna))
    portindices = signalnm.portindices
    portnumbers = signalnm.portnumbers
    portimpedanceindices = signalnm.portimpedanceindices
    vvn = signalnm.vvn
    modes = signalfreq.modes

    # find the indices associated with the components for which we will
    # calculate sensitivities
    sensitivityindices = zeros(Int,length(sensitivitynames))
    for i in eachindex(sensitivitynames)
        sensitivityindices[i] = componentnamedict[sensitivitynames[i]]
        if componenttypes[sensitivityindices[i]] == :S
            throw(ArgumentError(lazy"Sensitivities with respect to scattering block components are not supported; got $(sensitivitynames[i])."))
        end
    end

    # calculate the source currents
    Nports = length(portindices)

    # calculate the source terms in the branch basis
    bbm = zeros(Complex{Float64},Nbranches*Nsignalmodes,Nsignalmodes*Nports)

    # add a current source for each port and mode
    for (i,val) in enumerate(portindices)
        key = (nodeindices[1,val],nodeindices[2,val])
        for j = 1:Nsignalmodes
            bbm[(edge2indexdict[key]-1)*Nsignalmodes+j,(i-1)*Nsignalmodes+j] = 1
            # bbm2[(i-1)*Nmodes+j,(i-1)*Nmodes+j] = Lmean*1/phi0
        end
    end

    # calculate the source terms in the node basis
    bnm = transpose(Rbnm)*bbm
    if Nauxmna > 0
        bnm = vcat(bnm, zeros(eltype(bnm), Nauxmna, size(bnm, 2)))
    end
    # return bnm
    # if there is a symbolic frequency variable, then we need to redo the noise
    # port calculation because calcnoiseportimpedanceindices() can't tell if a
    # symbolic expression is complex. 
    noiseportimpedanceindices = if isnothing(symfreqvar)
        signalnm.noiseportimpedanceindices
    else
        calcnoiseportimpedanceindices(
            componenttypes, nodeindices,
            mutualinductorbranchnames,
            Symbolics.substitute.(vvn, symfreqvar => wmodes[1]))
    end

    # find the indices at which there are symbolic variables so we can
    # perform a substitution on only those. 
    Cnmfreqsubstindices  = symbolicindices(Cnm)
    Gnmfreqsubstindices  = symbolicindices(Gnm)
    invLnmfreqsubstindices  = symbolicindices(invLnm)


    Cnmcopy = freqsubst(Cnm,wmodes,symfreqvar)
    Gnmcopy = freqsubst(Gnm,wmodes,symfreqvar)
    invLnmcopy = freqsubst(invLnm,wmodes,symfreqvar)

    # Build the linearized system object: the sparsity structure of the
    # system matrix Asparse = (AoLjnm + invLnmcopy + Gnmcopy + Cnmcopy) with
    # stored numerical zeros, a plan for scattering the Josephson (pump
    # modulation) contribution AoLjnm = Rbnm'*AoLjbm*Rbnm into it directly
    # from the Fourier coefficients of cos(phi(t)) of the pump, index maps
    # for the frequency dependent linear terms, and the precomputed pump
    # modulation contribution and its complex conjugate. This shares the
    # machinery used for the Jacobians of the nonlinear system, see
    # HBLinearizedSystem and plancomplexjacobian. The per-frequency system
    # matrices are assembled from this object with assemblesystemmatrix!.
    # the hybrid (wave to modified nodal analysis) contribution of the
    # scattering block components: constant Kirchhoff current law couplings
    # of the auxiliary port currents (folded into the constant augmentation
    # matrix) and the frequency dependent constitutive equations, assembled
    # per frequency alongside the conductance term
    ssys = scatteringstampsystem(psc, Nsignalmodes;
        auxoffset = Nnodalmna + Nauxmna - Nauxscattering,
        Ntotal = Nnodalmna + Nauxmna, scale = 1.0)
    if !isnothing(ssys)
        Amna0 = spaddkeepzeros(Amna0, ssys.kcl)
    end
    lsys = HBLinearizedSystem(Amatrixindices, signalnm.Ljb, Rbnmmna,
        Nsignalmodes, cg.Nbranches, phimatrix, invLnmcopy, Gnmcopy, Cnmcopy,
        invLnm, Gnm, Cnm, invLnmfreqsubstindices, Gnmfreqsubstindices,
        Cnmfreqsubstindices, Amna0, AmnaG, symfreqvar, wpumpmodes, Nnodes;
        scattering = ssys)
    Asparse = lsys.Asparse

    # the vacuum noise channels of the dissipative scattering blocks, which
    # follow the noise ports of the dissipative lumped components in the rows
    # of the noise scattering matrix
    # a block which declares itself lossless is taken at its word for the
    # noise, so it is worth looking at the frequencies actually being solved
    checklosslessblocks(ssys, w, wpumpmodes)
    noiseplan = if returnSnoise || returnQE || returnCM
        planscatteringnoise(ssys)
    else
        nothing
    end
    Nnoisechannels = length(noiseportimpedanceindices) +
        (isnothing(noiseplan) ? 0 : noiseplan.Nchannels)
    # the temperature of each noise channel: the analysis default, overridden
    # by name, and overridden again by a block which states its own
    channeltemperatures = noisechanneltemperatures(psc,
        noiseportimpedanceindices, noiseplan, ssys, temperature)

    portimpedances = [vvn[i] for i in portimpedanceindices]
    noiseportimpedances = [vvn[i] for i in noiseportimpedanceindices]

    # assemble Asparse once at the first frequency so we have something
    # reasonable to factorize.
    assemblesystemmatrix!(Asparse, lsys, wmodes)

    # precompute the derivative of the linearized system matrix with respect
    # to a relative perturbation of each component in sensitivitynames.
    sensitivitystamps = if returnSsensitivity
        calcsensitivitystamps(sensitivityindices, psc, cg, signalnm, lsys,
            phimatrix, mnaindices, coupledbranches, Nnodalmna, Nsignalmodes,
            Nnodes)
    else
        SensitivityStamp[]
    end

    # the contribution of the shift of the pump operating point, when the
    # operating point derivatives are supplied, either as the derivatives of
    # the operating point itself (sensitivitynodeflux) or as the residual
    # derivatives (sensitivityresidual), from which the former are computed
    # here only if the forward contraction order needs them.
    # The operating point contribution can be contracted in either order.
    # The forward order costs one product against the full sparsity structure
    # of the system matrix per component and per frequency; the reverse order
    # costs (Nports*Nmodes)^2 transposed solves per frequency and a sparse
    # inner product per component, so it wins once there are more components
    # than output port mode pairs.
    if !(sensitivitymode == :auto || sensitivitymode == :forward ||
            sensitivitymode == :reverse)
        throw(ArgumentError(lazy"sensitivitymode must be :auto, :forward or :reverse, not $(sensitivitymode)."))
    end
    useoperatingpoint = returnSsensitivity &&
        !(isnothing(sensitivitynodeflux) && isnothing(sensitivityresidual))
    if useoperatingpoint && (isnothing(nonlinear) ||
            isnothing(nonlinear.operatingpoint))
        throw(ArgumentError("Including the operating point shift in the sensitivities requires a nonlinear solution with an operating point. Call hbnlsolve with returnoperatingpoint = true."))
    end
    # Validate the low level operating point inputs at the public boundary,
    # before anything sizes itself from them: the contractions index the
    # output arrays (sized by sensitivitynames) and the transposed solutions
    # (sized by the operating point) from these matrices inside @inbounds
    # loops, so malformed inputs must be rejected here, not discovered as
    # memory corruption. The two contraction orders read different inputs
    # (forward the node flux derivatives, reverse the residual derivatives),
    # so supplying both would let inconsistent inputs give order dependent
    # results: reject that too.
    if useoperatingpoint
        op = nonlinear.operatingpoint
        if !isnothing(sensitivitynodeflux) && !isnothing(sensitivityresidual)
            throw(ArgumentError("Supply either sensitivitynodeflux or sensitivityresidual, not both: the forward contraction order reads the former and the reverse order the latter, so inconsistent inputs would make the result depend on the contraction order."))
        end
        if !isnothing(sensitivitynodeflux) &&
                size(sensitivitynodeflux) != (length(op.x), length(sensitivitynames))
            throw(DimensionMismatch(lazy"sensitivitynodeflux must be (operating point state length) x (number of sensitivity components) = ($(length(op.x)), $(length(sensitivitynames))), got $(size(sensitivitynodeflux))."))
        end
        if !isnothing(sensitivityresidual) &&
                size(sensitivityresidual) != (size(op.jacobian, 1), length(sensitivitynames))
            throw(DimensionMismatch(lazy"sensitivityresidual must be (real Jacobian rows) x (number of sensitivity components) = ($(size(op.jacobian, 1)), $(length(sensitivitynames))), got $(size(sensitivityresidual))."))
        end
    end
    # Without Josephson junctions the linearized system matrix does not
    # depend on the operating point, so the contribution is identically
    # zero: drop the operating point inputs rather than build the junction
    # shaped contractions (whose reverse constructor indexes the first
    # junction) on an empty junction set.
    if useoperatingpoint && isempty(nonlinear.operatingpoint.sys.Ljb.nzind)
        useoperatingpoint = false
    end
    usereverse = if !useoperatingpoint
        false
    elseif sensitivitymode == :forward
        false
    elseif sensitivitymode == :reverse
        isnothing(sensitivityresidual) && throw(ArgumentError("The reverse mode sensitivity contraction requires the residual derivatives of calcresidualsensitivity."))
        true
    else # :auto
        # The forward order costs one product against the sparsity structure
        # of the linearized system per component; the reverse order costs one
        # transposed solve through the pump Jacobian per output port mode
        # pair, plus a transform, whatever the component count. So the
        # crossover is at a fixed multiple of the number of output port mode
        # pairs, and only the multiple has to be measured.
        #
        # A solve is not a product, and counting them as if they were put the
        # crossover an order of magnitude early, which cost a factor of ten
        # to anyone who crossed it: on a two hundred cell line at nine modes,
        # a hundred and twenty eight components took 0.42 s in the forward
        # order and 4.17 s in the reverse one. Measured on travelling wave
        # amplifiers of eighty and two hundred cells at five and nine modes,
        # the crossover was 6 to 13 times the pair count, so it is taken as
        # eight. Near a crossover the choice costs little either way, which
        # is what makes one constant serviceable across that spread.
        Ncomponents = isnothing(sensitivityresidual) ?
            size(sensitivitynodeflux, 2) : size(sensitivityresidual, 2)
        !isnothing(sensitivityresidual) && Ncomponents >
            8*(length(signalnm.portindices)*Nsignalmodes)^2
    end

    sensitivitydAop = if useoperatingpoint && !usereverse
        # the forward order contracts against the operating point shift
        # itself, so compute it from the residual derivatives if it was not
        # supplied. the reverse order works from the residual derivatives
        # directly and skips these per-component solves entirely.
        dx = if isnothing(sensitivitynodeflux)
            calcnodefluxsensitivity(nonlinear.operatingpoint,
                sensitivityresidual; factorization = factorization)
        else
            sensitivitynodeflux
        end
        calcoperatingpointstamps(nonlinear.operatingpoint, lsys, dx)
    else
        Vector{Complex{Float64}}[]
    end

    sensitivityreverse = if usereverse
        ReverseSensitivity(nonlinear.operatingpoint, lsys,
            sensitivityresidual)
    else
        nothing
    end


    # use this for debugging purposes to return the linearized system along
    # with the ingredients from which its per-frequency matrices and right
    # hand sides are built, so reference implementations can be constructed,
    # eg. in the tests. mirrors the debugJacobian keyword of hbnlsolve.
    if debuglsys
        return (lsys=lsys, bnm=bnm, wpumpmodes=wpumpmodes,
            phimatrix=phimatrix, Ljb=signalnm.Ljb, Nnodes=Nnodes,
            Nmodes=Nsignalmodes, Nnodalmna=Nnodalmna, Nauxmna=Nauxmna,
            Nauxmnar=Nauxmnar, mnaindices=mnaindices,
            coupledbranches=coupledbranches, vvn=vvn,
            portimpedanceindices=portimpedanceindices,
            noiseportimpedanceindices=noiseportimpedanceindices,
            nodeindices=nodeindices, componenttypes=componenttypes,
            symfreqvar=symfreqvar)
    end


    # make the output arrays for the scattering parameters, noise,
    # sensitivities, quantum efficiency, commutation relations, and node
    # quantities. an output which was not requested is a zero size array,
    # which is the signal that it should not be computed. the scattering
    # parameter sensitivities are scaled by the input waves and their port
    # impedance term depends on S itself, so S is computed whenever the
    # sensitivities are requested even if it is not returned.
    outputarrays = LinearizedArrays(;
        requestS = returnS,
        requestSnoise = returnSnoise,
        requestCnoise = returnCnoise,
        requestSsensitivity = returnSsensitivity,
        requestQE = returnQE, requestCM = returnCM,
        requestnodeflux = returnnodeflux,
        requestnodefluxadjoint = returnnodefluxadjoint,
        requestvoltage = returnvoltage,
        requestvoltageadjoint = returnvoltageadjoint,
        Nports = Nports, Nmodes = Nsignalmodes,
        Nnoisechannels = Nnoisechannels,
        Ncomponents = length(sensitivitynames), Nnodes = Nnodes,
        Nfrequencies = length(w))

    # generate the mode indices and find the signal index
    signalindex = 1

    # solve the linear system for the specified frequencies. the response for
    # each frequency is independent so it can be done in parallel; however
    # we want to reuse the factorization object and other input arrays. 
    # perform array allocations and factorization "nbatches" times.
    # parallelize using tasks
    # On a backend the sweep is solved there instead: the system matrices of
    # a batch of frequencies share one sparsity pattern, so they are assembled
    # by one kernel and factorized and solved by cuDSS as a uniform batch. The
    # frequency loop and every output it computes are the host ones, with only
    # the solve replaced, so both paths compute their outputs identically.
    #
    # It falls back to the host when the component values depend on the
    # symbolic frequency variable, so the assembly is not a constant quadratic
    # in the signal frequency.
    sensitivitytuple = (stamps = sensitivitystamps, dAop = sensitivitydAop,
        reverse = sensitivityreverse)
    usedevice = !(backend isa CPU) && cansweepondevice(lsys)
    if usedevice
        # what each direction's solutions are read for decides whether the
        # whole solution comes back or only some of its rows. The scattering
        # parameters read the forward solution at the port rows and the
        # adjoint one at the noise port rows; the node flux, voltage and
        # sensitivity outputs read all of both.
        fullforward = !isempty(outputarrays.nodeflux) ||
            !isempty(outputarrays.voltage) ||
            !isempty(outputarrays.Ssensitivity)
        fulladjoint = !isempty(outputarrays.nodefluxadjoint) ||
            !isempty(outputarrays.voltageadjoint) ||
            !isempty(outputarrays.Ssensitivity)
        # The noise scattering parameters are computed where the adjoint
        # solution is, so they ask nothing of the host copy. When they are its
        # only reader, which is the usual case on a lossy line, the adjoint
        # solution is never copied back at all.
        # The noise channels of the dissipative scattering blocks are formed
        # on the backend alongside those of the lumped noise ports whenever
        # the blocks' scattering parameters can be evaluated there, which is
        # the same condition their stamps are evaluated there under. When one
        # cannot be, the whole adjoint solution comes back instead and every
        # channel is formed on the host.
        blocknoiseondevice = isnothing(noiseplan) ||
            candeviceevaluate(lsys.scattering)
        wantsnoise = !isempty(outputarrays.Snoise) ||
            !isempty(outputarrays.QE) || !isempty(outputarrays.CM)
        if wantsnoise && !blocknoiseondevice
            @warn "A scattering block of this circuit has scattering parameters which cannot be evaluated on the backend, so its noise channels are formed on the host and the whole adjoint solution is copied back at every frequency, which is the largest transfer in the sweep. Give the block's data as a constant matrix, as tabulated data, or as a callable of one entry at a time (`ScatteringBlock(f; nports, form = :entry)`) to keep the noise on the backend." maxlog=1
        end
        devicenoiseplan = if wantsnoise && blocknoiseondevice &&
                (!isempty(noiseportimpedanceindices) || !isnothing(noiseplan))
            plandevicenoise(nodeindices, componenttypes,
                noiseportimpedanceindices, noiseportimpedances,
                Nsignalmodes, backend)
        else
            nothing
        end
        deviceblocknoiseplan = isnothing(devicenoiseplan) ? nothing :
            plandeviceblocknoise(lsys.scattering, noiseplan, Nsignalmodes,
                backend)
        adjointspec = if needsadjointsolve(outputarrays,
                noiseportimpedanceindices, noiseplan)
            (full = fulladjoint ||
                    (!isnothing(noiseplan) && isnothing(devicenoiseplan)),
                rows = isnothing(devicenoiseplan) ?
                    portsolutionrows(nodeindices, noiseportimpedanceindices,
                        Nsignalmodes) : Int[])
        else
            nothing
        end
        solutions = devicesolutions(lsys, bnm, w, backend,
            (full = fullforward, rows = portsolutionrows(nodeindices,
                portimpedanceindices, Nsignalmodes)),
            adjointspec)
        # One pass over the whole sweep, so the frequency loop's buffers are
        # made once; the callbacks solve the batch a frequency belongs to when
        # they first reach it. Splitting the loop across threads was measured
        # and it loses: a batch holds about a dozen frequencies, so each
        # thread would remake those buffers for a frequency or two of work.
        # The device solves a batch of frequencies, then the host computes
        # their outputs, reusing one set of scratch across every batch.
        #
        # Spreading those outputs across workers was tried and it loses. The
        # host work here is about a quarter of the call, not the half a
        # measurement which also counted the plans, the cuDSS analysis and the
        # output arrays had suggested; and giving each worker its own scratch
        # costs more in working set than the parallelism wins, because the
        # solution buffer alone is the size of the state of the whole circuit.
        # Measured on a thirteen mode travelling wave amplifier, the outputs
        # took 0.132 s on one worker and 0.189 s on eight.
        # How much host work a frequency carries decides how many workers
        # are worth their working set, because a worker's solution buffer is
        # the size of the state of the whole circuit. Without the
        # sensitivities the outputs are cheap and one worker wins: measured
        # on a thirteen mode travelling wave amplifier, they took 0.132 s on
        # one worker and 0.189 s on eight. The sensitivities are not cheap,
        # costing a product against the sparsity structure of the system
        # matrix per component and per frequency; measured on a two hundred
        # cell line with sixty four components, they took 20.7 ms per
        # component on one worker and 2.9 ms on eight.
        nworkers = isempty(outputarrays.Ssensitivity) ? 1 :
            max(1, min(nbatches, Base.Threads.nthreads()))
        wss = [LinearizedWorkspace(outputarrays, sensitivitytuple, lsys,
            Nports, Nsignalmodes, Nnoisechannels,
            length(wpumpmodes), factorization; assembles = false)
            for _ in 1:nworkers]
        # the noise scattering parameters are reduced on the backend against
        # scratch of their own, so a worker which computes them needs its own
        noisecbs = isnothing(devicenoiseplan) ? nothing :
            [devicenoise(devicenoiseplan,
                (isnothing(deviceblocknoiseplan) || t == 1) ?
                    deviceblocknoiseplan : withfactors(deviceblocknoiseplan),
                solutions.providers, solutions.adjointdevice,
                size(bnm, 2), wpumpmodes, w, !isempty(outputarrays.Snoise),
                channeltemperatures)
             for t in 1:nworkers]
        nb = solutions.batchsize
        inner!(t, batch) = hblinsolve_inner!(wss[t], outputarrays,
            sensitivitytuple, lsys, bnm,
            portindices, portimpedanceindices, noiseportimpedanceindices,
            portimpedances, noiseportimpedances, nodeindices,
            componenttypes, w, wpumpmodes, Nsignalmodes, Nnodes,
            symfreqvar, batch, factorization;
            presolved = solutions.forward,
            presolvedadjoint = solutions.adjoint,
            presolvednoise = isnothing(noisecbs) ? nothing : noisecbs[t],
            noiseplan = noiseplan,
            channeltemperatures = channeltemperatures)
        for lo in 1:nb:length(w)
            hi = min(lo + nb - 1, length(w))
            # only this touches the backend, and it is serial; the outputs of
            # the batch it staged are host work on disjoint frequencies
            solutions.solvebatch!(lo)
            if nworkers == 1
                inner!(1, lo:hi)
            else
                chunks = collect(Base.Iterators.partition(lo:hi,
                    1 + (hi - lo) ÷ nworkers))
                Base.Threads.@sync for (t, batch) in enumerate(chunks)
                    Base.Threads.@spawn inner!(t, batch)
                end
            end
        end
    else
        batches = Base.Iterators.partition(1:length(w),1+(length(w)-1)÷nbatches)
        Threads.@sync for batch in batches
            Base.Threads.@spawn hblinsolve_inner!(
                LinearizedWorkspace(outputarrays, sensitivitytuple, lsys,
                    Nports, Nsignalmodes, Nnoisechannels,
                    length(wpumpmodes), factorization),
                outputarrays, sensitivitytuple,
                lsys, bnm,
                portindices, portimpedanceindices, noiseportimpedanceindices,
                portimpedances, noiseportimpedances, nodeindices, componenttypes,
                w, wpumpmodes, Nsignalmodes, Nnodes, symfreqvar, batch,
                factorization; noiseplan = noiseplan,
                channeltemperatures = channeltemperatures)
        end
    end

    # calculate the `ideal` quantum efficiency based on the gain assuming an
    # ideal two mode amplifier
    QEideal = outputarrays.QEideal

    # turn all of the array outputs into keyed arrays if the
    # keyedarrays = true

    # if keyword argument keyedarrays = true then generate keyed arrays
    # for the scattering parameters
    Sout = if returnS && keyedarrays
        Stokeyed(outputarrays.S, modes, portnumbers, modes, portnumbers, w)
    elseif returnS
        outputarrays.S
    else
	zeros(Complex{Float64},0,0,0)
    end

    # if keyword argument keyedarrays = true then generate keyed arrays
    # for the noise scattering parameters
    Snoiseout = if returnSnoise && keyedarrays
        Snoisetokeyed(outputarrays.Snoise, modes,
            noisechannelnames(componentnames, noiseportimpedanceindices,
                noiseplan, ssys), modes, portnumbers, w)
    else
        outputarrays.Snoise
    end

    # the added noise covariance is indexed by output port mode on both
    # sides, like the scattering matrix
    Cnoiseout = if returnCnoise && keyedarrays
        Stokeyed(outputarrays.Cnoise, modes, portnumbers, modes, portnumbers)
    else
        outputarrays.Cnoise
    end

    # if keyword argument keyedarrays = true then generate keyed arrays
    # for Ssensitivity
    Ssensitivityout = if returnSsensitivity && keyedarrays
        Ssensitivitytokeyed(outputarrays.Ssensitivity, modes, portnumbers, modes,
            portnumbers, sensitivitynames, w)
    else
        outputarrays.Ssensitivity
    end

    # if keyword argument keyedarrays = true then generate keyed arrays
    # for the quantum efficiency
    QEout = if returnQE && keyedarrays
        Stokeyed(outputarrays.QE, modes, portnumbers, modes, portnumbers, w)
    else
        outputarrays.QE
    end

    QEidealout = if returnQE && keyedarrays
        Stokeyed(QEideal, modes, portnumbers, modes, portnumbers, w)
    else
        QEideal
    end

    # if keyword argument keyedarrays = true then generate keyed arrays
    # for the commutation relations
    CMout = if returnCM && keyedarrays
        CMtokeyed(outputarrays.CM, modes, portnumbers, w)
    else
        outputarrays.CM
    end

    # if keyword argument keyedarrays = true then generate keyed arrays
    nodefluxout = if returnnodeflux && keyedarrays
        nodevariabletokeyed(outputarrays.nodeflux, modes, nodenames, modes,
            portnumbers, w)
    else
        outputarrays.nodeflux
    end

    # if keyword argument keyedarrays = true then generate keyed arrays
    nodefluxadjointout = if returnnodefluxadjoint && keyedarrays
        nodevariabletokeyed(outputarrays.nodefluxadjoint, modes,
            nodenames, modes, portnumbers, w)
    else
        outputarrays.nodefluxadjoint
    end

    # if keyword argument keyedarrays = true then generate keyed arrays
    voltageout = if returnvoltage && keyedarrays
        nodevariabletokeyed(outputarrays.voltage, modes,
            nodenames, modes, portnumbers, w)
    else
        outputarrays.voltage
    end

    # if keyword argument keyedarrays = true then generate keyed arrays
    voltageadjointout = if returnvoltageadjoint && keyedarrays
        nodevariabletokeyed(outputarrays.voltageadjoint, modes,
            nodenames, modes, portnumbers, w)
    else
        outputarrays.voltageadjoint
    end

    return LinearizedHB(w, modes, Sout, Snoiseout, Cnoiseout, Ssensitivityout,
        QEout,
        QEidealout, CMout, nodefluxout, nodefluxadjointout, voltageout,
        voltageadjointout, nodenames, nodeindices, componentnames,
        componenttypes, componentnamedict, mutualinductorbranchnames,
        portnumbers, portindices, portimpedanceindices,
        noiseportimpedanceindices, sensitivitynames, sensitivityindices,
        Nsignalmodes, Nnodes, Nbranches, Nports, signalindex)
end

"""
    LinearizedArrays(; requestS, requestSnoise, requestSsensitivity,
        requestQE, requestCM, requestnodeflux, requestnodefluxadjoint,
        requestvoltage, requestvoltageadjoint, Nports, Nmodes,
        Nnoisechannels, Ncomponents, Nnodes, Nfrequencies)

The preallocated output arrays of one [`hblinsolve`](@ref) run, filled per
frequency by [`hblinsolve_inner!`](@ref). An output which was not requested
is a zero size array of the same dimensionality, which is the signal, via
`isempty`, that it should not be computed. This is the single place the
output shapes and the request conditions are defined.
"""
struct LinearizedArrays
    S::Array{Complex{Float64},3}
    Snoise::Array{Complex{Float64},3}
    # the added noise covariance at the output ports, `Y` of the Gaussian
    # channel whose `X` is `S`
    Cnoise::Array{Complex{Float64},3}
    Ssensitivity::Array{Complex{Float64},4}
    QE::Array{Float64,3}
    QEideal::Array{Float64,3}
    CM::Array{Float64,2}
    nodeflux::Array{Complex{Float64},3}
    nodefluxadjoint::Array{Complex{Float64},3}
    voltage::Array{Complex{Float64},3}
    voltageadjoint::Array{Complex{Float64},3}
end

function LinearizedArrays(; requestS::Bool, requestSnoise::Bool,
    requestCnoise::Bool = false,
    requestSsensitivity::Bool, requestQE::Bool, requestCM::Bool,
    requestnodeflux::Bool, requestnodefluxadjoint::Bool,
    requestvoltage::Bool, requestvoltageadjoint::Bool, Nports::Integer,
    Nmodes::Integer, Nnoisechannels::Integer, Ncomponents::Integer,
    Nnodes::Integer, Nfrequencies::Integer)

    NPM = Nports*Nmodes
    Nnodal = Nmodes*(Nnodes-1)
    za(request, T, dims...) = request ? zeros(T, dims...) :
        zeros(T, ntuple(_ -> 0, length(dims))...)
    return LinearizedArrays(
        za(requestS, Complex{Float64}, NPM, NPM, Nfrequencies),
        za(requestSnoise, Complex{Float64}, Nnoisechannels*Nmodes, NPM,
            Nfrequencies),
        za(requestCnoise, Complex{Float64}, NPM, NPM, Nfrequencies),
        za(requestSsensitivity, Complex{Float64}, NPM, NPM, Ncomponents,
            Nfrequencies),
        za(requestQE, Float64, NPM, NPM, Nfrequencies),
        za(requestQE, Float64, NPM, NPM, Nfrequencies),
        za(requestCM, Float64, NPM, Nfrequencies),
        za(requestnodeflux, Complex{Float64}, Nnodal, NPM, Nfrequencies),
        za(requestnodefluxadjoint, Complex{Float64}, Nnodal, NPM,
            Nfrequencies),
        za(requestvoltage, Complex{Float64}, Nnodal, NPM, Nfrequencies),
        za(requestvoltageadjoint, Complex{Float64}, Nnodal, NPM,
            Nfrequencies))
end

"""
    LinearizedWorkspace

The scratch of one worker's pass over a range of signal frequencies in
[`hblinsolve_inner!`](@ref).

There are fourteen buffers here, the largest of them the size of the solution
and of the system matrix, so making them is not free. They are a separate
object so that a worker can be given one once and reuse it across every range
it is handed, rather than remaking them per call. That is what lets the device
path hand out a batch of frequencies at a time and still keep the host work
parallel: without it, each worker would remake all of this for the frequency
or two of a batch that fell to it, which was measured and cost more than the
parallelism won.
"""
struct LinearizedWorkspace{TA,TC,TR}
    phin::Matrix{Complex{Float64}}
    inputwave::Matrix{Complex{Float64}}
    outputwave::Matrix{Complex{Float64}}
    noiseoutputwave::Matrix{Complex{Float64}}
    phinforward::Matrix{Complex{Float64}}
    dAsparse::TA
    dAphin::Matrix{Complex{Float64}}
    sensitivitycontraction::Matrix{Complex{Float64}}
    sensitivitycache::TC
    sensitivityrevbufs::TR
    sensitivitygamma::Vector{Complex{Float64}}
    sensitivitybeta::Vector{Complex{Float64}}
    wmodes::Vector{Float64}
    Asparsecopy::SparseMatrixCSC{Complex{Float64},Int}
    cache::FactorizationCache
    Sworking::Matrix{Complex{Float64}}
    Snoiseworking::Matrix{Complex{Float64}}
    # the scratch of the scattering block vacuum noise channels
    scatteringnoisework::ScatteringNoiseWorkspace
    # the occupation of each noise channel mode, rebuilt per frequency
    occupation::Vector{Float64}
end

"""
    LinearizedWorkspace(arrays::LinearizedArrays, sensitivity, lsys, Nports,
        Nmodes, Nnoisechannels, Nwpumpmodes, factorization;
        assembles::Bool = true)

Make the scratch of one worker. `assembles` is false when the solutions are
supplied from elsewhere and nothing writes into a copy of the system matrix,
which on a large problem is the biggest allocation here.
"""
function LinearizedWorkspace(arrays::LinearizedArrays, sensitivity, lsys,
    Nports::Integer, Nmodes::Integer, Nnoisechannels::Integer,
    Nwpumpmodes::Integer, factorization; assembles::Bool = true)

    n = size(lsys.Asparse, 1)
    np = Nports*Nmodes
    cplx(a, b) = zeros(Complex{Float64}, a, b)
    wantssensitivity = !isempty(arrays.Ssensitivity)
    # the forward solution, which the adjoint solve overwrites, the derivative
    # of the system matrix with respect to one component, and the contraction.
    # The operating point contribution is dense on the sparsity structure of
    # the system matrix, so it needs a matrix and a product; the stamps of the
    # individual components do not.
    wantsop = wantssensitivity && !isempty(sensitivity.dAop)
    phinforward = wantssensitivity ? cplx(n, np) : cplx(0, 0)
    # one factorization of the pump Jacobian per worker: a sparse
    # factorization is not safe to solve against from several threads at once
    sensitivitycache = if isnothing(sensitivity.reverse)
        nothing
    else
        c = FactorizationCache()
        tryfactorize!(c, factorization, sensitivity.reverse.op.jacobian)
        c
    end
    return LinearizedWorkspace(
        cplx(n, np), cplx(np, np), cplx(np, np),
        cplx(Nnoisechannels*Nmodes, np),
        phinforward,
        wantsop ? copy(lsys.Asparse) : lsys.Asparse,
        wantsop ? cplx(size(phinforward)...) : cplx(0, 0),
        cplx(np, np), sensitivitycache,
        isnothing(sensitivity.reverse) ? nothing :
            ReverseSensitivityBuffers(sensitivity.reverse, np),
        zeros(Complex{Float64}, np), zeros(Complex{Float64}, np),
        zeros(Float64, Nwpumpmodes),
        # a copy of the system matrix, because it is modified per frequency,
        # potentially by several workers at once
        assembles ? copy(lsys.Asparse) : lsys.Asparse,
        FactorizationCache(), cplx(np, np), cplx(Nnoisechannels*Nmodes, np),
        ScatteringNoiseWorkspace(),
        zeros(Float64, Nnoisechannels*Nmodes))
end

"""
    hblinsolve_inner!(arrays::LinearizedArrays, sensitivity, lsys, bnm,
        portindices, portimpedanceindices, noiseportimpedanceindices,
        portimpedances, noiseportimpedances, nodeindices, componenttypes,
        w, wpumpmodes, Nmodes, Nnodes, symfreqvar, wi, factorization)

Solve the linearized harmonic balance problem for a subset of the frequencies
given by `wi`, assembling the per-frequency system matrices from the
[`HBLinearizedSystem`](@ref) `lsys` with [`assemblesystemmatrix!`](@ref),
and writing into the output arrays of the [`LinearizedArrays`](@ref)
`arrays`. An empty output array means that output was not requested;
requested outputs are written per frequency through views, and small
working matrices stand in for the outputs which are computed but not stored
(for example `S` when only the quantum efficiency needs it). `sensitivity`
is a named tuple with the fixed operating point `stamps`, the operating
point `dAop` stamps of the forward contraction order, and the
[`ReverseSensitivity`](@ref) of the reverse order (or `nothing`). This
function is thread safe in that different frequencies can be computed in
parallel on separate threads; `lsys`, `arrays` (through disjoint per
frequency views) and `sensitivity` are shared.
"""
function hblinsolve_inner!(ws::LinearizedWorkspace, arrays::LinearizedArrays,
    sensitivity, lsys, bnm,
    portindices, portimpedanceindices, noiseportimpedanceindices,
    portimpedances, noiseportimpedances, nodeindices,
    componenttypes, w, wpumpmodes, Nmodes, Nnodes, symfreqvar, wi, factorization;
    noiseplan = nothing, channeltemperatures = nothing,
    presolved = nothing, presolvedadjoint = nothing,
    presolvednoise = nothing)

    # `presolved` replaces the assemble, factorize and solve of each frequency
    # with a callback which fills the solution some other way, and
    # `presolvedadjoint` likewise replaces the transposed solve. That is how
    # the device sweep hands back solutions it computed a batch at a time (see
    # [`devicesolutions`](@ref)). Everything downstream of the solve is
    # unchanged, so the two paths compute their outputs identically.
    isnothing(presolved) || !isnothing(presolvedadjoint) ||
        !needsadjointsolve(arrays, noiseportimpedanceindices, noiseplan) ||
        throw(ArgumentError(
            "a supplied solution needs a supplied adjoint solution too: the transposed solve otherwise needs the factorization of the forward system, which a supplied solution does not carry."))

    Nports = length(portindices)
    Nnoiseports = length(noiseportimpedanceindices)
    phin = ws.phin
    inputwave = ws.inputwave
    outputwave = ws.outputwave
    noiseoutputwave = ws.noiseoutputwave
    phinforward = ws.phinforward
    dAsparse = ws.dAsparse
    dAphin = ws.dAphin
    sensitivitycontraction = ws.sensitivitycontraction
    sensitivitycache = ws.sensitivitycache
    sensitivityrevbufs = ws.sensitivityrevbufs
    sensitivitygamma = ws.sensitivitygamma
    sensitivitybeta = ws.sensitivitybeta
    wmodes = ws.wmodes
    Asparsecopy = ws.Asparsecopy
    cache = ws.cache
    Sworking = ws.Sworking
    Snoiseworking = ws.Snoiseworking
    scatteringnoisework = ws.scatteringnoisework
    occupation = ws.occupation

    # whether the transposed (adjoint) system must be solved at each
    # frequency: for the sensitivities always, and otherwise when a
    # consumer of the adjoint solution (the noise scattering parameters,
    # the quantum efficiency, the commutation relations, or the adjoint
    # node outputs) is requested together with a source of it. This does
    # not depend on the frequency.
    needsadjoint = needsadjointsolve(arrays, noiseportimpedanceindices,
        noiseplan)

    # loop over the frequencies
    for i in wi

        Sview = isempty(arrays.S) ? Sworking : view(arrays.S, :, :, i)
        Snoiseview = isempty(arrays.Snoise) ? Snoiseworking : view(arrays.Snoise, :, :, i)

        # calculate the mode frequencies
        ws = w[i]
        wmodes .= ws .+ wpumpmodes

        # assemble the linearized system matrix at this frequency,
        # Asparsecopy = (AoLjnm + invLnm + im.*Gnm.*w - Cnm.*w.^2) with the
        # per column mode frequency, in a way that doesn't allocate
        # significant memory, taking the complex conjugates of the negative
        # frequency mode entries of the linear term matrices and
        # substituting any symbolic frequency variables.
        if isnothing(presolved)
            assemblesystemmatrix!(Asparsecopy, lsys, wmodes)

            # factor the sparse matrix
            # factorklu!(cache, Asparsecopy)
            tryfactorize!(cache, factorization, Asparsecopy)

            # solve the linear system
            trysolve!(phin, cache.factorization, bnm)
        else
            presolved(i, phin)
        end

        # convert to node voltages. node flux is defined as the time integral
        # of node voltage so node voltage is derivative of node flux which can
        # be accomplished in the frequency domain by multiplying by j*w.
        if !isempty(arrays.voltage)
            vv = view(arrays.voltage, :, :, i)
            @inbounds for t in axes(vv, 1)
                wm = im*wmodes[(t-1) % Nmodes + 1]
                for kc in axes(vv, 2)
                    vv[t, kc] = wm*phin[t, kc]
                end
            end
        end

        # copy the nodeflux for output. the auxiliary variables of the
        # modified nodal analysis augmentation are internal.
        if !isempty(arrays.nodeflux)
            copy!(view(arrays.nodeflux,:,:,i), view(phin, 1:size(arrays.nodeflux,1), :))
        end

        # calculate the scattering parameters.
        if !isempty(arrays.S) || !isempty(arrays.QE) || !isempty(arrays.QEideal) ||
                !isempty(arrays.CM) || !isempty(arrays.Ssensitivity)
            calcinputoutput!(inputwave, outputwave, phin, bnm,
                portimpedanceindices, portimpedanceindices, portimpedances,
                portimpedances, nodeindices, componenttypes, wmodes, symfreqvar)
            calcscatteringmatrix!(Sview, inputwave, outputwave)

            # the scalars which convert the adjoint contraction into the
            # scattering parameter derivatives read the input waves computed
            # just above, so they are formed here, before the Z parameter and
            # adjoint output calculations below overwrite inputwave with the
            # differently normalized input currents.
            if !isempty(arrays.Ssensitivity)
                calcsensitivityscaling!(sensitivitygamma, sensitivitybeta,
                    inputwave, bnm, portimpedanceindices, portimpedances,
                    componenttypes, nodeindices, wmodes, Nmodes)
            end
        end

        if needsadjoint

            # retain the forward solution, which the adjoint solve overwrites
            if !isempty(arrays.Ssensitivity)
                copy!(phinforward, phin)
            end

            # Solve the transposed linearized system, reusing the
            # factorization of the forward system. The adjoint solutions the
            # noise and quantum efficiency calculations require are the
            # solutions of the transposed system: by the adjoint identity the
            # response at an output port to a source anywhere in the circuit
            # is that source contracted against the transposed solution
            # driven at the port. For a circuit without scattering blocks it
            # is also the solution of the system with the complex conjugate
            # of the pump modulation matrix, to which it is related by the
            # diagonal similarity transformation of
            # [`assemblesystemmatrix!`](@ref); the hybrid scattering rows
            # break that similarity, so with blocks the two are different
            # systems and it is the transposed one, whose auxiliary port
            # current rows the block noise channels read, that is meant.
            if isnothing(presolvedadjoint)
                trysolvetranspose!(phin, cache.factorization, bnm)
            else
                presolvedadjoint(i, phin)
            end

            # copy the nodeflux adjoint for output
            if !isempty(arrays.nodefluxadjoint)
                copy!(view(arrays.nodefluxadjoint,:,:,i), view(phin, 1:size(arrays.nodefluxadjoint,1), :))
            end

            if !isempty(arrays.voltageadjoint)
                vv = view(arrays.voltageadjoint, :, :, i)
                @inbounds for t in axes(vv, 1)
                    wm = im*wmodes[(t-1) % Nmodes + 1]
                    for kc in axes(vv, 2)
                        vv[t, kc] = wm*phin[t, kc]
                    end
                end
            end

            # calculate the noise scattering parameters. `presolvednoise`
            # computes them where the adjoint solution was computed and
            # returns only what the quantum efficiency and the commutation
            # relations read of them, which is a vector rather than a matrix
            # with a row per noise port mode.
            noiseterm = transpose(Snoiseview)
            wantscnoise = !isempty(arrays.Cnoise)
            if !isempty(arrays.Snoise) || wantscnoise ||
                    !isempty(arrays.QE) || !isempty(arrays.CM)
                if isnothing(presolvednoise)
                    calcinputoutputnoise!(inputwave, noiseoutputwave, phin, bnm,
                        portimpedanceindices, noiseportimpedanceindices,
                        portimpedances, noiseportimpedances, nodeindices,
                        componenttypes, wmodes, symfreqvar)
                    # the vacuum noise the dissipative scattering blocks add,
                    # in channels which follow the lumped noise ports
                    if !isnothing(noiseplan)
                        scatteringnoisewaves!(noiseoutputwave, noiseplan,
                            lsys.scattering, phin, wmodes,
                            Nnoiseports*Nmodes, scatteringnoisework)
                    end
                    calcscatteringmatrix!(Snoiseview, inputwave, noiseoutputwave)
                else
                    noiseterm = presolvednoise(i, inputwave, Snoiseview,
                        wantscnoise ? view(arrays.Cnoise,:,:,i) : nothing)
                end
            end

            # calculate the scattering parameter sensitivities. phin now
            # holds the solution of the transposed system, which is the
            # adjoint solution the contraction needs.
            if !isempty(arrays.Ssensitivity)
                calcSsensitivity!(view(arrays.Ssensitivity,:,:,:,i),
                    sensitivity.stamps, sensitivity.dAop, dAsparse, dAphin,
                    phinforward, phin, Sview, sensitivitygamma,
                    sensitivitybeta, sensitivitycontraction, wmodes,
                    Nmodes, symfreqvar)
                if !isnothing(sensitivity.reverse)
                    calcSsensitivityreverse!(view(arrays.Ssensitivity,:,:,:,i),
                        sensitivity.reverse, lsys, phinforward, phin,
                        sensitivitygamma, sensitivitybeta,
                        sensitivitycache, sensitivityrevbufs)
                end
            end

            # calculate the quantum efficiency
            if !isempty(arrays.QE) || wantscnoise
                # the occupation of each channel, which is where a
                # temperature enters: the noise scattering parameters are a
                # scattering matrix and carry none
                noiseoccupation!(occupation, channeltemperatures, wmodes,
                    Nmodes)
            end
            if wantscnoise && isnothing(presolvednoise)
                calcnoisecovariance!(view(arrays.Cnoise,:,:,i), Snoiseview,
                    occupation)
            end
            if !isempty(arrays.QE)
                calcqe!(view(arrays.QE,:,:,i), Sview, noiseterm, occupation)
            end

            # calculate the commutation relations (Manley-Rowe relations)
            if !isempty(arrays.CM)
                calccm!(view(arrays.CM,:,i), Sview, noiseterm, wmodes)
            end
        else
            # calculate the quantum efficiency
            if !isempty(arrays.QE)
                calcqe!(view(arrays.QE,:,:,i), Sview)
            end

            # calculate the commutation relations (Manley-Rowe relations)
            if !isempty(arrays.CM)
                calccm!(view(arrays.CM,:,i), Sview, wmodes)
            end
        end

        # calculate the ideal QE
        if !isempty(arrays.QEideal)
            calcqeideal!(view(arrays.QEideal,:,:,i), Sview)
        end
    end
    return nothing
end

"""
    HBOperatingPoint(sys, x, jacobian, modelayout, Nnodal, Lmean, wmodes,
        Amna, mnaindices, coupledbranches, Nmodes, Nnodes)

The converged pump operating point of [`hbnlsolve`](@ref) together with
everything needed to propagate a component perturbation through it: the
[`HBSystem`](@ref) evaluation object, the converged augmented state, the
exact Jacobian of the equivalent real system assembled there, and the
scaled matrices and layout of the augmented system.

Requested with `returnoperatingpoint = true`. The Jacobian is the exact
Jacobian of the equivalent real system, assembled with
[`assemblerealjacobian!`](@ref), rather than the complex holomorphic
Jacobian of the `:quasinewton` method, which is only an approximation: the
harmonic balance residual is not complex differentiable, so the implicit
function theorem does not hold with the holomorphic Jacobian, while in the
real representation it applies directly.
"""
struct HBOperatingPoint
    # the HBSystem evaluation object and the ModeLayout; untyped because
    # hbsolve.jl is included before hbsystem.jl and realcomplexconv.jl
    sys
    x::Vector{Complex{Float64}}
    jacobian::SparseMatrixCSC{Float64,Int}
    modelayout
    Nnodal::Int
    Lmean::Complex{Float64}
    wmodes::Vector{Float64}
    Amna::SparseMatrixCSC
    mnaindices::Vector{Int}
    coupledbranches::Vector{Int}
    Nmodes::Int
    Nnodes::Int
end

"""
    componentlookups(mnaindices, coupledbranches, Ljb)

Constant time lookups for [`componentstamp`](@ref), built once per stamp
table rather than searched per component: the ordinal of a promoted
resistor within `mnaindices`, the set of mutually coupled branches, and the
ordinal of a junction branch within `Ljb.nzind`. Without these the
classification repeats linear searches per component, which becomes
quadratic over a large sensitivity set.
"""
function componentlookups(mnaindices, coupledbranches, Ljb)
    return (mnaordinal = Dict(idx => r for (r, idx) in enumerate(mnaindices)),
        coupled = Set(coupledbranches),
        junctionordinal = Dict(b => j for (j, b) in enumerate(Ljb.nzind)))
end

"""
    componentstamp(idx::Integer, psc::ParsedSortedCircuit, cg::CircuitGraph,
        nm::CircuitMatrices, lookups, Nmodes::Integer, Nnodes::Integer)

Classify the component at index `idx` for sensitivity analysis and build its
raw one-component matrix, without any solver scaling, negative frequency
conjugation, or padding, which the callers apply for their own grids. The
component matrices are built with the same functions which build the system
matrices, [`calcCn`](@ref), [`calcGn`](@ref), [`calcLb`](@ref) and
[`calcinvLn`](@ref), applied to the single component, so the node and mode
conventions agree by construction. Returns one of

- `(:C, M)`: the component's capacitance matrix,
- `(:G, M)`: the component's conductance matrix, for a resistor which is not
    promoted by the modified nodal analysis formulation,
- `(:Gpromoted, r)`: the ordinal `r` of a promoted resistor within
    `mnaindices`, whose conductance appears only in its constitutive
    equation rows of the caller's augmentation,
- `(:Lj, j)`: the ordinal `j` of a Josephson junction within the junction
    branch vector `nm.Ljb`,
- `(:invL, M)`: the component's inverse inductance matrix.

This is the single definition of which components are supported: `:C`, `:L`,
`:R` and `:Lj` with numeric values. Mutually coupled inductors and
components with symbolic (frequency dependent) values throw, with the same
message from both the fixed operating point stamps
([`calcsensitivitystamps`](@ref)) and the residual derivatives
([`calcresidualsensitivity`](@ref)). `lookups` are the constant time
tables of [`componentlookups`](@ref).
"""
function componentstamp(idx::Integer, psc::ParsedSortedCircuit,
    cg::CircuitGraph, nm::CircuitMatrices, lookups,
    Nmodes::Integer, Nnodes::Integer)

    componenttypes = psc.componenttypes
    nodeindices = psc.nodeindices
    vvn = nm.vvn
    componenttype = componenttypes[idx]
    value = vvn[idx]
    if !(value isa Number)
        throw(ArgumentError(lazy"Sensitivities require a numeric component value, but the value of $(psc.componentnames[idx]) is $(value). Components with symbolic frequency dependent values are not supported."))
    end
    n1 = nodeindices[1, idx]
    n2 = nodeindices[2, idx]
    if componenttype == :C
        return (:C, calcCn(componenttypes[[idx]], nodeindices[:,[idx]],
            vvn[[idx]], Nmodes, Nnodes))
    elseif componenttype == :R
        r = get(lookups.mnaordinal, idx, nothing)
        if isnothing(r)
            return (:G, calcGn(componenttypes[[idx]], nodeindices[:,[idx]],
                vvn[[idx]], Nmodes, Nnodes))
        else
            return (:Gpromoted, r)
        end
    elseif componenttype == :L
        b = cg.edge2indexdict[(n1, n2)]
        if b in lookups.coupled
            throw(ArgumentError(lazy"Sensitivities are not supported for the mutually coupled inductor $(psc.componentnames[idx])."))
        end
        Lb = calcLb(componenttypes[[idx]], nodeindices[:,[idx]],
            vvn[[idx]], cg.edge2indexdict, 1, cg.Nbranches)
        return (:invL, calcinvLn(Lb, cg.Rbn, Nmodes))
    elseif componenttype == :Lj
        b = cg.edge2indexdict[(n1, n2)]
        j = get(lookups.junctionordinal, b, nothing)
        if isnothing(j)
            throw(ArgumentError(lazy"The Josephson junction $(psc.componentnames[idx]) was not found in the branch inductance vector."))
        end
        return (:Lj, j)
    else
        throw(ArgumentError(lazy"Sensitivities are only supported for C, L, R, and Lj components, not $(componenttype), the type of $(psc.componentnames[idx])."))
    end
end

"""
    calcresidualsensitivity(op::HBOperatingPoint, psc, cg, nm,
        sensitivityindices)

Calculate the derivative of the harmonic balance residual with respect to a
relative (logarithmic) perturbation of each component value, at the
operating point. Combined with the implicit function theorem applied to
`F(x, r) = 0` in the equivalent real representation,

    dx/dr = -inv(J)*(dF/dr),

with `J` the exact real Jacobian retained by [`hbnlsolve`](@ref), this gives
the derivative of the operating point itself
(see [`calcnodefluxsensitivity`](@ref)). Returns a sparse matrix whose
columns are `dF/dr` for each component, in the real representation of the
augmented residual: each component touches only its own rows (its nodes and
modes, a promoted resistor's constitutive rows, or a junction branch's
Kirchhoff rows), so the storage scales with the touched entries rather than
with `Nstate*Ncomponents`.

The residual is affine in `C`, `1/R` and `1/L`, so those parameter
derivatives are that component's own contribution to the linear term applied
to the converged state, with a sign, built from the shared classification of
[`componentstamp`](@ref) with the same nondimensionalization and negative
frequency conjugation the solver applies. The Josephson junction term is
the residual's own sine contribution restricted to that junction.

Note that the auxiliary branch currents of the modified nodal analysis
formulation are scaled by the solver scale (see [`calcsolverscale`](@ref)),
which itself depends on the port impedances, so for a port resistor the
auxiliary rows of the returned derivative differ from a finite difference of
a re-solve by that change of normalization. The node flux rows, which are
the physical quantity and the only rows the linearized system depends on,
are unaffected.
"""
function calcresidualsensitivity(op::HBOperatingPoint,
    psc::ParsedSortedCircuit, cg::CircuitGraph, nm::CircuitMatrices,
    sensitivityindices)

    if isnothing(op.jacobian)
        throw(ArgumentError("The operating point does not contain a Jacobian."))
    end
    # evaluate at the operating point regardless of what the shared
    # evaluation object was last used for: the Josephson term below reads
    # the cached time domain branch fluxes.
    setpoint!(op.sys, op.x)
    Ntot = length(op.x)
    Nmodes = op.Nmodes
    Nnodes = op.Nnodes
    wmodes = op.wmodes

    # the component's own matrix, nondimensionalized, negative frequency
    # conjugated and padded exactly as hbnlsolve does with the full matrices
    function scaledpadded(M)
        Ms = SparseMatrixCSC{Complex{Float64},Int}(copy(M))
        conjnegfreq!(Ms, wmodes)
        rmul!(Ms, op.Lmean)
        return mnapadto(Ms, Ntot)
    end

    # The residual derivatives are sparse by nature: each component touches
    # only its own rows of the state (its nodes and modes, one promoted
    # resistor's constitutive rows, or one junction branch's Kirchhoff
    # rows), so they are accumulated as triplets of the real representation
    # and returned as a sparse matrix. The dense alternative peaks at
    # O(Nstate*Ncomponents) memory, which is exactly the many component
    # regime the reverse contraction order exists for. Duplicate triplets
    # (a component matrix with several entries in one row) are summed by
    # sparse, and the real representation is linear, so per entry
    # realification distributes over the sums.
    isrealmode = op.modelayout.isreal
    nmd = length(isrealmode)
    # the offset of complex entry r in the real representation (one real
    # row for the self conjugate modes, a real and an imaginary row
    # otherwise), matching complex_to_real! with the default scale.
    realindexmap = zeros(Int, Ntot)
    k = 1
    for j in eachindex(realindexmap)
        realindexmap[j] = k
        k += isrealmode[(j-1) % nmd + 1] ? 1 : 2
    end
    Nreal = realdim(Ntot, isrealmode)
    Ir = Int[]; Jc = Int[]; Vr = Float64[]
    function pushentry!(r, comp, v)
        kr = realindexmap[r]
        push!(Ir, kr); push!(Jc, comp); push!(Vr, real(v))
        if !isrealmode[(r-1) % nmd + 1]
            push!(Ir, kr+1); push!(Jc, comp); push!(Vr, imag(v))
        end
        return nothing
    end

    # work vector and branch support for the Josephson terms, which are
    # evaluated through the full residual machinery and then read off on
    # the Kirchhoff rows of that junction's branch only.
    cwork = Complex{Float64}[]
    branchsupport = Vector{Tuple{Int,Int}}[]

    x = op.x
    Nm = length(wmodes)
    lookups = componentlookups(op.mnaindices, op.coupledbranches,
        op.sys.Ljb)
    for (comp, idx) in enumerate(sensitivityindices)
        kind, info = componentstamp(idx, psc, cg, nm, lookups,
            Nmodes, Nnodes)
        if kind == :C || kind == :G || kind == :invL
            # dF_comp = c * Ms * Diagonal(w.^power) * x, accumulated per
            # stored entry of the component's own (tiny) matrix.
            Ms = scaledpadded(info)
            c, power = kind == :C ? (-1.0 + 0im, 2) :
                kind == :G ? (0.0 - 1im, 1) : (-1.0 + 0im, 0)
            rows = rowvals(Ms)
            vals = nonzeros(Ms)
            for col in axes(Ms, 2)
                xc = x[col]
                iszero(xc) && continue
                w = wmodes[(col-1) % Nm + 1]
                scale = power == 0 ? c*xc : power == 1 ? c*w*xc : c*w^2*xc
                for pp in nzrange(Ms, col)
                    pushentry!(rows[pp], comp, vals[pp]*scale)
                end
            end
        elseif kind == :Gpromoted
            # a promoted resistor appears only in the conductance entries
            # of its constitutive equations, which are the auxiliary rows
            # of that resistor with node flux columns.
            rowlo = op.Nnodal + (info-1)*Nmodes + 1
            rowhi = op.Nnodal + info*Nmodes
            rows = rowvals(op.Amna)
            vals = nonzeros(op.Amna)
            for col in 1:op.Nnodal
                xc = x[col]
                iszero(xc) && continue
                for pp in nzrange(op.Amna, col)
                    r = rows[pp]
                    if rowlo <= r <= rowhi
                        pushentry!(r, comp, -vals[pp]*xc)
                    end
                end
            end
        else # :Lj
            if isempty(cwork)
                cwork = zeros(Complex{Float64}, Ntot)
                branchsupport = branchnodesandsigns(op.sys.Rbnm, Nmodes,
                    size(op.sys.Rbnm, 1) ÷ Nmodes)
            end
            residualjosephsonterm!(cwork, op.sys, info)
            branch = op.sys.Ljb.nzind[info]
            for (node, _) in branchsupport[branch]
                for m in 1:Nmodes
                    r = (node-1)*Nmodes + m
                    pushentry!(r, comp, -cwork[r])
                end
            end
        end
    end
    return sparse(Ir, Jc, Vr, Nreal, length(sensitivityindices))
end

"""
    calcnodefluxsensitivity(op::HBOperatingPoint, dFr::AbstractMatrix;
        factorization = KLUfactorization())

Solve `dx/dr = -inv(J)*(dF/dr)` for the residual sensitivities `dFr` of
[`calcresidualsensitivity`](@ref), with one factorization of the exact real
Jacobian of the operating point. Returns a matrix whose columns are `dx/dr`
for each component, in the complex representation of the augmented state.
"""
function calcnodefluxsensitivity(op::HBOperatingPoint, dFr::AbstractMatrix;
    factorization = KLUfactorization())

    Ntot = length(op.x)
    cache = FactorizationCache()
    tryfactorize!(cache, factorization, op.jacobian)
    # dFr is sparse from construction; the factorization wants a dense
    # right hand side, so densify one column at a time rather than the
    # whole state by component matrix.
    rhs = zeros(Float64, size(dFr, 1))
    dxr = zeros(Float64, size(dFr, 1))
    dx = zeros(Complex{Float64}, Ntot, size(dFr, 2))
    for k in axes(dx, 2)
        rhs .= view(dFr, :, k)
        trysolve!(dxr, cache.factorization, rhs)
        rmul!(dxr, -1)
        real_to_complex!(view(dx,:,k), dxr, op.modelayout.isreal)
    end
    return dx
end

# the Josephson contribution of the residual, restricted to junction j: the
# node vector of the Fourier coefficients of sin(phi_b(t))/Lj of that
# junction alone, which is the derivative of the residual with respect to
# minus the logarithm of that junction inductance.
function residualjosephsonterm!(out, sys, j::Integer)
    _ensuresin!(sys)
    copyto!(sys.worktd, sys.sintd)
    applyfft!(sys.phimatrix, sys.worktd, sys.rfftplan)
    for i in axes(sys.phimatrix, ndims(sys.phimatrix))
        if i != j
            selectdim(sys.phimatrix, ndims(sys.phimatrix), i) .= 0
        end
    end
    # the Josephson contribution alone, without the linear term
    applybackwardterm!(out, sys.nonlineartermplan, sys.phimatrix, sys.x;
        addlinearterm = false)
    return out
end

# the entries of M in the given row and column ranges, with every other entry
# dropped, keeping the size of M.
function maskrowscolumns(M::SparseMatrixCSC, firstrow::Integer,
    lastrow::Integer, firstcol::Integer, lastcol::Integer)
    I, J, V = findnz(M)
    keep = (firstrow .<= I .<= lastrow) .& (firstcol .<= J .<= lastcol)
    return SparseMatrixCSC{Complex{Float64},Int}(sparse(I[keep], J[keep],
        Complex{Float64}.(V[keep]), size(M,1), size(M,2)))
end




"""
    ReverseSensitivity(op, dFr, T, Em, nzrow, nzcol, realindexmap,
        branchnodes)

Everything the reverse mode contraction of [`calcSsensitivityreverse!`](@ref)
needs, precomputed once and shared read only across the signal frequencies.
"""
struct ReverseSensitivity
    op::HBOperatingPoint
    # the residual derivatives, sparse: each component touches only its own
    # nodes, so the inner product per component is over a handful of entries
    # rather than over the whole state.
    dFr::SparseMatrixCSC{Float64,Int}
    T::Matrix{Float64}
    # the transpose of the frequency to time domain map is a forward
    # transform again (the discrete Fourier matrix is symmetric), executed
    # with this plan by applyffttranspose! against per-thread work arrays.
    fftplan
    nzrow::Vector{Int}
    nzcol::Vector{Int}
    realindexmap::Vector{Int}
    branchnodes
end

"""
    ReverseSensitivityBuffers(rev::ReverseSensitivity)

The mutable work arrays of one invocation of
[`calcSsensitivityreverse!`](@ref), allocated once per batch of signal
frequencies rather than at every frequency. Each thread of
[`hblinsolve`](@ref) owns its own set; the [`ReverseSensitivity`](@ref)
itself is shared read only.
"""
struct ReverseSensitivityBuffers
    P::Vector{Complex{Float64}}
    Q::Vector{Complex{Float64}}
    # the output functional covectors of a chunk of output pairs, and their
    # solutions through the transposed pump Jacobian, batched so the sparse
    # solver amortizes its per-call overhead over many right hand sides.
    G::Matrix{Complex{Float64}}
    Psi::Matrix{Complex{Float64}}
    eta::Matrix{Complex{Float64}}
    c::Matrix{Complex{Float64}}
    # the zero padded input and the single output grid of the transposed
    # transform: the holomorphic and antiholomorphic halves are folded into
    # eta one after the other, so their transforms need not coexist.
    padded::Array{Complex{Float64}}
    tgrid::Array{Complex{Float64}}
end

# the byte budget of the two nr x chunk complex right hand side and
# solution buffers of one frequency batch of the reverse contraction
# (together they cost 32*nr*chunk bytes, so at this budget a pump system of
# nr = 8000 real unknowns gets a chunk of about 128 columns, and a system a
# hundred times larger degrades gracefully toward single column solves).
# batching amortizes the per-call overhead of the sparse triangular solves,
# measured about 3x for tens of right hand sides, and the budget rather
# than a fixed column count keeps the per batch memory bounded for large
# pump systems; nbatches controls how many batches (and therefore budgets)
# are active at once.
const REVERSESENSITIVITYCHUNKBYTES = 32*2^20

function ReverseSensitivityBuffers(rev::ReverseSensitivity, NPM::Integer)
    sys = rev.op.sys
    NLj = size(sys.phimatrix)[end]
    NF = length(sys.phimatrix)
    nr = size(rev.op.jacobian, 1)
    chunk = clamp(REVERSESENSITIVITYCHUNKBYTES ÷ (32*nr), 1, NPM^2)
    return ReverseSensitivityBuffers(
        zeros(Complex{Float64}, NF), zeros(Complex{Float64}, NF),
        zeros(Complex{Float64}, nr, chunk),
        zeros(Complex{Float64}, nr, chunk),
        zeros(Complex{Float64}, size(rev.T, 1), NLj),
        zeros(Complex{Float64}, size(rev.T, 2), NLj),
        zeros(Complex{Float64}, size(sys.phitd)),
        zeros(Complex{Float64}, size(sys.phitd)))
end

"""
    calcbranchtimedomainmap(sys, Nmodes, NLj)

The matrix of the map from the branch fluxes of one Josephson junction to its
physical time domain branch flux, with the real and the imaginary part of
each mode as separate columns. The map is the same for every junction,
because the packing and the inverse transform act on each junction
independently, so it is built once by transforming unit branch fluxes with
[`phivectortomatrix!`](@ref) and `applyifft!`, which keeps the conjugate mode
bookkeeping inside the functions which define it. This is the transpose of
the linear map from the unknowns to the time domain branch fluxes, restricted
to one junction.
"""
function calcbranchtimedomainmap(sys, Nmodes::Integer, NLj::Integer)
    branch = sys.Ljb.nzind[1]
    Nbranches = size(sys.Rbnm, 1) ÷ Nmodes
    T = zeros(Float64, length(sys.phitd) ÷ NLj, 2*Nmodes)
    fd = similar(sys.phimatrix)
    td = similar(sys.phitd)
    tdflat = reshape(td, :, NLj)
    bv = zeros(Complex{Float64}, Nbranches*Nmodes)
    for m in 1:Nmodes
        for (q, part) in enumerate((one(Complex{Float64}), im))
            fill!(bv, 0)
            bv[(branch-1)*Nmodes+m] = part
            fill!(fd, 0)
            phivectortomatrix!(bv[sys.Ljbm.nzind], fd, sys.freqindexmap,
                sys.conjsourceindices, sys.conjtargetindices, NLj)
            applyifft!(td, fd, sys.irfftplan)
            T[:, 2*(m-1)+q] .= view(tdflat, :, 1)
        end
    end
    return T
end

"""
    ReverseSensitivity(op::HBOperatingPoint, lsys, dFr)

Precompute the reverse mode contraction data: the branch flux map
([`calcbranchtimedomainmap`](@ref)), the transform of the pump harmonic grid,
the row and the column of each nonzero of the linearized system matrix, and
the offset of each entry of the augmented state in its real representation.
"""
function ReverseSensitivity(op::HBOperatingPoint, lsys, dFr)
    # the operating point, and therefore the branch flux map and the
    # incidence lists, live on the pump mode grid, not the signal mode grid
    Nmodes = op.Nmodes
    sys = op.sys
    # the per-frequency contraction reads the cached time domain sine of the
    # branch fluxes, so pin it to the operating point here, in this serial
    # constructor: the threads of hblinsolve only read it, and updating a
    # shared cache from them would be a race. without this the cache would
    # hold whatever point the shared evaluation object was last used at.
    setpoint!(sys, op.x)
    _ensuresin!(sys)
    NLj = size(sys.phimatrix)[end]
    A = lsys.Asparse
    nzrow = zeros(Int, nnz(A))
    nzcol = zeros(Int, nnz(A))
    rows = rowvals(A)
    for j in axes(A, 2)
        for p in nzrange(A, j)
            nzrow[p] = rows[p]
            nzcol[p] = j
        end
    end
    # the offset of complex entry j in the real representation. note that
    # isreal is indexed by mode, not by entry of the augmented state.
    isrealmode = op.modelayout.isreal
    nmd = length(isrealmode)
    realindexmap = zeros(Int, length(op.x))
    k = 1
    for j in eachindex(realindexmap)
        realindexmap[j] = k
        k += isrealmode[(j-1) % nmd + 1] ? 1 : 2
    end
    Nbranches = size(sys.Rbnm, 1) ÷ Nmodes
    return ReverseSensitivity(op, SparseMatrixCSC{Float64,Int}(dFr),
        calcbranchtimedomainmap(sys, Nmodes, NLj),
        plan_applyffttranspose(sys.phimatrix, sys.phitd), nzrow, nzcol,
        realindexmap, branchnodesandsigns(sys.Rbnm, Nmodes, Nbranches))
end

"""
    calcSsensitivityreverse!(Ssensitivity, rev::ReverseSensitivity, lsys,
        phin, phinadjoint, gamma, beta, cache,
        bufs::ReverseSensitivityBuffers)

Add the contribution of the shift of the pump operating point to the
scattering parameter sensitivities, contracting in the reverse order so that
the cost per component is a sparse inner product rather than a product
against a matrix which is dense on the sparsity structure of the linearized
system.

For each pair of output and input port modes `(a,b)`, the transpose of the
Josephson scatter of [`addjosephsonterm!`](@ref) gives

    transpose(lam_a)*dAop_k*phi_b
        = sum_s P[s]*dcos_k[s] + Q[s]*conj(dcos_k[s]),

with `P` and `Q` accumulated over the scatter lists of the plan. Since
`dcos_k` is the directional derivative of the Fourier coefficients of
`cos(phi_b(t))` along the operating point shift
([`cosdirectionalderivative!`](@ref)), that is a linear functional of the
shift, and with

    alpha = Em*P,  gam = conj(Em)*Q,  eta = -sin(phi_b(t)).*(alpha + gam)

its covector is the transpose of the branch flux map applied to `eta`.
Finally the implicit function theorem gives `dx_k = -inv(J)*dF_k`, so pushing
that covector through the transposed Jacobian once per output pair leaves a
sparse inner product with `dF_k` for each component. The cost per signal
frequency is `(Nports*Nmodes)^2` transposed solves, independent of the number
of components, instead of one product against the full sparsity structure per
component. The solves are batched into multi right hand side calls whose
column count is set by the `REVERSESENSITIVITYCHUNKBYTES` budget, which
amortizes the per-call overhead of the sparse triangular solves while
keeping the per batch work matrices memory bounded.

The transform of the pump harmonic grid is applied one dimension at a time,
so this supports any number of pump tones.
"""
function calcSsensitivityreverse!(Ssensitivity, rev::ReverseSensitivity,
    lsys, phin, phinadjoint, gamma, beta, cache,
    bufs::ReverseSensitivityBuffers)

    op = rev.op
    # the pump mode count: the operating point shift lives on the pump grid
    Nmodes = op.Nmodes
    sys = op.sys
    plan = lsys.complexjacobianplan
    NPM = size(phin, 2)
    NLj = size(sys.phimatrix)[end]
    Ncomponents = size(rev.dFr, 2)
    isrealmode = op.modelayout.isreal
    nmd = length(isrealmode)

    # the cached time domain sine of the pump branch fluxes, pinned to the
    # operating point by the ReverseSensitivity constructor and read only
    # here.
    sintd = reshape(sys.sintd, :, NLj)
    P = bufs.P
    Q = bufs.Q
    G = bufs.G
    Psi = bufs.Psi
    eta = bufs.eta
    c = bufs.c
    padded = bufs.padded
    tgrid = bufs.tgrid
    wsize = size(sys.phimatrix)
    dFrows = rowvals(rev.dFr)
    dFvals = nonzeros(rev.dFr)

    # The output pairs are independent, so their solves through the
    # transposed pump Jacobian are batched: the covectors of a chunk of
    # pairs are accumulated as the columns of G and pushed through the
    # factorization in one multi right hand side call, which amortizes the
    # per-call overhead of the sparse triangular solves over the chunk.
    wcov = zeros(eltype(P), plan.josephson.n)
    pairs = vec(CartesianIndices((NPM, NPM)))
    for chunk in Iterators.partition(eachindex(pairs), size(G, 2))
        for (col, pi) in enumerate(chunk)
            a, b = Tuple(pairs[pi])
            # the transpose of the Josephson scatter
            fill!(P, 0)
            fill!(Q, 0)
            # the covector per destination, then the transpose of the
            # Josephson map applied to it. The plain and conjugated halves
            # separate by the sign of the mode coupling index.
            @inbounds for p in eachindex(wcov)
                wcov[p] = phinadjoint[rev.nzrow[p],a]*phin[rev.nzcol[p],b]
            end
            josephsonadjoint!(P, Q, plan.josephson, wcov)

            # the transpose of the cos directional derivative: a forward
            # transform of the zero padded coefficients through the plan of
            # the operating point grid (the transform matrix is symmetric),
            # and the antiholomorphic half by conjugating around the same
            # transform. the two halves are folded into eta one after the
            # other through the same output grid.
            applyffttranspose!(tgrid, reshape(P, wsize), padded, rev.fftplan)
            tf = reshape(tgrid, :, NLj)
            @inbounds for i in eachindex(eta)
                eta[i] = -sintd[i]*tf[i]
            end
            @inbounds for i in eachindex(Q)
                Q[i] = conj(Q[i])
            end
            applyffttranspose!(tgrid, reshape(Q, wsize), padded, rev.fftplan)
            @inbounds for i in eachindex(eta)
                eta[i] -= sintd[i]*conj(tf[i])
            end
            mul!(c, transpose(rev.T), eta)

            # the transpose of the branch flux map, into the real
            # representation of the augmented state, as one column of G
            g = view(G, :, col)
            fill!(g, 0)
            @inbounds for (jj, branch) in enumerate(sys.Ljb.nzind)
                for m in 1:Nmodes
                    cre = c[2*(m-1)+1, jj]
                    cim = c[2*(m-1)+2, jj]
                    holo = (cre - im*cim)/2
                    anti = (cre + im*cim)/2
                    for (node, sgn) in rev.branchnodes[branch]
                        j = (node-1)*Nmodes + m
                        k = rev.realindexmap[j]
                        g[k] += sgn*(holo + anti)
                        if !isrealmode[(j-1) % nmd + 1]
                            g[k+1] += sgn*im*(holo - anti)
                        end
                    end
                end
            end
        end

        # through the transposed Jacobian once for the whole chunk
        ncols = length(chunk)
        Gc = view(G, :, 1:ncols)
        Psic = view(Psi, :, 1:ncols)
        trysolvetranspose!(Psic, cache.factorization, Gc)

        # the sparse inner products with the residual derivatives
        for (col, pi) in enumerate(chunk)
            a, b = Tuple(pairs[pi])
            @inbounds for k in 1:Ncomponents
                acc = zero(Complex{Float64})
                for r in nzrange(rev.dFr, k)
                    acc += Psi[dFrows[r], col]*dFvals[r]
                end
                Ssensitivity[a,b,k] += gamma[a]*beta[b]*acc
            end
        end
    end
    return nothing
end

"""
    calcoperatingpointstamps(op::HBOperatingPoint, lsys, dx)

Calculate the contribution of a shift of the pump operating point to the
derivative of the linearized system matrix, for each column of `dx`. The
linearized system matrix depends on the operating point only through the
Fourier coefficients of `cos(phi_b(t))` of the Josephson junction branch
fluxes, so the contribution is the directional derivative of those
coefficients along the operating point shift
([`cosdirectionalderivative!`](@ref)) scattered into the system matrix
through the same plan which assembles it ([`addjosephsonterm!`](@ref)), so
the mode coupling and its truncation agree exactly. Returns a vector of
nonzero value vectors aligned with the sparsity structure of the system
matrix, which are frequency independent.
"""
function calcoperatingpointstamps(op::HBOperatingPoint, lsys, dx)
    setpoint!(op.sys, op.x)
    dcos = similar(op.sys.phimatrix)
    stamps = Vector{Vector{Complex{Float64}}}(undef, size(dx,2))
    for k in axes(dx, 2)
        cosdirectionalderivative!(dcos, op.sys,
            Vector{Complex{Float64}}(view(dx,:,k)))
        nzval = zeros(Complex{Float64}, nnz(lsys.Asparse))
        addjosephsonterm!(nzval, lsys.complexjacobianplan, dcos)
        stamps[k] = nzval
    end
    return stamps
end

"""
    SensitivityStamp(kind, M, indexmap, freqsubstindices, nzval, portindex)

The derivative of the linearized harmonic balance system matrix with respect
to a relative (logarithmic) perturbation of one component value, `p -> r*p`
evaluated at `r = 1`. The system matrix is affine in `C`, `1/R`, `1/L` and
`1/Lj`, so the derivative is that component's own contribution to the system
matrix, with a sign: positive for the capacitance and negative for the
quantities which enter inversely.

`kind` selects how the stamp is assembled at each signal frequency, mirroring
[`assemblesystemmatrix!`](@ref):

- `:C`: `-M*wmodes2m`, the component's capacitance matrix,
- `:G`: `-im*M*wmodesm`, the component's conductance matrix, either the nodal
    stamp or, for a resistor promoted by the modified nodal analysis
    formulation, the conductances of its constitutive equations,
- `:invL`: `-M`, the component's inverse inductance matrix,
- `:Lj`: the constant nonzero values `nzval`, the negative of the pump
    modulated Josephson contribution of that junction alone.

`indexmap` maps the nonzeros of `M` into the nonzeros of the system matrix
and `freqsubstindices` locates any symbolic entries. `portindex` is the
index of the port whose impedance this component is, or zero, which selects
the additional wave normalization term of [`calcSsensitivity!`](@ref).
"""
struct SensitivityStamp
    kind::Symbol
    rows::Vector{Int}
    cols::Vector{Int}
    vals::Vector{Complex{Float64}}
    portindex::Int
end
"""
    calcsensitivitystamps(sensitivityindices, psc, cg, nm, lsys, phimatrix,
        mnaindices, coupledbranches, Nnodalmna, Nmodes, Nnodes)

Build the [`SensitivityStamp`](@ref) of each component in
`sensitivityindices`. The classification and the raw one-component matrices
come from [`componentstamp`](@ref), which is shared with the residual
derivatives of [`calcresidualsensitivity`](@ref), so the two grids cannot
disagree on which components are supported or how they are built. The
Josephson junction stamp is the pump modulated contribution of that
junction alone, obtained by scattering the Fourier coefficients of
`cos(phi(t))` of that junction through the same plan
([`addjosephsonterm!`](@ref)) which assembles the system matrix, so the mode
coupling and its truncation agree exactly.
"""
function calcsensitivitystamps(sensitivityindices, psc::ParsedSortedCircuit,
    cg::CircuitGraph, nm::CircuitMatrices, lsys,
    phimatrix, mnaindices, coupledbranches, Nnodalmna, Nmodes, Nnodes)

    Ntot = size(lsys.Asparse, 1)
    stamps = Vector{SensitivityStamp}(undef, length(sensitivityindices))
    lookups = componentlookups(mnaindices, coupledbranches, nm.Ljb)
    portordinal = Dict(idx => p
        for (p, idx) in enumerate(nm.portimpedanceindices))

    for (k, idx) in enumerate(sensitivityindices)
        portindex = get(portordinal, idx, 0)
        kind, info = componentstamp(idx, psc, cg, nm, lookups,
            Nmodes, Nnodes)
        if kind == :C
            stamps[k] = tripletstamp(:C, mnapadto(info, Ntot), portindex)
        elseif kind == :G
            stamps[k] = tripletstamp(:G, mnapadto(info, Ntot), portindex)
        elseif kind == :Gpromoted
            # the promoted resistors keep their constitutive equations as
            # explicit rows, so the conductance appears there instead.
            M = maskrowscolumns(lsys.AmnaG,
                Nnodalmna+(info-1)*Nmodes+1, Nnodalmna+info*Nmodes,
                1, size(lsys.AmnaG, 2))
            stamps[k] = tripletstamp(:G, M, portindex)
        elseif kind == :invL
            stamps[k] = tripletstamp(:invL, mnapadto(info, Ntot), portindex)
        else # :Lj
            # the contribution of this junction alone, by zeroing the Fourier
            # coefficients of every other junction before the scatter.
            onejunction = zero(phimatrix)
            selectdim(onejunction, ndims(onejunction), info) .=
                selectdim(phimatrix, ndims(phimatrix), info)
            nzval = zeros(Complex{Float64}, nnz(lsys.Asparse))
            addjosephsonterm!(nzval, lsys.complexjacobianplan, onejunction)
            rmul!(nzval, -1)
            stamps[k] = tripletstamp(:Lj,
                SparseMatrixCSC(size(lsys.Asparse,1), size(lsys.Asparse,2),
                    copy(lsys.Asparse.colptr), copy(lsys.Asparse.rowval),
                    nzval), portindex)
        end
    end
    return stamps
end

# pad a matrix with empty rows and columns for the auxiliary variables of the
# modified nodal analysis formulation, and convert to the element type of the
# system matrix.
function mnapadto(M::SparseMatrixCSC, Ntot::Integer)
    padded = size(M,1) == Ntot ? M : mnapad(M, Ntot - size(M,1))
    return SparseMatrixCSC{Complex{Float64},Int}(padded)
end

# convert a component matrix to the compact triplet form of a
# SensitivityStamp, dropping structural zeros. The stamps of the individual
# components are extremely sparse compared with the system matrix (a
# capacitor to ground touches one node), so the contraction is driven by
# these entries rather than by the sparsity structure of the system matrix.
function tripletstamp(kind::Symbol, M::SparseMatrixCSC, portindex::Integer)
    I, J, V = findnz(M)
    keep = .!iszero.(V)
    return SensitivityStamp(kind, I[keep], J[keep],
        Complex{Float64}.(V[keep]), portindex)
end

"""
    sensitivitystampvalue(stamp::SensitivityStamp, t::Integer, wmodes,
        Nmodes)

The value of entry `t` of the derivative of the linearized harmonic balance
system matrix with respect to a relative perturbation of the component of
`stamp`, at the mode frequencies `wmodes`. Applies the same per mode
frequency scaling and negative frequency mode conjugation as
[`assemblesystemmatrix!`](@ref), which are indexed by the column.
"""
@inline function sensitivitystampvalue(stamp::SensitivityStamp, t::Integer,
    wmodes, Nmodes)

    stamp.kind == :Lj && return stamp.vals[t]
    m = (stamp.cols[t] - 1) % Nmodes + 1
    w = wmodes[m]
    v = modevalue(stamp.vals[t], w)
    if stamp.kind == :C
        return -v*w^2
    elseif stamp.kind == :G
        return -im*v*w
    else
        return -v
    end
end

"""
    calcsensitivityscaling!(gamma, beta, inputwave, bnm,
        portimpedanceindices, portimpedances, componenttypes, nodeindices,
        wmodes, Nmodes)

Calculate the scalars which convert the derivative of the node fluxes into
the derivative of the scattering parameters. Writing the output wave of
[`calcinputoutput!`](@ref) as a linear functional of the node fluxes, the
part which depends on them is
`(1/2)*kval*(1 + conj(Z)/Z)*im*w_n*(phi_n1 - phi_n2)`, and the scattering
parameters divide by the input wave, so

    dS[(j,n),(i,m)] = gamma[(j,n)]*beta[(i,m)]
                      *(dphi[node1,(j,n)] - dphi[node2,(j,n)])

with `gamma = (1/2)*kval*(1 + conj(Z)/Z)*im*w_n/s_{(j,n)}` and
`beta = 1/inputwave[(i,m),(i,m)]`. The node flux difference is contracted
with the adjoint solution by [`calcSsensitivity!`](@ref), which is why the
`im*w_n` factor of the port voltage is folded into `gamma` here: it is
exactly the factor relating the adjoint source vector to the source vector
of the forward problem.

`s_{(j,n)}` is the source current of the port's own unit drive
([`calcsourcecurrent`](@ref) on the diagonal), `±1` depending on whether the
canonical orientation of the port branch in the incidence matrix agrees with
the node order of the port component. The adjoint solution is the solve
against the source columns of `bnm`, which carry the canonical branch
orientation, while the output functional differences the node fluxes in the
component node order, so their ratio enters the contraction. Without it the
sensitivities of any output at a port written with its nodes in the opposite
order of the branch orientation would have the wrong sign, even though the
scattering parameters themselves, which use `s` consistently in both the
input and the output waves, would be correct.
"""
function calcsensitivityscaling!(gamma, beta, inputwave, bnm,
    portimpedanceindices, portimpedances, componenttypes, nodeindices,
    wmodes, Nmodes)

    for i in eachindex(portimpedanceindices)
        for j in 1:Nmodes
            row = (i-1)*Nmodes + j
            portimpedance = calcimpedance(portimpedances[i],
                componenttypes[portimpedanceindices[i]], wmodes[j], nothing)
            kval = portwavescale(portimpedance, wmodes[j])
            # the orientation of the port branch relative to the node order
            # of the port component, from the source current of the port's
            # own unit drive.
            sourcecurrent = calcsourcecurrent(
                nodeindices[1,portimpedanceindices[i]],
                nodeindices[2,portimpedanceindices[i]], bnm, Nmodes, j, row)
            gamma[row] = iszero(sourcecurrent) ? 0 :
                1/2*kval*(1 + conj(portimpedance)/portimpedance)*
                im*wmodes[j]/sourcecurrent
            beta[row] = iszero(inputwave[row,row]) ? 0 :
                1/inputwave[row,row]
        end
    end
    return nothing
end

"""
    calcSsensitivity!(Ssensitivity, stamps, dA, dAphin, phin, phinadjoint,
        S, gamma, beta, wmodes, Nmodes, symfreqvar)

Calculate the derivative of the scattering matrix with respect to a relative
(logarithmic) perturbation of each component value, `p -> r*p` evaluated at
`r = 1`, at the pump operating point, with the adjoint method. Overwrites
`Ssensitivity`.

Differentiating the linearized system `A*phi = b`, whose source terms do not
depend on any component value, gives `dphi = -inv(A)*dA*phi`, so with the
adjoint solutions `lam` of the transposed system driven by the output
functionals of [`calcsensitivityscaling!`](@ref),

    dS[(j,n),(i,m)] = -gamma[(j,n)]*beta[(i,m)]
                      *transpose(lam[:,(j,n)])*dA*phi[:,(i,m)].

The adjoint source vectors are the source vectors of the forward problem
scaled by `im*w_n`, which is folded into `gamma`, so `phinadjoint`, the
solution of the transposed system already computed for the noise and quantum
efficiency calculations, is used directly.

When the perturbed component is itself a port impedance the wave
normalization of [`calcinputoutput!`](@ref) moves as well, contributing the
additional closed form term `-(1/2)*(P*(S+I) + (S+I)*P)` with `P` the
projector onto that port, which is exact for constant real port impedances.
"""
function calcSsensitivity!(Ssensitivity, stamps, dAop, dA, dAphin, phin,
    phinadjoint, S, gamma, beta, contraction, wmodes, Nmodes, symfreqvar)

    NPM = size(phin, 2)
    for (k, stamp) in enumerate(stamps)
        # the component's own contribution, driven by the entries of its
        # stamp rather than by the sparsity structure of the system matrix,
        # of which the stamp of a single component is a tiny part.
        fill!(contraction, 0)
        @inbounds for t in eachindex(stamp.rows)
            i = stamp.rows[t]
            j = stamp.cols[t]
            v = sensitivitystampvalue(stamp, t, wmodes, Nmodes)
            iszero(v) && continue
            for b in 1:NPM
                vphi = v*phin[j,b]
                iszero(vphi) && continue
                for a in 1:NPM
                    contraction[a,b] += phinadjoint[i,a]*vphi
                end
            end
        end

        # the contribution of the shift of the pump operating point, which is
        # frequency independent but dense on the sparsity structure of the
        # system matrix, so it goes through a sparse matrix vector product.
        if !isempty(dAop)
            copyto!(nonzeros(dA), dAop[k])
            mul!(dAphin, dA, phin)
            mul!(contraction, transpose(phinadjoint), dAphin, 1, 1)
        end

        for b in 1:NPM
            for a in 1:NPM
                Ssensitivity[a,b,k] = -gamma[a]*beta[b]*contraction[a,b]
            end
        end

        # the wave normalization term when the component is a port impedance
        if stamp.portindex > 0
            p = stamp.portindex
            for b in 1:NPM
                for a in 1:NPM
                    sI = S[a,b] + (a == b ? 1 : 0)
                    correction = zero(Complex{Float64})
                    if (a-1) ÷ Nmodes + 1 == p
                        correction += sI/2
                    end
                    if (b-1) ÷ Nmodes + 1 == p
                        correction += sI/2
                    end
                    Ssensitivity[a,b,k] -= correction
                end
            end
        end
    end
    return nothing
end


"""
    StagedStageInfo

One attempted stage of [`stagedhbnlsolve`](@ref), stored in
`solverinfo.stages` of the returned solution -- every attempt appears, in
the order it ran, including stalled steps and growth retreats, so the
whole continuation walk can be examined afterwards.

# Fields
- `label`: `"staged"`.
- `converged`: whether this attempt's inner solve converged.
- `iterations`: total inner Newton iterations of the attempt.
- `grid`: the harmonic truncation the attempt solved on.
- `sfrom`: the last accepted drive fraction before the attempt.
- `starget`: the drive fraction the attempt targeted.
- `ds`: `starget - sfrom` (negative for a growth retreat).
- `action`: `:advance` (drive step on the current grid), `:grow` (first
    solve on a larger grid after carrying a converged point up), or
    `:final` (the full-drive solve on the finest grid).
- `accepted`: whether the attempt's result was kept as the new operating
    point (a stalled attempt is recorded but not accepted).
- `seconds`: wall time of the attempt, including the stage's system
    assembly.
- `finalresidual`: the residual norm the attempt ended at.
- `inner`: the inner solver's own stage records ([`IterationInfo`](@ref)),
    with their Krylov linear-solve diagnostics when the inner method is
    `:newtonkrylov`.
"""
struct StagedStageInfo <: AbstractStageInfo
    label::String
    converged::Bool
    iterations::Int
    grid::Tuple
    sfrom::Float64
    starget::Float64
    ds::Float64
    action::Symbol
    accepted::Bool
    seconds::Float64
    finalresidual::Float64
    inner::Vector
end

function Base.show(io::IO, ::MIME"text/plain", r::StagedStageInfo)
    print(io, "StagedStageInfo: grid=", r.grid, " s ", round(r.sfrom, digits = 4),
        " -> ", round(r.starget, digits = 4), " ", r.action,
        r.accepted ? " accepted" : " stalled",
        " newton=", r.iterations,
        " |F|=", round(r.finalresidual, sigdigits = 2),
        " (", round(r.seconds, digits = 2), " s)")
end

"""
    defaultgridladder(Nharmonics::NTuple{N,Int})

The default coarse-to-fine harmonic grid ladder for [`stagedhbnlsolve`](@ref):
repeated halving of every dimension down to two harmonics, finest last.
"""
function defaultgridladder(Nharmonics::NTuple{N,Int}) where {N}
    grids = [Nharmonics]
    g = Nharmonics
    while any(x -> x > 2, g)
        g = map(x -> max(2, cld(x, 2)), g)
        g == grids[end] && break
        push!(grids, g)
    end
    return reverse(grids)
end

# embed a converged solution into a larger grid's initial guess by matching
# mode tuples; modes the small grid does not carry start at zero
function stagedembed(out, bigmodes)
    small = reshape(Array(out.nodeflux), out.Nmodes, :)
    pos = Dict(m => i for (i, m) in enumerate(bigmodes))
    X = zeros(ComplexF64, length(bigmodes), size(small, 2))
    for (i, m) in enumerate(out.modes)
        haskey(pos, m) && (X[pos[m], :] = small[i, :])
    end
    return vec(X)
end

"""
    stagedhbnlsolve(w, Nharmonics, sources, circuit, circuitdefs; kwargs...)

Source continuation on an adaptively grown harmonic grid, reached through
`hbnlsolve(...; method = :staged)`.

Near a critical drive the Newton basin is small and the iteration count
large, so those iterations are spent where they are cheap: the drive is
climbed in warm started steps on a small harmonic truncation, and each
larger grid is warm started from the last by matching mode tuples. The
schedule is adaptive in both directions because every truncation has its own
solvability boundary, and the boundaries are not monotone in the grid
(measured: the (2,2) truncation of a 64 junction RPM at strong two tone
drive converges at 98.4% of the drive while (4,2) does not): a stalled drive
step is first halved, a stall at the minimum step grows the grid at the
current converged drive, and a carried point that fails to reconverge after
growth retreats the drive on the new grid until it converges. Interior
points converge only to `interiorftol` under a small iteration budget --
they exist to keep the iterate inside the basin -- and the single expensive
solve, the finest grid at full drive, starts inside the basin and gets the
caller's `ftol` and `iterations`.

A point carried to the finest grid at full drive that stalls there without
ever having converged on that grid is not diagnosed as a fold: the drive is
retreated on the finest grid and climbed back, exactly as a mid-ladder
growth retreats, since the stalled drive value was established on different
physics. Only a stall from a point converged on the finest grid itself
brackets a fold: between the
last converged drive fraction and the stalled one the solution branch ends
(the self oscillation threshold), and the error says so, with the bracket.
Distinguishing that from a merely hard problem is something no cold started
method can do. Measured on a 128 junction RPM at equal 1.6 uA two tone
drive -- reported "unreachable" by every direct method -- the walk brackets
the fold at 91.9-93.4% of the drive: the operating point does not exist.

# Keywords (through `stagedkwargs`)
- `grids = defaultgridladder(Nharmonics)`: the coarse-to-fine ladder;
    finest entry must equal `Nharmonics`.
- `s0 = 0.5`: the first drive fraction attempted.
- `smin = 0.02`: the minimum drive step; a stall below it grows the grid.
- `interiorftol = 1e-7`, `interioriterations = 60`: tolerance and Newton
    budget of the interior points. The budget is deliberately small: a
    stalled probe is evident within tens of iterations, and interior stalls
    are the pure overhead of the walk (measured accepted interior points
    take 4-24 iterations, 51 near a pole).
- `innermethod = :newtonkrylov`: the solver for every stage.
- `interiorescalation = false`: whether interior stage solves may escalate
    their preconditioner to the full Jacobian. Off by default: an interior
    probe exists only to produce a cheap warm start, escalating inside one
    pays a full factorization for a throwaway point, and at high tone
    counts that factorization may not fit in memory at all (measured: a
    527 mode three tone grid ran out of 30 GB the moment an interior probe
    escalated, while the whole staged walk to that point took 19 seconds).
    A probe that fails without escalation is simply a stall, which the
    schedule answers with a smaller step. The final solve keeps the
    caller's escalation behavior.
- `maxattempts = 60`: bound on the total number of stage solves.
- `verbose = false`: print one line per stage solve.
"""
function stagedhbnlsolve(w::NTuple{N,Number}, Nharmonics::NTuple{N,Int},
    sources, circuit, circuitdefs;
    grids::Vector{NTuple{N,Int}} = defaultgridladder(Nharmonics),
    s0 = 0.5, smin = 0.02, interiorftol = 1e-7,
    interioriterations::Integer = 60, innermethod::Symbol = :newtonkrylov,
    interiorescalation::Bool = false,
    maxattempts::Integer = 60, verbose::Bool = false,
    iterations = 1000, maxharmonics::NTuple{N,Number} = Nharmonics,
    maxintermodorder = Inf, dc::Bool = false, odd::Bool = true,
    even::Bool = false, ftol = 1e-8, symfreqvar = nothing,
    sorting = :number, keyedarrays::Bool = true,
    sensitivitynames::Vector{String} = String[],
    returnoperatingpoint::Bool = false, krylovcouplingmodes = :none,
    krylovrecycle::Integer = 0, krylovharvest::Integer = 8,
    krylovkwargs::NamedTuple = (;), factorization = nothing,
    backend = CPU(), precision::Type{<:AbstractFloat} = Float64) where {N}

    isempty(grids) && throw(ArgumentError("`grids` must not be empty."))
    grids[end] == Nharmonics || throw(ArgumentError(
        lazy"the finest grid $(grids[end]) must equal `Nharmonics` = $(Nharmonics)."))
    0 < s0 <= 1 || throw(ArgumentError(lazy"`s0` = $(s0) must be in (0, 1]."))
    innermethod === :staged && throw(ArgumentError(
        "`innermethod` must be a non-staged method."))

    modesof(grid) = removeconjfreqs(truncfreqs(calcfreqsrdft(grid);
        dc = dc, odd = odd, even = even,
        maxintermodorder = maxintermodorder,
        maxharmonics = map(min, maxharmonics, grid))).modes
    scaled(s) = [(mode = t.mode, port = t.port, current = s*t.current)
        for t in sources]
    solve(grid, s, x0, final) = hbnlsolve(w, grid, scaled(s), circuit,
        circuitdefs; dc = dc, odd = odd, even = even,
        maxintermodorder = maxintermodorder,
        maxharmonics = map(min, maxharmonics, grid),
        method = innermethod, x0 = x0, symfreqvar = symfreqvar,
        sorting = sorting,
        keyedarrays = final ? keyedarrays : false,
        sensitivitynames = final ? sensitivitynames : String[],
        returnoperatingpoint = final ? returnoperatingpoint : false,
        krylovcouplingmodes = krylovcouplingmodes,
        krylovrecycle = krylovrecycle, krylovharvest = krylovharvest,
        krylovkwargs = (final || interiorescalation) ? krylovkwargs :
            merge((; krylovescalate = typemax(Int)), krylovkwargs),
        factorization = factorization,
        backend = backend, precision = precision,
        ftol = final ? ftol : interiorftol,
        iterations = final ? iterations : interioriterations)

    gi = 1
    s = 0.0             # last converged drive fraction on the current grid
    ds = s0
    x = nothing         # its solution, in the raw nodeflux layout
    out = nothing
    attempts = 0
    stagerecords = AbstractStageInfo[]
    pendinggrow = false
    record = function (cand, grid, sfrom, starget, action, accepted, secs)
        si = cand.solverinfo
        push!(stagerecords, StagedStageInfo("staged", si.converged,
            sum(st -> st.iterations, si.stages; init = 0), grid,
            Float64(sfrom), Float64(starget), Float64(starget - sfrom),
            action, accepted, secs, Float64(si.finalresidual), si.stages))
        return nothing
    end
    while true
        attempts += 1
        attempts > maxattempts && error(
            lazy"the staged schedule did not converge in `maxattempts` = $(maxattempts) stage solves; last converged drive fraction $(s) on grid $(grids[gi]).")
        starget = min(1.0, s + ds)
        final = gi == length(grids) && starget >= 1.0
        t0 = time_ns()
        cand = solve(grids[gi], starget, x, final)
        ok = cand.solverinfo.converged
        # the first solve after carrying a full-drive point to a larger
        # grid is a growth, not a drive advance
        record(cand, grids[gi], s, starget,
            final ? :final : (pendinggrow ? :grow : :advance), ok,
            (time_ns() - t0)/1e9)
        ok && (pendinggrow = false)
        if verbose
            st = cand.solverinfo.stages[end]
            println("staged: grid=", grids[gi], " s=", round(starget, digits = 4),
                " newton=", st.iterations,
                " |F|=", round(cand.solverinfo.finalresidual, sigdigits = 2),
                ok ? "" : " STALL")
        end
        if ok
            out = cand
            s = starget
            final && break
            if s >= 1.0
                # full drive reached on a coarse grid: grow
                x = stagedembed(out, modesof(grids[gi+1]))
                gi += 1
                pendinggrow = true
            else
                x = vec(reshape(Array(out.nodeflux), out.Nmodes, :))
                ds = min(2*ds, 1.0 - s)
            end
        elseif starget - s > smin
            # the effective step, not `ds`: a target capped at full drive
            # must not be re-attempted identically while `ds` halves above
            # the cap
            ds = (starget - s)/2
        elseif gi < length(grids)
            # the current grid's own pole: grow at the converged drive,
            # retreating the drive on the new grid if the carried point does
            # not reconverge there (the boundaries are not monotone)
            isnothing(x) && error(
                "the first stage stalled at its first drive step; lower `s0`.")
            bigmodes = modesof(grids[gi+1])
            bigx = stagedembed(out, bigmodes)
            gi += 1
            reconverged = false
            for f in (1.0, 0.9, 0.8, 0.65, 0.5)
                t0 = time_ns()
                re = solve(grids[gi], f*s, bigx, false)
                record(re, grids[gi], s, f*s, :grow,
                    re.solverinfo.converged, (time_ns() - t0)/1e9)
                verbose && println("staged: grow -> ", grids[gi], " at s=",
                    round(f*s, digits = 4), " |F|=",
                    round(re.solverinfo.finalresidual, sigdigits = 2),
                    re.solverinfo.converged ? "" : " STALL")
                if re.solverinfo.converged
                    out = re
                    s = f*s
                    x = vec(reshape(Array(out.nodeflux), out.Nmodes, :))
                    reconverged = true
                    break
                end
            end
            reconverged || error(
                lazy"the carried point did not reconverge on grid $(grids[gi]) even at half its drive; the truncation boundaries of the ladder are too far apart. Add an intermediate grid.")
            ds = s0/2
        elseif pendinggrow
            # The stalled point was CARRIED to the finest grid from a
            # coarser one and has never converged here, so a fold diagnosis
            # would rest on a drive value established on different physics.
            # Retreat the drive on the finest grid and walk back up, exactly
            # as a mid-ladder growth retreats on its new grid; only a stall
            # from a point converged on the finest grid itself brackets a
            # fold.
            isnothing(x) && error(
                "the first stage stalled at its first drive step; lower `s0`.")
            reconverged = false
            for f in (0.9, 0.8, 0.65, 0.5)
                t0 = time_ns()
                re = solve(grids[gi], f*s, x, false)
                record(re, grids[gi], s, f*s, :grow,
                    re.solverinfo.converged, (time_ns() - t0)/1e9)
                verbose && println("staged: retreat on ", grids[gi], " at s=",
                    round(f*s, digits = 4), " |F|=",
                    round(re.solverinfo.finalresidual, sigdigits = 2),
                    re.solverinfo.converged ? "" : " STALL")
                if re.solverinfo.converged
                    out = re
                    s = f*s
                    x = vec(reshape(Array(out.nodeflux), out.Nmodes, :))
                    reconverged = true
                    break
                end
            end
            reconverged || error(
                lazy"the carried point did not reconverge on the finest grid $(grids[gi]) even at half its drive; the truncation boundaries of the ladder are too far apart. Add an intermediate grid.")
            pendinggrow = false
            ds = s0/2
        else
            error(lazy"no harmonic balance solution was found at the requested drive: the source continuation converged at $(round(s, digits = 4)) of the requested amplitudes on the finest grid $(grids[end]) and stalled at $(round(starget, digits = 4)). The solution branch ends between them (a fold, the self oscillation threshold); the operating point does not exist at full drive.")
        end
    end
    # the returned diagnostics are the whole walk: one StagedStageInfo per
    # attempt, each carrying its inner solver records, in place of just the
    # final solve's stages
    si = out.solverinfo
    F0 = try
        stagerecords[1].inner[1].normresidual[1]
    catch
        si.initialresidual
    end
    newsi = SolverInfo(stagerecords, F0, si.finalresidual, si.converged,
        si.sourcefold)
    vals = Any[getfield(out, f) for f in fieldnames(typeof(out))]
    vals[findfirst(==(:solverinfo), fieldnames(typeof(out)))] = newsi
    out = typeof(out)(vals...)
    return out
end

"""
    hbnlsolve(w::NTuple{N,Number}, Nharmonics::NTuple{N,Int}, sources,
        circuit, circuitdefs; iterations = 1000,
        maxintermodorder = Inf, dc = false, odd = true, even = false,
        x0 = nothing, ftol = 1e-8, symfreqvar = nothing, sorting= :number)

Harmonic balance solver supporting an arbitrary number of large signals
(strong tones or pumps) and arbitrary numbers of ports, sources, and drives
including direct current (zero frequency) or flux pumping using a current
source and a mutual inductor. Use `hblinsolve` to linearize the system of
equations about the operating point found with `hbnlsolve`.

The system is solved in a modified nodal analysis (MNA) formulation in the
node flux basis: resistors with constant real values (including complex
storage with zero imaginary part) and mutually coupled inductor branches
are assigned auxiliary branch current variables with their constitutive
relations kept as explicit equations, which is algebraically equivalent to
the nodal formulation wherever the latter is well posed. Promoting the
coupled inductors keeps the system matrix entries bounded as the coupling
coefficient approaches one, where the nodal inverse inductance entries
diverge as `1/(1-k^2)`. The system is nondimensionalized by the solver
inductance scale `Z0/w0` (see [`calcsolverscale`](@ref)), the geometric
mean port impedance over the geometric mean nonzero drive frequency, so
the residual tolerance `ftol` is independent of the unit system and the
auxiliary variables have magnitudes comparable to the node fluxes. One gauge fixing equation per floating
inductive/Josephson subnetwork and zero-frequency mode makes circuits
which are structurally singular at DC in a purely nodal formulation
(nodes or subnetworks with no inductive path to ground) exactly solvable
without workaround inductors; if the net direct current injected into
such a subnetwork is nonzero, no periodic solution exists and an
`ArgumentError` is thrown. The reported residual norms are those of the
augmented system, a supplied `x0` may have either the node flux length or
the augmented length, and the returned structure contains only the node
fluxes and the original incidence matrix. Commensurate drive frequencies
whose retained intermodulation products reach (numerically) zero
frequency are rejected with an `ArgumentError`. See `src/mna.jl`.
# Arguments
- `w::NTuple{N,Number}`: a tuple containing the angular frequencies of the
    strong tones (or pumps) such as (2\\*pi\\*5.0e9,) for a single tone at 5
    GHz and (2\\*pi\\*5.0e9,2\\*pi\\*6.0e9) for a tone at 5 GHz and a tone at
    6 GHz. The frequencies should be non-commensurate. For commensurate
    frequencies, the lowest frequency should be provided here, and the other
    added to `sources` with a mode index equal to the ratio.
- `Nharmonics::NTuple{N,Int}`: a tuple of integers describing how many
    harmonics to simulate for each of the tones. The length of the tuple must
    equal the number of non-commensurate tones.
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
- `circuit`: vector of tuples each of which contain the component name, the
    first node, the second node, and the component value. The first three must
    be strings.
- `circuitdefs`: a dictionary where the keys are symbols or symbolic
    variables for component values and the values are the numerical values
    for the components.

# Keywords
- `iterations = 1000`: the number of iterations before the nonlinear solver
    returns an error.
- `maxintermodorder = Inf`: the maximum intermod order as defined by the sum of
    the absolute values of the integers multiplying each of the frequencies
    being less than or equal to `maxintermodorder`. This performs a diamond
    truncation of the discrete Fourier space.
- `dc = false`: include 0 frequency terms in the harmonic balance analysis.
- `odd = true`: include odd terms in the harmonic balance analysis.
- `even = false`: include even terms in the harmonic balance analysis.
- `x0 = nothing`: initial value for the nodeflux.
- `ftol = 1e-8`: the residual tolerance `norm(F) <= ftol` at which the
    nonlinear solution is considered converged. `F` is  is scaled by `Z0/w0`.
    Aee [`calcsolverscale`](@ref)).
- `method = :newtonkrylov`: the nonlinear solver. `:quasinewton`
    uses the complex holomorphic Jacobian `Jx` only, an approximation to the
    full Jacobian. `:newton` solves the equivalent real system with the full
    Jacobian. `:newtonkrylov` uses the matrix-free real Jacobian.
    `:staged` runs [`stagedhbnlsolve`](@ref): source continuation on an
    adaptively grown harmonic grid, spending the many Newton iterations a
    near critical drive needs on small cheap truncations and warm starting
    each larger grid from the last, with `:newtonkrylov` solving every
    stage. It is the strategy for operating points the direct methods fail
    outright, and the one that distinguishes a hard operating point from a
    nonexistent one: a stall on the finest grid with the minimum drive step
    brackets a fold -- the self oscillation threshold -- and errors with
    the bracketed drive fraction rather than returning garbage.
- `stagedkwargs::NamedTuple = (;)`: options for `method = :staged`; see
    [`stagedhbnlsolve`](@ref) (`grids`, `s0`, `smin`, `interiorftol`,
    `interioriterations`, `innermethod`, `verbose`).
- `andersondepth::Integer = method == :quasinewton ? 5 : 0`: the depth of the
    Anderson acceleration of the Newton fixed point iteration, the maximum
    number of previous iterates used for the extrapolation. Values less than
    one disable the acceleration.
- `symfreqvar = nothing`: the symbolic frequency variable, eg `w`.
- `sorting = :number`: sort the nodes by:
    `:name`: Sort the vector of strings. This always works but leads
    to results like "101" comes before "11".
    `:number`: Convert the node strings to integer and sort by these
    (this errors if the nodes names cannot be converted to integers).
    `:none`: Don't perform any sorting except to place the ground node
    first. In other words, order the nodes in the order they are found in
    `circuit`.

# Returns
- `NonlinearHB`: A simple structure to hold the harmonic balance solutions.
    See [`NonlinearHB`](@ref).

# Examples
```jldoctest
circuit = Tuple{String,String,String,Union{Complex{Float64},Symbol,Int64}}[]
push!(circuit,("P1","1","0",1))
push!(circuit,("R1","1","0",:Rleft))
push!(circuit,("L1","1","0",:Lm)) 
push!(circuit,("K1","L1","L2",:K1))
push!(circuit,("C1","1","2",:Cc)) 
push!(circuit,("L2","2","3",:Lm)) 
push!(circuit,("Lj3","3","0",:Lj)) 
push!(circuit,("Lj4","2","0",:Lj)) 
push!(circuit,("C2","2","0",:Cj))
circuitdefs = Dict{Symbol,Complex{Float64}}(
    :Lj =>2000e-12,
    :Lm =>10e-12,
    :Cc => 200.0e-15,
    :Cj => 900e-15,
    :Rleft => 50.0,
    :Rright => 50.0,
    :K1 => 0.9,
)

Idc = 50e-5
Ip=0.0001e-6
wp=2*pi*5e9
Npumpmodes = 2
out=hbnlsolve(
    (wp,),
    (Npumpmodes,),
    [
        (mode=(0,),port=1,current=Idc),
        (mode=(1,),port=1,current=Ip),
    ],
    circuit,circuitdefs;dc=true,odd=true,even=false)
isapprox(out.nodeflux[:],
    ComplexF64[15.190314040027522 - 8.56492651167657e-24im, 2.991103820177504e-6 - 1.8501001011477133e-8im, -6.835392148510984 - 1.0356102442254259e-14im, 7.396422335315908e-6 - 4.5749403967992827e-8im, 6.835392148539885 - 1.0356102451770844e-14im, 1.008026285172782e-5 - 6.23498762664213e-8im],
    atol = 1e-6)

# output
true
```
"""
function hbnlsolve(w::NTuple{N,Number}, Nharmonics::NTuple{N,Int}, sources,
    circuit, circuitdefs; iterations = 1000,
    maxharmonics::NTuple{N,Int} = Nharmonics,
    maxintermodorder = Inf, dc::Bool = false, odd::Bool = true,
    even::Bool = false, x0 = nothing, ftol = 1e-8,
    switchofflinesearchtol = nothing, alphamin = nothing,
    method = :newtonkrylov, andersondepth::Integer = method == :quasinewton ? 5 : 0,
    symfreqvar = nothing, sorting = :number, keyedarrays::Bool = true,
    sensitivitynames::Vector{String} = String[],
    returnoperatingpoint::Bool = false,
    krylovcouplingmodes = :none,
    krylovrecycle::Integer = 0, krylovharvest::Integer = 8,
    krylovkwargs::NamedTuple = (;),
    stagedkwargs::NamedTuple = (;),
    factorization = nothing, backend = CPU(),
    precision::Type{<:AbstractFloat} = Float64, debugJacobian = false,
    ) where {N}

    if method === :staged
        return stagedhbnlsolve(w, Nharmonics, sources, circuit, circuitdefs;
            iterations = iterations, maxharmonics = maxharmonics,
            maxintermodorder = maxintermodorder, dc = dc, odd = odd,
            even = even, ftol = ftol, symfreqvar = symfreqvar,
            sorting = sorting, keyedarrays = keyedarrays,
            sensitivitynames = sensitivitynames,
            returnoperatingpoint = returnoperatingpoint,
            krylovcouplingmodes = krylovcouplingmodes,
            krylovrecycle = krylovrecycle, krylovharvest = krylovharvest,
            krylovkwargs = krylovkwargs, factorization = factorization,
            backend = backend, precision = precision, stagedkwargs...)
    end

    # calculate the frequency struct
    freq = removeconjfreqs(
        truncfreqs(
            calcfreqsrdft(Nharmonics); dc = dc, odd = odd, even = even,
            maxintermodorder = maxintermodorder, maxharmonics = maxharmonics,
        )
    )

    indices = fourierindices(freq)

    Nmodes = length(freq.modes)

    # parse and sort the circuit
    psc = parsesortcircuit(circuit, sorting = sorting)

    # calculate the circuit graph
    cg = calccircuitgraph(psc)

    # calculate the numeric matrices
    nm = numericmatrices(psc, cg, circuitdefs, Nmodes = Nmodes)


    return hbnlsolve(w, sources, freq, indices, psc, cg, nm;
        iterations = iterations, x0 = x0, ftol = ftol,
        switchofflinesearchtol = switchofflinesearchtol, alphamin = alphamin,
        method = method, andersondepth = andersondepth,
        symfreqvar = symfreqvar, keyedarrays = keyedarrays,
        sensitivitynames = sensitivitynames,
        returnoperatingpoint = returnoperatingpoint,
        krylovcouplingmodes = krylovcouplingmodes,
        krylovrecycle = krylovrecycle, krylovharvest = krylovharvest,
        krylovkwargs = krylovkwargs,
        factorization = factorization, backend = backend,
        precision = precision, debugJacobian = debugJacobian,
        )
end

"""
    hbnlsolve(w::NTuple{N,Number}, sources, frequencies::Frequencies{N},
        indices::FourierIndices{N}, psc::ParsedSortedCircuit, cg::CircuitGraph,
        nm::CircuitMatrices; iterations = 1000, x0 = nothing,
        ftol = 1e-8, symfreqvar = nothing)

New version of the nonlinear harmonic balance solver suitable for arbitrary
numbers of ports, sources, and drives including direct current (zero frequency)
or flux pumping using a current source and a mutual inductor.

# Examples
```jldoctest
circuit = Tuple{String,String,String,Union{Complex{Float64},Symbol,Int64}}[]
push!(circuit,("P1","1","0",1))
push!(circuit,("R1","1","0",:Rleft))
push!(circuit,("L1","1","0",:Lm)) 
push!(circuit,("K1","L1","L2",:K1))
push!(circuit,("C1","1","2",:Cc)) 
push!(circuit,("L2","2","3",:Lm)) 
push!(circuit,("Lj3","3","0",:Lj)) 
push!(circuit,("Lj4","2","0",:Lj)) 
push!(circuit,("C2","2","0",:Cj))
circuitdefs = Dict{Symbol,Complex{Float64}}(
    :Lj =>2000e-12,
    :Lm =>10e-12,
    :Cc => 200.0e-15,
    :Cj => 900e-15,
    :Rleft => 50.0,
    :Rright => 50.0,
    :K1 => 0.9,
)

Idc = 50e-5
Ip=0.0001e-6
wp=2*pi*5e9
Nharmonics = (2,)
frequencies = JosephsonCircuits.removeconjfreqs(
    JosephsonCircuits.truncfreqs(
        JosephsonCircuits.calcfreqsrdft(Nharmonics),
        dc=true, odd=true, even=false, maxintermodorder=Inf,
    )
)
fi = JosephsonCircuits.fourierindices(frequencies)
Nmodes = length(frequencies.modes)
psc = JosephsonCircuits.parsesortcircuit(circuit)
cg = JosephsonCircuits.calccircuitgraph(psc)
nm = JosephsonCircuits.numericmatrices(psc, cg, circuitdefs, Nmodes = Nmodes)

out=hbnlsolve(
    (wp,),
    [
        (mode=(0,),port=1,current=Idc),
        (mode=(1,),port=1,current=Ip),
    ],
    frequencies, fi, psc, cg, nm)
isapprox(out.nodeflux[:],
    ComplexF64[15.190314040027522 - 8.56492651167657e-24im, 2.991103820177504e-6 - 1.8501001011477133e-8im, -6.835392148510984 - 1.0356102442254259e-14im, 7.396422335315908e-6 - 4.5749403967992827e-8im, 6.835392148539885 - 1.0356102451770844e-14im, 1.008026285172782e-5 - 6.23498762664213e-8im],
    atol = 1e-6)

# output
true
```

The system is solved in a modified nodal analysis formulation in the node
flux basis; see the primary [`hbnlsolve`](@ref) docstring and `src/mna.jl`
for details.
"""
function hbnlsolve(w, sources, frequencies::Frequencies,
    indices::FourierIndices, psc::ParsedSortedCircuit, cg::CircuitGraph,
    nm::CircuitMatrices; iterations = 1000, x0 = nothing,
    ftol = 1e-8, switchofflinesearchtol = nothing, alphamin = nothing,
    method = :newtonkrylov, andersondepth::Integer = method == :quasinewton ? 5 : 0,
    symfreqvar = nothing, keyedarrays::Bool = true,
    sensitivitynames::Vector{String} = String[],
    returnoperatingpoint::Bool = false,
    krylovcouplingmodes = :none,
    krylovrecycle::Integer = 0, krylovharvest::Integer = 8,
    krylovkwargs::NamedTuple = (;),
    factorization = nothing, backend = CPU(),
    precision::Type{<:AbstractFloat} = Float64, debugJacobian = false,
    )

    # the factorization defaults to the one matching the backend: KLU on the
    # host, cuDSS on a device, where a host factorization could not be applied
    # to the device vectors of the Krylov iteration anyway. An explicit
    # factorization is always honored.
    factorization = if !isnothing(factorization)
        factorization
    elseif backend isa CPU
        KLUfactorization()
    else
        CUDSSFactorization()
    end

    # deprecation warnings for switchofflinesearchtol and alphamin.
    if !isnothing(switchofflinesearchtol)
        Base.depwarn(lazy"The `switchofflinesearchtol` kwarg is deprecated and no longer used (and no longer necessary). Please remove it to avoid errors in future versions.", :hbnlsolve; force=true)
    end

    if !isnothing(alphamin)
        Base.depwarn(lazy"The `alphamin` kwarg is deprecated and no longer used (and no longer necessary). Please remove it to avoid errors in future versions.", :hbnlsolve; force=true)
    end

    # reject non-finite physical inputs before mode construction and
    # source assembly: the frequency canonicalization and source
    # compatibility diagnostics compare against bounds built from these
    # values, and non-finite inputs would otherwise reach them (an
    # infinite value compares as within an infinite bound, and a NaN
    # silently fails every comparison) or fail later without explanation.
    all(isfinite, w) || throw(ArgumentError("All drive frequencies must be finite."))
    for source in sources
        if !isfinite(abs(source[:current]))
            throw(ArgumentError("All source currents must be finite."))
        end
    end

    Nharmonics = frequencies.Nharmonics
    Nw = frequencies.Nw
    Nt = frequencies.Nt
    coords = frequencies.coords
    modes = frequencies.modes

    conjsymdict = indices.conjsymdict
    freqindexmap = indices.vectomatmap
    conjsourceindices = indices.conjsourceindices
    conjtargetindices = indices.conjtargetindices
    Amatrixmodes = indices.hbmatmodes
    Amatrixindices = indices.hbmatindices
    Amatrixconjindices = indices.hbconjmatindices
    # The harmonic balance residual is computed with cyclic Fourier
    # transforms, so its exact Jacobian couples modes whose differences
    # alias back onto the sampled grid; hbmatind with `alias = true`
    # includes those couplings (the sum couplings of hbconjmatind always
    # alias). The exact real Jacobian of method = :newton is assembled from
    # the aliased difference indices, so it is the exact derivative of the
    # residual for multi-tone problems as well; this was established with
    # the matrix-free Jacobian-vector products of HBSystem as the ground
    # truth, which the assembled Jacobian matches to machine precision. The
    # complex holomorphic Jacobian of method = :quasinewton deliberately
    # keeps the truncated (non-aliased) indices: it is an approximation to
    # the exact Jacobian either way, the truncation does not change its
    # iteration counts in practice, and the aliased couplings would densify
    # it and slow its factorization.
    Amatrixindicesaliased = hbmatind(frequencies; alias = true)[2]

    # generate the frequencies of the modes
    Nmodes = length(modes)
    wmodes = calcmodefreqs(w,modes)

    # only the all-zero mode tuple may represent zero frequency: its
    # frequency is exactly 0.0 by construction, which makes the exact
    # iszero classification in the gauge fixing and source compatibility
    # machinery sound. any other retained tuple whose physical frequency
    # cancels to zero, or to within roundoff of the magnitudes being
    # combined (commensurate drive frequencies), would duplicate the DC
    # coordinate with different conjugacy assumptions and produce
    # vanishing capacitor and resistor stamps, so it is rejected with an
    # explanation rather than allowed to form a (nearly) singular system.
    for m in 1:Nmodes
        if any(!iszero, modes[m])
            terms = [float(real(modes[m][j]*w[j])) for j in eachindex(w)]
            if isnumericallyzero(wmodes[m], terms)
                throw(ArgumentError("The mode tuple $(modes[m]) has a "*
                    "(numerically) zero physical frequency for the drive "*
                    "frequencies $(w./(2*pi)) Hz, so it would duplicate "*
                    "the DC coordinate: the drive frequencies are "*
                    "commensurate at the retained intermodulation order. "*
                    "Choose incommensurate drive frequencies (for example "*
                    "by a small offset) or reduce the number of "*
                    "harmonics."))
            end
        end
    end

    # extract the elements we need
    Nnodes = psc.Nnodes
    componentnames = psc.componentnames
    componentnamedict = psc.componentnamedict
    componenttypes = psc.componenttypes
    nodenames = psc.nodenames
    nodeindices = psc.nodeindices
    Nbranches = cg.Nbranches
    edge2indexdict = cg.edge2indexdict
    Ljb = nm.Ljb
    Ljbm = nm.Ljbm
    Rbnm = nm.Rbnm
    Cnmcopy = nm.Cnm
    Gnmcopy = nm.Gnm
    invLnmcopy = nm.invLnm
    portindices = nm.portindices
    portnumbers = nm.portnumbers
    portimpedanceindices = nm.portimpedanceindices
    vvn = nm.vvn
    # if there are no inductors, then Lmean will be zero so set it to be one
    Lmean = if iszero(nm.Lmean)
        one(eltype(nm.Lmean))
    else
        nm.Lmean
    end
    # fail immediately, with the actual cause, if any component value
    # contains symbolic variables which were not assigned numerical values
    # in circuitdefs (values depending only on the symbolic frequency
    # variable are frequency dependent components and are accepted).
    checkcomponentvaluesdefined(componentnames, vvn, symfreqvar)

    # nondimensionalize the system with the solver inductance scale Z0/w0
    # (see calcsolverscale) instead of the mean inductance: the scaled
    # entries are then of order one for circuits driven near their
    # characteristic impedance and frequency, the auxiliary branch currents
    # of the modified nodal analysis formulation have magnitudes comparable
    # to the node fluxes even in circuits without inductors, and the
    # residual tolerance ftol is independent of the unit system. The scale
    # multiplies rows only, so the node fluxes and all physical outputs are
    # unchanged. The local name Lmean is kept for continuity with the
    # matrix and plan interfaces.
    Lmean = calcsolverscale(w, componenttypes, vvn, portimpedanceindices,
        Lmean)
    Lb = nm.Lb

    # find the indices associated with the components for which we will
    # calculate sensitivities
    sensitivityindices = zeros(Int,length(sensitivitynames))
    for i in eachindex(sensitivitynames)
        sensitivityindices[i] = componentnamedict[sensitivitynames[i]]
        if componenttypes[sensitivityindices[i]] == :S
            throw(ArgumentError(lazy"Sensitivities with respect to scattering block components are not supported; got $(sensitivitynames[i])."))
        end
    end

    # calculate the diagonal frequency matrices
    wmodesm = Diagonal(repeat(wmodes, outer = Nnodes-1))
    wmodes2m = Diagonal(repeat(wmodes.^2, outer = Nnodes-1))

    # calculate the source terms in the branch basis
    bbm = calcsources(modes, sources, portindices, portnumbers,
        nodeindices, edge2indexdict, Lmean, Nnodes, Nbranches, Nmodes)

    # convert from the node basis to the branch basis
    bnm = transpose(Rbnm)*bbm

    # calculate the dimensions of the array which holds the frequency
    # domain information for the fourier transform
    Nwtuple = NTuple{length(Nw)+1,Int}((Nw..., length(Ljb.nzval)))

    # create an array to hold the frequency domain data for the
    # fourier transform
    # on the backend: the time domain array and every workspace of the system
    # derive from it with similar, so they follow it there.
    phimatrix = tobackend(backend, zeros(Complex{precision}, Nwtuple))

    # create an array to hold the time domain data for the RFFT. also generate
    # the plans, which come from the package extension on a device backend.
    phimatrixtd, irfftplan, rfftplan = plan_applynl(phimatrix, backend)

    # the number of frequency entries per Josephson junction in phimatrix
    Nfreq = prod(Nwtuple[1:end-1])

    x = if isnothing(x0)
        zeros(Complex{Float64}, (Nnodes-1)*Nmodes)
    else
        copy(x0)
    end
    F = zeros(Complex{Float64}, (Nnodes-1)*Nmodes)
    # substitute in the mode frequencies for components which have frequency
    # defined symbolically.
    Cnm = freqsubst(Cnmcopy, wmodes, symfreqvar)
    Gnm = freqsubst(Gnmcopy, wmodes, symfreqvar)
    # Inductor branches which participate in mutual coupling are promoted
    # to auxiliary branch current variables by the modified nodal analysis
    # formulation below, instead of being eliminated through the inverse of
    # the branch inductance matrix: the inverse inductance entries of the
    # nodal formulation diverge as 1/(1-k^2) as the coupling coefficient k
    # approaches one (and do not exist at |k| = 1), while the branch
    # inductance entries of the promoted constitutive equations remain
    # bounded. numericmatrices already excludes the coupled branches from
    # the inverse inductance matrix.
    Mb = nm.Mb
    coupledbranches = mnacoupledbranches(Mb)
    invLnm = freqsubst(invLnmcopy, wmodes, symfreqvar)


    # take the complex conjugate of the terms associated with modes with
    # negative frequencies. this is the same operation hblinsolve performs
    # with sparseaddconjsubst!. these matrices are used in both the residual
    # and the Jacobian in calcfj2!.
    conjnegfreq!(Cnm, wmodes)
    conjnegfreq!(Gnm, wmodes)
    conjnegfreq!(invLnm, wmodes)

    # scale the matrices for numerical reasons
    rmul!(Cnm,Lmean)
    rmul!(Gnm,Lmean)
    rmul!(invLnm,Lmean)


    # Set up the modified nodal analysis (MNA) formulation, which assigns
    # auxiliary branch current variables to the port resistors, keeping their
    # constitutive relations as explicit equations, and adds one gauge
    # fixing equation per floating component of the static flux-stiffness
    # graph and zero-frequency mode. Eliminating the auxiliary variables
    # recovers the nodal equations exactly, so the results are identical
    # whenever the nodal formulation is well posed. Because the MNA
    # equations are linear, the whole augmentation reduces to padding the
    # incidence matrix with empty columns and folding the constant matrix
    # Amna into the static linear term, after which the plan machinery,
    # HBSystem, and the solvers operate on the augmented system unchanged.
    # See mna.jl for details.
    Rbnmout = Rbnm
    Nnodal = (Nnodes-1)*Nmodes
    # reject inductor and junction values for which the static
    # classification and the matrix stamps are ill defined
    checkstaticstiffnessvalues(componenttypes, vvn)
    mnaindices = mnaportresistorindices(componenttypes, nodeindices,
        psc.mutualinductorbranchnames, vvn)
    Nauxr = length(mnaindices)*Nmodes
    Nauxscattering = countscatteringports(psc)*Nmodes
    Naux = Nauxr + length(coupledbranches)*Nmodes + Nauxscattering
    # remove the promoted resistors from the conductance matrix, applying
    # the same conjugation and scaling as for Gnm above.
    if !isempty(mnaindices)
        Gnmp = calcGn(componenttypes[mnaindices],
            nodeindices[:, mnaindices], vvn[mnaindices], Nmodes, Nnodes)
        conjnegfreq!(Gnmp, wmodes)
        rmul!(Gnmp, Lmean)
        Gnm = mnasubtractpromoted(Gnm, Gnmp)
    end
    # gauge fixing is based on the connected components of the static
    # flux-stiffness graph (edges are inductors and Josephson junctions),
    # which handles floating inductive and Josephson subnetworks as well
    # as individually isolated nodes. before adding the gauge equations,
    # check that the net direct current injected into each floating
    # component is zero, because otherwise no periodic solution exists and
    # the gauge equation would silently absorb the incompatible source
    # into the flux reference.
    floatingcomponents = calcstaticfluxcomponents(componenttypes,
        nodeindices, vvn, Nnodes)
    checkdcsourcecompatibility(floatingcomponents, bnm, wmodes, Nmodes,
        nodenames)
    gaugeindices = calcdcgaugeindices(floatingcomponents, wmodes, Nmodes)
    Amna = calcAmna(mnaindices, nodeindices, vvn, gaugeindices, wmodes,
        Nmodes, Nnodes, Lmean)
    # pad with the coupled inductor auxiliary variables and add their
    # constitutive equations and Kirchhoff current law couplings, which are
    # real and frequency independent (see calcAmnaind).
    Amna = mnapad(Amna, length(coupledbranches)*Nmodes + Nauxscattering)
    AmnaL = calcAmnaind(coupledbranches, Lb, Mb, cg.Rbn, Nmodes,
        Nnodal + Nauxr, Nnodal + Naux, Lmean)
    Amna = spaddkeepzeros(Amna, AmnaL)
    # the hybrid (wave to modified nodal analysis) contribution of the
    # scattering block components: the pump mode frequencies are fixed, so
    # the constitutive equations im*w_m*Lmean*B(w_m)*phi - C(w_m)*i = 0 and
    # the Kirchhoff current law couplings of the auxiliary port currents
    # form a constant matrix, folded into the augmentation exactly like the
    # promoted resistor equations (same signed frequency conjugation
    # convention and Lmean scaling; see scatteringlinearterm). The
    # residual, Jacobian, and solver machinery operate on the augmented
    # system unchanged.
    if Nauxscattering > 0
        # The pump mode frequencies are fixed, so this is a constant matrix
        # folded into the augmentation before the solve begins, exactly as
        # the promoted resistor equations are. Nothing about it is per
        # iteration, so it needs nothing of the backend: the augmented
        # system it produces is what every backend already solves.
        Amna = spaddkeepzeros(Amna, scatteringlinearterm(psc, wmodes,
            Nmodes; auxoffset = Nnodal + Naux - Nauxscattering,
            Ntotal = Nnodal + Naux, scale = Lmean))
    end
    # pad the system matrices and vectors with the auxiliary variables.
    # the incidence matrix gains empty columns so the branch fluxes and
    # the Josephson junction terms are unaffected, and the constant MNA
    # equations are folded into the static linear term, which is added
    # with unit coefficient and is its own contribution to the Jacobian.
    Rbnm = hcat(Rbnm, spzeros(eltype(Rbnm), size(Rbnm,1), Naux))
    invLnm = spaddkeepzeros(mnapad(invLnm, Naux), Amna)
    Cnm = mnapad(Cnm, Naux)
    Gnm = mnapad(Gnm, Naux)
    bnm = vcat(bnm, zeros(eltype(bnm), Naux))
    wmodesm = Diagonal(repeat(wmodes,
        outer = Nnodes-1+length(mnaindices)+length(coupledbranches)+
            countscatteringports(psc)))
    wmodes2m = Diagonal(repeat(wmodes.^2,
        outer = Nnodes-1+length(mnaindices)+length(coupledbranches)+
            countscatteringports(psc)))
    if length(x) == Nnodal
        # accept a nodal initial guess, materializing keyed arrays or
        # other array types into a plain vector. transform it into the
        # selected gauge (a physically irrelevant common DC shift of a
        # floating component would otherwise enter the gauge fixing rows)
        # and initialize the auxiliary currents consistently with the node
        # fluxes, so the initial augmented residual equals the initial
        # nodal Kirchhoff current law residual.
        x = vcat(Vector{Complex{Float64}}(vec(x)),
            zeros(Complex{Float64}, Naux))
        mnagaugenormalize!(x, floatingcomponents, wmodes, Nmodes)
        mnainitialauxall!(x, Amna, Nnodal, Nauxr, coupledbranches, Lb, Mb,
            cg.Rbn, Nmodes, Lmean)
    elseif length(x) == Nnodal + Naux
        # accept a full augmented state, with the layout documented in
        # calcAmna, applying the same gauge normalization and reconciling
        # the auxiliary currents with the constitutive relations.
        x = Vector{Complex{Float64}}(vec(x))
        mnagaugenormalize!(x, floatingcomponents, wmodes, Nmodes)
        mnainitialauxall!(x, Amna, Nnodal, Nauxr, coupledbranches, Lb, Mb,
            cg.Rbn, Nmodes, Lmean)
    else
        throw(DimensionMismatch(lazy"The initial value x0 has length $(length(x)) but the solver expects $(Nnodal) node flux unknowns, optionally followed by $(Naux) auxiliary current unknowns."))
    end
    F = zeros(Complex{Float64}, Nnodal + Naux)

    modelayout = ModeLayout(selfconjmodes(frequencies), Nnodal + Naux)
    Fr = complex_to_real(F, modelayout.isreal)
    xr = complex_to_real(x, modelayout.isreal)

    # Build the Jacobian matrix and precomputed assembly plan for the chosen
    # method. The plans map the Fourier coefficients of cos(phi(t)) and the
    # frequency dependent linear term matrices directly into the nonzeros of
    # the Jacobian, so no complex branch matrices, incidence matrix
    # multiplications, or real conversions are performed while iterating.
    # Only the machinery for the chosen method is built, unless the debug
    # output is requested, in which case both are built so they can be
    # compared.
    # the complex Jacobian, as :newton's real one below: on a backend the
    # structure is stored transposed, because a device factorization is
    # compressed by rows, and the values live on the backend too. The host
    # `Jx` is still built either way, because `debugJacobian` compares
    # against it and the sensitivity calculation reads its structure.
    Jx, complexjacobianplan = if method == :quasinewton || debugJacobian
        plancomplexjacobian(Amatrixindices, Ljb, Lmean, Rbnm, Nmodes,
            Nbranches, Nfreq, invLnm, Gnm, Cnm)
    else
        nothing, nothing
    end
    devicex = method == :quasinewton && !(backend isa CPU)
    Jxb, complexjacobianplan = if devicex
        Jxt = sparse(transpose(Jx))
        nodesandsignsx = branchnodesandsigns(Rbnm, Nmodes, Nbranches)
        cjp = planstructurecomplexjacobian(Jxt, Float64, Amatrixindices, Ljb,
            Lmean, nodesandsignsx, invLnm, Gnm, Cnm, wmodesm, wmodes2m,
            Nmodes, Nfreq, backend; transposed = true)
        # the pattern goes to the backend as well, so the sparse products of
        # the line search read it there rather than from the host
        patt = DeviceSparsePattern(
            tobackend(backend, SparseArrays.getcolptr(Jxt)),
            tobackend(backend, rowvals(Jxt)), size(Jxt, 1), size(Jxt, 2))
        DeviceValuedSparseMatrix(patt,
            tobackend(backend, zeros(eltype(Jx), nnz(Jxt)))), cjp
    else
        Jx, complexjacobianplan
    end

    # method = :newtonkrylov deliberately does not appear here. Its steps come
    # from the matrix-free Jacobian-vector product and its preconditioner
    # builds and assembles its own, much sparser, restricted plan, so neither
    # the full Jacobian nor its plan is ever used. On a multi-tone problem
    # that plan is the single largest object in the solve (millions of stored
    # entries), so not building it is a substantial saving in both time and
    # memory.
    Jr, realjacobianplan = if method == :newton || debugJacobian ||
            returnoperatingpoint
        # on a backend the structure is built there, transposed, because a
        # device factorization is compressed by rows, and the assembly writes
        # its values there too; on a host it is the matrix KLU factorizes
        devicej = !(backend isa CPU)
        Jrs, nodesandsigns = realjacobianstructure(Amatrixindicesaliased,
            Amatrixconjindices, Ljb, Rbnm, Nmodes, Nbranches, invLnm, Gnm,
            Cnm, modelayout, modelayout, Float64;
            transposed = devicej, backend = backend)
        rjp = planstructurerealjacobian(Jrs, Float64, Amatrixindicesaliased,
            Amatrixconjindices, Ljb, Lmean, nodesandsigns, invLnm, Gnm, Cnm,
            wmodesm, wmodes2m, modelayout, modelayout, Nmodes, Nfreq,
            backend; transposed = devicej)
        (devicej ? DeviceValuedSparseMatrix(Jrs,
            tobackend(backend, zeros(Float64, nnz(Jrs)))) : Jrs), rjp
    else
        nothing, nothing
    end

    # the unified evaluation object for the nonlinear system: the residual,
    # the matrix-free Jacobian-vector and Hessian-vector products, and the
    # assembled Jacobians are all evaluated through it, in the complex or
    # the equivalent real representation.
    # the real form of the linear term is needed exactly when the solve
    # applies the real representation entry points. That is stated directly
    # rather than derived from realjacobianplan: :newtonkrylov solves in the
    # real representation but takes its preconditioner from a
    # ModeCouplingPreconditioner, which builds its own restricted plan, so it
    # needs no full realjacobianplan and the two conditions are not the same.
    # method = :quasinewton stays in the complex representation and skips it.
    realrepresentation = method == :newton || method == :newtonkrylov ||
        debugJacobian || returnoperatingpoint
    sys = HBSystem(Rbnm, invLnm, Gnm, Cnm, wmodesm, wmodes2m, bnm,
        Ljb, Ljbm, Lmean, Nbranches, freqindexmap, conjsourceindices,
        conjtargetindices, phimatrix, phimatrixtd, irfftplan, rfftplan,
        modelayout, realjacobianplan, complexjacobianplan, backend;
        realbackward = realrepresentation)

    # the residual and the complex (holomorphic) Jacobian, an approximation
    # to the exact Jacobian, for method == :quasinewton
    function fj!(F, Jx, x)
        setpoint!(sys, x)
        isnothing(F) || residual!(F, sys)
        isnothing(Jx) || jacobian!(Jx, sys)
        return nothing
    end

    # the residual and the exact Jacobian of the equivalent real system, for
    # method == :newton
    function fjreal!(Fr, Jr, xr)
        setpoint!(sys, xr)
        isnothing(Fr) || residual!(Fr, sys)
        isnothing(Jr) || jacobian!(Jr, sys)
        # write back the canonical real representation of the point. through
        # the plan rather than complex_to_real!, whose serial cursor over a
        # BitVector is host only; the two agree on CPU().
        applycomplextoreal!(xr, sys.nonlineartermplan, sys.x)
        return nothing
    end

    # diagnostics for each invocation of the nonlinear solver, stored in
    # the output so the solution process can be assessed after the run
    solverstages = IterationInfo[]

    # the source scale at which the source stepping detected a fold in
    # the branch of solutions, if any
    sourcefold = Ref(NaN)

    # use this for debugging purposes to return the residual and Jacobian
    # functions along with the ingredients from which the Jacobians are
    # assembled, so reference implementations can be constructed, eg. in the
    # tests.
    if debugJacobian
        return (F=F, x=x, Fr=Fr, xr=xr, Jx=Jxb, Jr=Jr, fj=fj!, fjreal=fjreal!,
            sys=sys, Nnodal=Nnodal, mnaindices=mnaindices,
            gaugeindices=gaugeindices, floatingcomponents=floatingcomponents,
            coupledbranches=coupledbranches,
            complexjacobianplan=complexjacobianplan,
            realjacobianplan=realjacobianplan, phimatrix=phimatrix,
            cosphimatrix=(x -> (setpoint!(sys, x);
                JosephsonCircuits._updatecosphimatrix!(sys);
                sys.phimatrix)), modelayout=modelayout,
            Amatrixindices=Amatrixindices, Amatrixmodes=Amatrixmodes,
            modes=modes,
            Amatrixindicesaliased=Amatrixindicesaliased,
            Amatrixconjindices=Amatrixconjindices, Ljb=Ljb, Ljbm=Ljbm,
            Lmean=Lmean, Rbnm=Rbnm, invLnm=invLnm, Gnm=Gnm, Cnm=Cnm,
            wmodesm=wmodesm, wmodes2m=wmodes2m, Nmodes=Nmodes,
            Nbranches=Nbranches, Nfreq=Nfreq)
    end

    # solve the nonlinear system
    info = if method == :quasinewton

        solveonbackend!(fj!, F, Jxb, x, backend; iterations = iterations,
            ftol = ftol, andersondepth = andersondepth,
            factorization = factorization)

    elseif method == :newton

        # solve the equivalent real nonlinear system with the exact real
        # Jacobian then convert back to complex
        info = solveonbackend!(fjreal!, Fr, Jr, xr, backend;
            iterations = iterations, ftol = ftol,
            andersondepth = andersondepth, factorization = factorization)
        real_to_complex!(x,xr,modelayout.isreal)
        real_to_complex!(F,Fr,modelayout.isreal)
        info
    elseif method == :newtonkrylov

        # the matrix free real Jacobian
        jvpreal!(Jvr, vr) = jacobianvectorproduct!(Jvr, sys, vr)

        # The right preconditioner is the Jacobian materialized on a
        # coarser frequency grid: the mode coupling is kept in full only on
        # the modes where the linear circuit responds strongly (measured,
        # see modesoftness) and reduced to the mode diagonal elsewhere, so
        # its factorization is a fraction of the full Jacobian's. The setup
        # above was done once, on the user's harmonic selection; the coarse
        # grid exists only inside the preconditioner, so the residual, the
        # Jacobian-vector product and the solution are always those of the
        # requested truncation.
        # `krylovrecycle > 0` wraps the base preconditioner in a recycled
        # deflation subspace (see RecyclingPreconditioner). The intended
        # pairing is with `krylovcouplingmodes = :none`: a batch of small
        # per mode factorizations plus dense level 2 BLAS, which needs no
        # large sparse factorization and is the part of this solver that
        # ports to a GPU.
        #
        # The default is the mode block diagonal, whose factorization is a
        # batch of small independent per mode solves rather than one large
        # sparse factorization. On a strongly pumped line the block diagonal
        # alone stalls; what rescues it is escalation, which grows the base
        # only on repeated linear failures and in practice fires once or
        # twice.
        #
        # `krylovcouplingmodes = :band => p` restricts the retained coupling by
        # harmonic *offset* rather than by column (see `modebandmask`). That is
        # the restriction the Toeplitz structure of the nonlinear term asks for,
        # and at equal fill it is not close: on an eight mode chain driven to
        # max|phi| = 1.9 rad a bandwidth of one converges a Newton path in 118
        # GMRES iterations where a two column selection storing the same number
        # of nonzeros fails to converge in 1051. Its fill grows linearly in the
        # mode count where the full Jacobian's grows quadratically, so it
        # overtakes the full factorization once there are enough modes: measured
        # on a fixed circuit at max|phi| ~ 1.75 rad, the ratio of banded to full
        # total solve time is 1.28 at 8 modes, 0.75 at 16, 0.61 at 24 and 0.31
        # at 32, while the bandwidth that wins stays at one. It escalates by one
        # offset at a time rather than jumping to the full operator. Measured on a two tone line, this is the only configuration
        # whose standing improves with problem size: at 288 cells it is 3.72 s
        # against 8.19 s for mode selection and 5.48 s for the direct solve,
        # having been the slower of the two at 128 cells.
        #
        # `krylovrecycle > 0` additionally wraps the base in a recycled
        # deflation subspace. It also rescues the block diagonal, and needs no
        # sparse factorization at all, which makes it the natural fit for a
        # GPU; but it is a second answer to the same problem and measured a
        # net loss against escalation at 192 cells and above, so it is off by
        # default. `krylovcouplingmodes = 12, krylovrecycle = 0` recovers the
        # earlier frequency based mode selection, which is still the fastest
        # option on moderately sized lines.
        base = ModeCouplingPreconditioner(sys, Amatrixindicesaliased,
            Amatrixconjindices, Ljb, Lmean, Rbnm, Nmodes, Nbranches, Nfreq,
            invLnm, Gnm, Cnm, modelayout;
            couplingmodes = krylovcouplingmodes,
            factorization = factorization, precision = precision,
            Amatrixmodes = Amatrixmodes)

        # the Krylov workspaces are allocated `similar` to the vectors handed
        # in, so handing in device vectors is what puts the whole iteration on
        # the device. On CPU() tobackend adopts them and these are no-ops.
        xrb = tobackend(backend, convert(Vector{precision}, xr))
        Frb = tobackend(backend, convert(Vector{precision}, Fr))

        # built from the vector rather than from its length, so the deflation
        # subspace lives where the iteration does. The `n::Integer` method
        # allocates on the host, which on a device backend leaves the whole
        # recycling path handing host arrays to device kernels.
        pc = if krylovrecycle > 0
            RecyclingPreconditioner(base, jvpreal!, xrb;
                kmax = krylovrecycle, kharvest = krylovharvest)
        else
            base
        end

        # Defaults tuned for this preconditioner, overridable through
        # `krylovkwargs`. The preconditioner is refreshed eagerly, which the
        # generic default does not do: rebuilding the block diagonal costs a
        # fraction of a full factorization, while a preconditioner frozen at
        # zero flux is stale immediately, because at zero flux the Jacobian
        # has no harmonic coupling at all.
        #
        # The restart length cannot be set independently of the
        # preconditioner. A restricted but not diagonal preconditioner leaves
        # a small number of directions which the Krylov space must resolve
        # before a restart discards the progress made on them, and until it
        # can, the iteration count is pinned and does not respond to the
        # preconditioner at all. Measured on a two tone grid of 74 modes, a
        # banded preconditioner needs 1441 iterations without converging at a
        # restart of 30, 60 and 120, and converges in 216 at 240; the block
        # diagonal converges at none of them.
        #
        # The block diagonal has the same problem in a milder form, and it is
        # the reason this solver escalates to the full Jacobian: a short
        # Krylov space cannot resolve the near null directions the block
        # diagonal leaves, the linear solve stalls, and escalation is what
        # rescues it. A long cycle attacks the stall directly and so is the
        # cheaper answer where it works -- the basis is `krylovrestart + 1`
        # vectors of the system dimension, which is cheap next to the sparse
        # factorization escalation would build. The default is therefore a
        # long cycle for every preconditioner, restricted or not. An explicit
        # `krylovkwargs.krylovrestart` always wins.
        krylovdefaults = (; krylovrefreshiterations = 1, krylovrestart = 400)

        info = nlsolvekrylov!(fjreal!, jvpreal!, Frb, xrb, pc;
            iterations = iterations, ftol = ftol,
            merge(krylovdefaults, krylovkwargs)...)
        # back to the host for the complex representation returned to the
        # caller, whose conversion walks a BitVector serially
        copyto!(xr, tohost(xrb))
        copyto!(Fr, tohost(Frb))
        real_to_complex!(x,xr,modelayout.isreal)
        real_to_complex!(F,Fr,modelayout.isreal)
        info
    else
        throw(ArgumentError("Method $(method) is not defined."))
    end


    push!(solverstages, IterationInfo(info.label, 1.0,
        info.regularization, info.converged, info.iterations,
        info.normresidual, info.alpha, info.backtracks,
        info.andersonaccepted, info.krylov))

    if !info.converged
        @warn string(lazy"Solver did not converge.")
    end

    converged = info.converged

    # validate the original, ungauged Kirchhoff current law equations
    # directly by reconstructing their residuals from the augmented
    # residual and the state (the gauge fixing equations add x[g] to the
    # augmented residual of each gauge row g). the acceptance policy,
    # including its block-relative infinity-norm tolerance and its
    # rejection of non-finite residuals, is implemented and unit tested in
    # mnavalidatekcl.
    if converged && !isempty(gaugeindices)
        setpoint!(sys, x)
        # `sys` holds its state wherever the backend put it, while `F` and `x`
        # are host vectors, so the residual is evaluated into a buffer like the
        # system's own state and copied back for the host side validation.
        # Evaluating straight into `F` launches the backward term kernel with a
        # host array argument, which does not compile.
        Fgauge = similar(sys.x)
        residual!(Fgauge, sys)
        copyto!(F, Fgauge)
        kclok, normkcl, kcltol = mnavalidatekcl(F, x, gaugeindices,
            Nnodal, bnm, ftol)
        if !kclok
            @warn "The original (ungauged) Kirchhoff current law equations are violated beyond the solver resolution at the returned solution, indicating that a gauge fixing equation absorbed an incompatibility, such as a net direct current injected into a floating subnetwork, into the arbitrary flux reference. Marking the solution as not converged." normkcl kcltol
            converged = false
        end
    end

    # remove the auxiliary variables so the output contains only the node
    # fluxes
    nodeflux = x[1:Nnodal]

    # the norm of the residual at the initial value, for the diagnostics
    normF0 = info.normresidual[1]
    # the norm of the residual at the returned solution, for the
    # diagnostics
    normFfinal = info.normresidual[end]


    # assemble the stored diagnostics of the solution process
    solverinfo = SolverInfo(solverstages, normF0, normFfinal, converged,
        sourcefold[])

    # calculate the scattering parameters for the pump
    Nports = length(portindices)
    S = zeros(Complex{Float64}, Nports*Nmodes, Nports*Nmodes)
    inputwave = zeros(Complex{Float64}, Nports*Nmodes)
    outputwave = zeros(Complex{Float64}, Nports*Nmodes)
    portimpedances = [vvn[i] for i in portimpedanceindices]
    if !isempty(S)
        calcinputoutput!(inputwave, outputwave, nodeflux,
            bnm[1:Nnodal]/Lmean,
            portimpedanceindices, portimpedanceindices, portimpedances,
            portimpedances, nodeindices, componenttypes, wmodes, symfreqvar)
        calcscatteringmatrix!(S, inputwave, outputwave)
    end

    nodefluxout = if keyedarrays
        nodevariabletokeyed(nodeflux, modes, nodenames)
    else
        nodeflux
    end

    #
    Sout = if !isempty(S) && keyedarrays
        Stokeyed(S, modes, portnumbers, modes, portnumbers)
    else
        S
    end

    # the exact Jacobian of the equivalent real system at the converged
    # solution, and everything else needed to propagate a component
    # perturbation through the operating point with the implicit function
    # theorem. The Jacobian stored during the iteration is from the last
    # Newton step, so it is re-assembled here.
    operatingpoint = if returnoperatingpoint
        # retain the augmentation the solver assembled above (the promoted
        # resistor constitutive equations, the coupled inductor blocks, and
        # the gauge fixing rows, on the full augmented layout) rather than
        # reconstructing it, so the operating point cannot drift from the
        # system which was actually solved. calcresidualsensitivity masks
        # the constitutive equation rows of one promoted resistor out of it;
        # the coupled inductor and gauge entries lie outside that mask.
        # An operating point is a host object, whichever backend solved for
        # it. What differentiates it evaluates transforms of one direction at
        # a time, per component, so it needs a working system and not a
        # snapshot of arrays; a system on a backend would make every one of
        # those a transfer, spread across the callers, invisible to a test
        # suite which has no device. So a solve on a backend gets a host
        # twin, built once here from the same ingredients, and everything
        # downstream is the host code the tests already cover.
        #
        # The same reasoning gives the Jacobian: its values are on the
        # backend behind a structure held transposed, and what reads it are
        # sparse direct factorizations and a structural row mask, so it comes
        # home too. Both are the size of the state and are retained once per
        # pump solve rather than produced per signal frequency; the signal
        # sweep that follows still runs wherever it was asked to.
        opsys = if phimatrix isa Array
            setpoint!(sys, x)
            sys
        else
            hostphimatrix = zeros(eltype(phimatrix), size(phimatrix))
            hostphitd, hostirfftplan, hostrfftplan =
                plan_applynl(hostphimatrix, CPU())
            twin = HBSystem(Rbnm, invLnm, Gnm, Cnm, wmodesm, wmodes2m, bnm,
                Ljb, Ljbm, Lmean, Nbranches, freqindexmap, conjsourceindices,
                conjtargetindices, hostphimatrix, hostphitd, hostirfftplan,
                hostrfftplan, modelayout, nothing, nothing, CPU())
            setpoint!(twin, x)
            twin
        end
        # the Jacobian is the one the solver assembled, not the twin's: it is
        # the augmentation the solve actually used
        setpoint!(sys, x)
        jacobian!(Jr, sys)
        HBOperatingPoint(opsys, copy(x), hostsparse(Jr), modelayout, Nnodal,
            Lmean, wmodes, Amna,
            mnaindices, coupledbranches, Nmodes, Nnodes)
    else
        nothing
    end

    return NonlinearHB(w, frequencies, nodefluxout, Rbnmout, Ljb, Lb, Ljbm,
        Nmodes, Nbranches, nodenames, portnumbers, modes, Sout, solverinfo,
        operatingpoint)

end


"""
    calcsources(modes, sources, portindices, portnumbers, nodeindices,
        edge2indexdict, Lmean, Nnodes, Nbranches, Nmodes)

Calculate the source terms in the branch basis. See also [`addsources!`](@ref).

# Examples
```jldoctest
modes = [(0,), (1,)]
sources = [(mode = (0,), port = 1, current = 0.0005), (mode = (1,), port = 1, current = 1.0e-10)]
portindices = [1]
portnumbers = [1]
nodeindices = [2 2 2 2 0 2 3 4 3 3; 1 1 1 1 0 3 4 1 1 1]
edge2indexdict = Dict((1, 2) => 1, (3, 1) => 2, (1, 3) => 2, (4, 1) => 3, (2, 1) => 1, (1, 4) => 3, (3, 4) => 4, (4, 3) => 4)
Lmean = 1.005e-9 + 0.0im
Nnodes = 4
Nbranches = 4
Nmodes = 2
JosephsonCircuits.calcsources(modes, sources, portindices, portnumbers,
    nodeindices, edge2indexdict, Lmean, Nnodes, Nbranches, Nmodes)

# output
8-element Vector{ComplexF64}:
     1526.863796602709 + 0.0im
 0.0003053727593205418 + 0.0im
                   0.0 + 0.0im
                   0.0 + 0.0im
                   0.0 + 0.0im
                   0.0 + 0.0im
                   0.0 + 0.0im
                   0.0 + 0.0im
```
"""
function calcsources(modes, sources, portindices, portnumbers, nodeindices,
    edge2indexdict, Lmean, Nnodes, Nbranches, Nmodes)

    bbm = zeros(Complex{Float64}, Nbranches*Nmodes)

    addsources!(bbm, modes, sources, portindices, portnumbers,
        nodeindices, edge2indexdict, Lmean, Nnodes, Nbranches,
        Nmodes)

    return bbm
end

"""
    addsources!(bbm, modes, sources, portindices, portnumbers,
        nodeindices, edge2indexdict, Lmean, Nnodes, Nbranches, Nmodes)

Calculate the source terms in the branch basis. Overwrite bbm with the output.
See also [`calcsources`](@ref).
"""
function addsources!(bbm, modes, sources, portindices, portnumbers,
    nodeindices, edge2indexdict, Lmean, Nnodes, Nbranches, Nmodes)

    # fill the vector with zeros
    fill!(bbm,0)

    # make a dictionary of ports
    portdict = Dict{eltype(portnumbers), eltype(portindices)}()
    for i in eachindex(portindices)
        portdict[portnumbers[i]] = portindices[i]
    end

    # make a dictionary of modes
    modedict = Dict{eltype(modes), Int}()
    for i in eachindex(modes)
        modedict[modes[i]] = i
    end

    for source in sources
        # pull out the necessary values from the named tuple
        port = source[:port]
        mode = source[:mode]
        current = source[:current]

        # check if the port is in the dictionary of ports
        if haskey(portdict, port)
            portindex = portdict[port]
            # check if the mode is in the dictionary of modes
            if haskey(modedict, mode)
                # if we find the mode and the port, set that branch in bbm
                # equal to the current scaled by the mean inductance and the
                # flux quantum.
                modeindex = modedict[mode]
                key = (nodeindices[1, portindex], nodeindices[2, portindex])
                bbm[(edge2indexdict[key]-1)*Nmodes+modeindex] += Lmean*current/phi0
            else
                throw(ArgumentError(lazy"Source mode $(mode) not found."))
            end
        else
            throw(ArgumentError(lazy"Source port $(port) not found."))
        end
    end

    return nothing
end
