
# The result types of the harmonic balance solves, and the entry point
# `hbsolve` which runs the nonlinear pump solve and the linearized signal
# sweep in sequence.

"""
    NonlinearHB(w, frequencies, nodeflux, Rbnm, Ljb, Lb, Ljbm, Nmodes,
        Nbranches, nodes, ports, modes, S, solverinfo, operatingpoint,
        dcnodevoltage)

The solution of the nonlinear harmonic balance problem returned by
[`hbnlsolve`](@ref).

# Fields
- `w`: the tuple of pump angular frequencies in radians per second, one
    per non-commensurate pump.
- `frequencies`: the [`Frequencies`](@ref) describing the retained pump
    harmonics and intermodulation products.
- `nodeflux`: the node fluxes at each retained mode. With
    `keyedarrays = true` an `Nmodes` by `Nnodes - 1` keyed array with axes
    `:outputmode` and `:node`; otherwise a vector of length
    `Nmodes*(Nnodes - 1)` with the mode index varying fastest.
- `Rbnm`: the incidence matrix between the branch and node bases, with
    each entry repeated `Nmodes` times.
- `Ljb`: sparse vector of the Josephson junction inductances by branch.
- `Lb`: sparse vector of the linear inductances by branch.
- `Ljbm`: `Ljb` with each entry repeated `Nmodes` times.
- `Nmodes`: the number of retained modes.
- `Nbranches`: the number of branches in the circuit graph.
- `nodes`: the node names, ground (`"0"`) first.
- `ports`: the port numbers.
- `modes`: the retained modes as tuples of harmonic indices; `(1,)` is
    the pump of a single pump solve, `(1,0)` the first of two pumps.
- `S`: the scattering matrix at the pump frequencies, relating the inputs
    and outputs at each combination of port and mode. Its zero frequency
    entries are identically zero: the waves are in units of
    `sqrt(photons/second)`, whose normalization `1/sqrt(|w|)` (see
    [`portwavescale`](@ref)) has no limit at zero, so there is no direct
    current wave to report. The direct current operating point is in
    `dcnodevoltage`.
- `solverinfo`: diagnostics of the solution process; see
    [`SolverInfo`](@ref).
- `operatingpoint`: the converged operating point with the exact real
    Jacobian there, when the solver was called with
    `returnoperatingpoint = true`; otherwise `nothing`. See
    [`HBOperatingPoint`](@ref).
- `dcnodevoltage`: the average voltage of each non-ground node in volts,
    keyed by node name like `nodeflux`, when the analysis has a zero
    frequency mode; otherwise `nothing`. This is distinct from the zero
    mode of `nodeflux`, which is the static flux setting the inductor
    currents and junction phases. A vector of zeros is an answer (a node
    shorted to ground sits at zero volts, and so does every node of a
    circuit into which no direct current is injected) and differs from
    `nothing`.
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
    dcnodevoltage
end

function NonlinearHB(w, frequencies, nodeflux, Rbnm, Ljb, Lb, Ljbm, Nmodes,
    Nbranches, nodes, ports, modes, S, solverinfo, operatingpoint)
    return NonlinearHB(w, frequencies, nodeflux, Rbnm, Ljb, Lb, Ljbm,
        Nmodes, Nbranches, nodes, ports, modes, S, solverinfo,
        operatingpoint, nothing)
end

# Constructors with the trailing fields omitted: without the direct current
# voltages (above), and without the operating point and the diagnostics
# (below), which default to `nothing` and to an empty, converged record.
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
    LinearizedHB(w, modes, S, Snoise, Cnoise, Ssensitivity, QE, QEideal,
        CM, nodeflux, nodefluxadjoint, voltage, voltageadjoint, nodenames,
        nodeindices, componentnames, componenttypes, componentnamedict,
        mutualinductorbranchnames, portnumbers, portindices,
        portimpedances, noiseportimpedanceindices, sensitivitynames,
        sensitivityindices, Nmodes, Nnodes, Nbranches, Nports, signalindex)

The solution of the linearized harmonic balance problem returned by
[`hblinsolve`](@ref). An output which was not requested is an empty
array. The frequency dependent outputs are indexed as
`[outputmode, outputport, inputmode, inputport, frequency]` (keyed arrays
with those axis names when `keyedarrays = true`).

# Fields
- `w`: the signal angular frequencies in radians per second.
- `modes`: the retained signal modes as tuples of harmonic indices; `(0,)`
    is the signal itself and `(k,)` the idler offset by `k` pump harmonics.
- `S`: the scattering matrix relating the inputs and outputs at each
    combination of port, mode and signal frequency.
- `Snoise`: the scattering matrix from the noise channels of the
    dissipative elements to the ports. The channels are the dissipative
    lumped components (see [`noiseindices`](@ref)) followed by one per port
    of each dissipative [`ScatteringParameters`](@ref). Being a scattering
    matrix it describes a transformation and does not depend on
    temperature.
- `Cnoise`: the added noise covariance at the output ports,
    `sum_c occupation[c]*Snoise[c,i]*conj(Snoise[c,j])`, when
    `returnCnoise = true`. Together with `S` this is the Gaussian channel
    the circuit implements, taking an input covariance to
    `S*sigma*S' + Cnoise`; the temperature of each channel enters here
    through the occupation. At zero temperature its diagonal is the noise
    term in the denominator of the quantum efficiency.
- `Ssensitivity`: the derivative of `S` with respect to a relative
    (logarithmic) perturbation of each component in `sensitivitynames`, or
    of each design parameter when the sensitivity pair interface is used,
    at a fixed pump operating point or including the shift of the operating
    point when `sensitivityoperatingpoint = true`.
- `QE`: the quantum efficiency at each combination of port, mode and
    frequency.
- `QEideal`: the quantum efficiency of an ideal amplifier with the same
    gain.
- `CM`: the bosonic commutation relations of each output, `sum_j
    |S[i,j]|^2 sign(w_j)` over the input modes, which equal `+1` for an
    output at positive frequency and `-1` for one at negative frequency
    when the scattering matrix is complete.
- `nodeflux`: the node fluxes resulting from a unit input at each port and
    mode.
- `nodefluxadjoint`: the node fluxes of the adjoint (time reversed
    modulation) problem.
- `voltage`: the node voltages resulting from a unit input at each port
    and mode.
- `voltageadjoint`: the node voltages of the adjoint problem.
- `nodenames`: the node names, ground first.
- `nodeindices`: the 2 by `Ncomponents` matrix of component node indices
    from the [`CompiledCircuit`](@ref).
- `componentnames`, `componenttypes`, `componentnamedict`,
    `mutualinductorbranchnames`: the corresponding fields of the
    [`CompiledCircuit`](@ref).
- `portnumbers`: the port numbers, in the order the port axes use.
- `portindices`: the flat component index of each port.
- `portimpedances`: the reference impedance of each port, which the
    scattering parameters are normalized to.
- `noiseportimpedanceindices`: the flat component indices of the internal
    dissipative components, in the order of the noise channel axis of
    `Snoise`.
- `sensitivitynames`: the component names the sensitivities were taken
    with respect to.
- `sensitivityindices`: their flat component indices.
- `Nmodes`: the number of retained signal modes.
- `Nnodes`: the number of nodes, including ground.
- `Nbranches`: the number of branches in the circuit graph.
- `Nports`: the number of ports.
- `signalindex`: the position of the signal mode `(0,)` in `modes`, which
    is always 1.
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

# A constructor without `Cnoise`, which is an empty array when it was not
# requested.
LinearizedHB(w, modes, S, Snoise, Ssensitivity, QE, QEideal, CM, nodeflux, nodefluxadjoint, voltage, voltageadjoint, nodenames, nodeindices, componentnames, componenttypes, componentnamedict, mutualinductorbranchnames, portnumbers, portindices, portimpedances, noiseportimpedanceindices, sensitivitynames, sensitivityindices, Nmodes, Nnodes, Nbranches, Nports, signalindex) =
    LinearizedHB(w, modes, S, Snoise, Array{Complex{Float64},3}(undef, 0, 0, 0), Ssensitivity, QE, QEideal, CM, nodeflux, nodefluxadjoint, voltage, voltageadjoint, nodenames, nodeindices, componentnames, componenttypes, componentnamedict, mutualinductorbranchnames, portnumbers, portindices, portimpedances, noiseportimpedanceindices, sensitivitynames, sensitivityindices, Nmodes, Nnodes, Nbranches, Nports, signalindex)

"""
    HB(nonlinear, linearized)

The result of [`hbsolve`](@ref): the [`NonlinearHB`](@ref) solution for
the pump and its harmonics in `nonlinear`, and the [`LinearizedHB`](@ref)
solution for the signals in `linearized`.
"""
struct HB
    nonlinear
    linearized
end

# Docstring fragments shared by `hbsolve`, `hbnlsolve` and `hblinsolve`, so
# that one keyword is described in one place.
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

const _DOC_NLKWARGS = """
- `rtol = 0.0`: a relative residual tolerance; the solve is converged when
    `norm(F) <= max(ftol, rtol*norm(F0))` with `F0` the initial residual.
- `x0 = nothing`: an initial value for the node fluxes, either of the node
    flux length or of the full augmented length including the auxiliary
    variables of the modified nodal analysis formulation.
- `keyedarrays = true`: return `nodeflux` and `S` as keyed arrays with named
    axes rather than plain arrays.
- `sensitivitynames::Vector{String} = String[]`: the components whose
    indices are recorded for the sensitivity calculation.
- `returnoperatingpoint = false`: assemble and return the exact real
    Jacobian at the converged solution in the `operatingpoint` field, for
    sensitivities which include the shift of the operating point.
- `krylovcouplingmodes = :none`: the mode couplings the Newton-Krylov
    preconditioner keeps in full: `:none` (the mode block diagonal, the
    default), `:all`, `:band => p` for couplings up to the harmonic offset
    `p` of each tone (grown automatically on repeated linear failures),
    `:auto` or `:auto => tol` for a band whose width is measured from the
    solution, `:clusters` for clusters of modes measured from the coupling
    strengths of the operator (which, with a
    [`BlockFactorization`](@ref), is the form for three or more tones), a
    vector of mode indices, or an `Nmodes` by `Nmodes` `Bool` mask. See
    [`ModeCouplingPreconditioner`](@ref). For two strong tones
    `:auto` is the setting that matters: the block diagonal's error is
    then of high rank (every mode couples to every other through the
    pump mixing), and on 128- and 256-junction lines the measured band
    converges in 41 and 108 Arnoldi steps against about 3000 for the
    block diagonal, three to five times faster on a GPU, with no
    escalation to the full Jacobian.
- `krylovrecycle = 0`: when positive, wrap the preconditioner in a
    [`RecyclingPreconditioner`](@ref) keeping a deflation subspace of at
    most this many vectors across Newton steps and, through the reuse
    object of a cached sweep, across its points; `krylovharvest = 8` is the
    number of vectors harvested from each GMRES solve. This pays when the
    block diagonal leaves a deficiency of rank comparable to `krylovrecycle`
    (a two-tone line of about a hundred junctions converges without any
    factorization of the full Jacobian, in half the time) and costs when
    the deficiency is of high rank and the escalation does the work anyway,
    or when the base is strong enough on its own (`:auto`). Recycling is a
    correction layer for the few global directions a structural base
    leaves, not a substitute for the base; the directions it finds are
    the small singular directions of the preconditioned operator, the
    channels of large linear response, rather than its eigenmodes.
- `krylovdeflationform = :adef1`: how the recycled subspace is applied,
    `:adef1` (the projection on the input of the base solve) or `:adef2`
    (on its output, fused with the Jacobian product of the Arnoldi step);
    see [`RecyclingPreconditioner`](@ref). `:floquet` selects the
    residual-image form with physical candidates instead
    ([`FloquetPreconditioner`](@ref)), which harvests harmonic Ritz
    directions alongside the singular ones, from every restart cycle, and
    drops the directions the base already handles.
- `krylovdeflationkwargs::NamedTuple = (;)`: further keyword arguments of
    the deflation wrapper's constructor ([`RecyclingPreconditioner`](@ref)
    or [`FloquetPreconditioner`](@ref)), such as `nritz`, `kcandidate`,
    `benefittol` or `escalateafter`, which override its defaults.
- `krylovkwargs::NamedTuple = (;)`: further options of
    [`nlsolvekrylov!`](@ref), which override its defaults.
- `linearsolver = InternalGMRES()`: the linear solver of the Newton-Krylov
    step, [`InternalGMRES`](@ref) or [`KrylovJL`](@ref).
- `factorization = nothing`: the sparse factorization used by the direct
    methods and the preconditioner: [`KLUfactorization`](@ref) on the host
    and [`CUDSSFactorization`](@ref) on a CUDA device unless given.
- `backend = CPU()`: the KernelAbstractions backend the solve runs on.
- `precision = Float64`: the floating point type of the solve.
- `debugJacobian = false`: instead of solving, return a named tuple with
    the residual and Jacobian functions and the ingredients they are
    assembled from, for building reference implementations in tests.
- `returnsystem = false`: instead of solving, return a named tuple with
    the [`HBSystem`](@ref), the initial real state and residual, the real
    representation layout and (when `assemblejacobian = true`) the
    assembled real Jacobian, for driving an external solver.
- `switchofflinesearchtol`, `alphamin`: deprecated and ignored with a
    warning."""

const _DOC_RETURNS = """
- `returnS = true`: return the scattering parameters of the linearized
    solve.
- `returnSnoise = false`: return the noise scattering parameters.
- `returnCnoise = false`: return the added noise covariance at the output
    ports; see [`LinearizedHB`](@ref).
- `returnQE = true`: return the quantum efficiency.
- `returnCM = true`: return the commutation relations.
- `returnnodeflux = false`, `returnvoltage = false`: return the node fluxes
    and voltages of the linearized solve.
- `returnnodefluxadjoint = false`, `returnvoltageadjoint = false`: return
    the node fluxes and voltages of the adjoint (time reversed modulation)
    linearized solve.
- `keyedarrays = true`: return the outputs as keyed arrays with named,
    labeled axes rather than plain arrays."""

const _DOC_TEMPERATURE = """
- `temperature = 0.0`: the physical temperature in Kelvin of every
    dissipative element which does not state its own, and so of the noise
    it adds. A channel at temperature `T` carries `coth(hbar*w/(2*k*T))`
    times its vacuum noise, which at zero temperature is the vacuum noise
    itself. Raising it lowers the quantum efficiency and changes `Cnoise`
    but leaves `Snoise` and the commutation relations alone, since those
    describe the transformation rather than the state. The ports are
    vacuum by definition. A component may state its own temperature in
    the typed format (`Resistor(R; temperature = T)`, or a
    [`ScatteringParameters`](@ref) with `noise = ThermalEquilibrium(T)`);
    a tuple netlist cannot, and takes this default throughout."""

const _DOC_SENSNAMES = """
- `sensitivitynames::Vector{String} = String[]`: the names of the
    components to take sensitivities with respect to. Supported types are
    `C`, `L`, `R` and `Lj` with numeric values. A
    [`ScatteringParameters`](@ref) block has no scalar value to perturb and
    cannot be named; see [`designsensitivities`](@ref) for sensitivities
    with respect to the parameters of a block.
- `sensitivitypairs`, `sensitivityblockpairs`, `nsensitivityparameters`,
    `sensitivitylabels`: the design parameter interface used by
    [`designsensitivities`](@ref), which names physical parameters rather
    than components; each pair `(componentname, parameterindex, alpha)`
    gives the relative direction `alpha = (dv/dp)/v` of a component value
    under a parameter. Not intended to be passed directly."""

const _DOC_SENSMODE = """
- `sensitivitymode = :auto`: the order in which the operating point shift
    is contracted into the sensitivities. `:forward` costs one product
    against the linearized system per component and signal frequency;
    `:reverse` pushes the output functionals through the transposed pump
    Jacobian once per output port and mode pair, so its cost does not grow
    with the number of components. `:auto` chooses `:reverse` when there
    are more components than output port and mode pairs. Both support any
    number of pumps."""

const _DOC_SSENS = """
- `returnSsensitivity = false`: return `dS/dr`, the derivative of the
    scattering matrix with respect to a relative perturbation `r` of each
    named component value (`p -> r*p` at `r = 1`), computed by the adjoint
    method."""

const _DOC_NBATCHES = """
- `nbatches = Base.Threads.nthreads()`: the number of batches the signal
    frequencies are split into for multithreading; `1` runs single
    threaded."""

const _DOC_LINBACKEND = """
- `backend = CPU()`: the KernelAbstractions backend the sweep is solved
    on. On a device the system matrices of a batch of signal frequencies,
    which share one sparsity pattern, are assembled by one kernel and
    factorized and solved as a uniform batch (see
    [`CUDSSFactorization`](@ref)). The sweep falls back to the host when
    the component values depend on the symbolic frequency variable, or when
    an output needs the adjoint (transposed) solve: the noise scattering
    parameters, the scattering parameter sensitivities, or the adjoint node
    outputs."""

const _DOC_SORTING = """
- `sorting = :number`: how the nodes are ordered, with ground always first.
    `:number` parses the node names as integers and sorts numerically
    (an error if a name is not an integer); `:name` sorts the names as
    strings, so that "101" comes before "11"; `:none` keeps the order of
    first appearance. The methods taking a typed [`Circuit`](@ref) default
    to `:name`, since hierarchical net names are not integers."""

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
- `sourcefold`: always `NaN`. Reserved for a source stepping continuation
    which is not implemented; kept so the field layout does not change.
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
        circuit, circuitdefs; dc = false, threewavemixing = false,
        fourwavemixing = true, maxpumpintermodorder = Inf,
        maxmodulationintermodorder = Inf,
        Nevaluationharmonics = map(i -> 2i, Npumpharmonics),
        frequencywindow = (0, Inf),
        maxmodulationharmonics = Nmodulationharmonics,
        iterations = 1000, ftol = 1e-8, method = :newtonkrylov,
        andersondepth = method == :quasinewton ? 5 : 0, x0 = nothing,
        symfreqvar = nothing, nbatches = Base.Threads.nthreads(),
        sorting = :number, returnS = true, returnSnoise = false,
        returnQE = true, returnCM = true, returnnodeflux = false,
        returnvoltage = false, returnnodefluxadjoint = false,
        returnvoltageadjoint = false, keyedarrays = true,
        temperature = 0.0, returnCnoise = false,
        sensitivitynames::Vector{String} = String[],
        sensitivityoperatingpoint = true, sensitivitymode = :auto,
        returnSsensitivity = false, krylovcouplingmodes = :none,
        krylovrecycle = 0, krylovharvest = 8, krylovdeflationform = :adef1,
        krylovdeflationkwargs = (;), krylovkwargs = (;), stagedkwargs = (;),
        factorization = nothing,
        backend = CPU(), precision = Float64)

Solve a circuit driven by one or more strong pumps, then linearize about
that operating point and sweep the weak signal frequencies `ws`. This
calls [`hbnlsolve`](@ref) for the pump and [`hblinsolve`](@ref) for the
signals and returns both solutions in an [`HB`](@ref). Any number of
pumps, ports and modes is supported, for three and four wave mixing.

The system is solved in a modified nodal analysis formulation in the node
flux basis. Resistors with constant real values and mutually coupled
inductor branches are given auxiliary branch current variables, which
keeps the system matrix bounded as a coupling coefficient approaches one,
and one gauge fixing equation per floating inductive or Josephson
subnetwork and zero frequency mode makes circuits with no inductive path
to ground solvable without workaround inductors (see `src/mna.jl`). The
nonlinear system is nondimensionalized by the scale `Z0/w0` (see
[`calcsolverscale`](@ref)), so `ftol` does not depend on the unit system.
The returned node fluxes and voltages contain only the node coordinates,
not the auxiliary variables. The linearized solve throws an
`ArgumentError` when any signal plus pump mode frequency is numerically
zero; estimate a direct current limit from a sequence of decreasing
nonzero frequencies instead.

# Arguments
- `ws`: the signal angular frequency or frequencies in radians per second,
    such as `2*pi*5.0e9` or `2*pi*(4.5:0.001:5.0)*1e9`.
- `wp::NTuple{N,Number}`: the pump angular frequencies in radians per
    second, `(2*pi*5.0e9,)` for a single pump or `(2*pi*5.0e9, 2*pi*6.0e9)`
    for two. The pumps should be non-commensurate; for commensurate pumps
    give the lowest frequency here and add the others to `sources` with a
    mode index equal to the frequency ratio.
- `sources::Vector`: a vector of named tuples `(mode, port, current)`.
    `mode` is a tuple of harmonic indices, one per pump frequency, `port`
    the port number, and `current` the complex current amplitude in
    Amperes. `[(mode=(1,0), port=1, current=Ip1), (mode=(0,1), port=1,
    current=Ip2)]` applies a source at `1*wp[1] + 0*wp[2]` and one at
    `0*wp[1] + 1*wp[2]`, both at port 1. A source with `mode = (0,)` and
    `dc = true` is a direct current bias.
- `Nmodulationharmonics::NTuple{M,Int}`: how many harmonics of each pump to
    retain around the signal in the linearized solve, which sets the
    signal and idler modes.
- `Npumpharmonics::NTuple{N,Int}`: how many harmonics of each pump to
    retain as unknowns in the nonlinear solve. Its length is the number of
    non-commensurate pumps. The nonlinearity is evaluated on the larger
    `Nevaluationharmonics` grid.
- `circuit`: a typed [`Circuit`](@ref), a legacy netlist of
    `(name, node1, node2, value)` tuples, or a [`CompiledCircuit`](@ref).
- `circuitdefs`: a dictionary from the symbols or symbolic variables used
    as component values to their numerical values. Optional when every
    component value is numeric.

# Keywords
- `dc = false`: retain the zero frequency mode in the nonlinear solve.
- `threewavemixing = false`: retain the even pump harmonics, which are
    what three wave mixing processes couple through.
- `fourwavemixing = true`: retain the odd pump harmonics.
- `maxpumpintermodorder = Inf`: keep only the pump modes whose harmonic
    indices have an absolute sum of at most this order, a diamond
    truncation of the multi-pump Fourier space.
- `maxmodulationintermodorder = Inf`: the same truncation for the signal
    modes.
- `Nevaluationharmonics = map(i -> 2i, Npumpharmonics)`: the
    harmonics of each pump on the grid where the nonlinearity is sampled,
    at least `Npumpharmonics`; twice the retained set by default, which
    dealiases the products of order three, the leading ones of a junction.
    See [`hbnlsolve`](@ref).
- `maxpumpharmonics`: deprecated and ignored with a warning;
    `Npumpharmonics` is the retained set and `Nevaluationharmonics` the
    sampling grid.
- `maxmodulationharmonics = Nmodulationharmonics`: an upper bound on the
    absolute harmonic index retained for each pump in the signal modes,
    applied together with the intermodulation truncation; see
    [`truncfreqs`](@ref).
- `frequencywindow = (0, Inf)`: a lower and upper bound, in the units of
    `wp`, on the absolute frequency of the retained pump modes; see
    [`hbnlsolve`](@ref). The signal modes are not windowed.
- `iterations = 1000`: the maximum number of nonlinear solver iterations
    before it returns unconverged.
$(_DOC_FTOL)
$(_DOC_METHOD)
$(_DOC_STAGEDKWARGS)
- `krylovcouplingmodes = :none`: which mode couplings the Newton-Krylov
    preconditioner retains; forwarded to [`hbnlsolve`](@ref).
- `krylovrecycle = 0`, `krylovharvest = 8`, `krylovdeflationform = :adef1`:
    the recycled deflation subspace of the Newton-Krylov preconditioner;
    forwarded to [`hbnlsolve`](@ref).
- `krylovkwargs::NamedTuple = (;)`, `krylovdeflationkwargs::NamedTuple = (;)`:
    further options of the Newton-Krylov solver and of its deflation
    wrapper; forwarded to [`hbnlsolve`](@ref).
$(_DOC_ANDERSON)
- `x0 = nothing`: an initial value for the node fluxes of the nonlinear
    solve.
- `symfreqvar = nothing`: the symbolic frequency variable, such as `w`,
    when component values are expressions in the frequency.
$(_DOC_NBATCHES)
$(_DOC_SORTING)
$(_DOC_RETURNS)
$(_DOC_TEMPERATURE)
$(_DOC_SENSNAMES)
$(_DOC_SENSMODE)
- `sensitivityoperatingpoint = true`: include the shift of the pump
    operating point in the sensitivities, making them total derivatives
    rather than derivatives at a fixed operating point. Near the gain peak
    of a strongly pumped amplifier the operating point term is comparable
    to or larger than the fixed point term. Requires the exact real
    Jacobian of the nonlinear solution, which is assembled and factorized
    once. Without Josephson junctions the operating point contribution is
    identically zero and is skipped.
$(_DOC_SSENS)
- `factorization = nothing`: the sparse factorization used by the direct
    and the linearized solves. `nothing` selects [`KLUfactorization`](@ref)
    on the host and [`CUDSSFactorization`](@ref) on a CUDA device for the
    nonlinear solve, and `KLUfactorization` for the linearized solve.
- `backend = CPU()`: the KernelAbstractions backend both solves run on. The
    nonlinear solve assembles, factorizes and iterates there; the
    linearized sweep solves batches of signal frequencies there and falls
    back to the host for what it cannot serve (see [`hblinsolve`](@ref)).
- `precision = Float64`: the floating point type of the nonlinear solve.
- `switchofflinesearchtol`, `alphamin`: deprecated and ignored with a
    warning.
- `returnZ`, `returnZadjoint`, `returnZsensitivity`,
    `returnZsensitivityadjoint`: removed; passing any of them warns.
    Compute impedances from the scattering parameters instead.

# Returns
- `HB`: the nonlinear and linearized solutions; see [`HB`](@ref).

"""
function hbsolve(ws, wp::NTuple{N,Number}, sources::Vector,
    Nmodulationharmonics::NTuple{M,Int}, Npumpharmonics::NTuple{N,Int},
    circuit, circuitdefs;dc::Bool = false, threewavemixing::Bool = false,
    fourwavemixing::Bool = true, maxpumpintermodorder=Inf,
    maxmodulationintermodorder=Inf,
    Nevaluationharmonics::NTuple{N,Int} = map(i -> 2i, Npumpharmonics),
    maxpumpharmonics = nothing,
    frequencywindow = (0, Inf),
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
    krylovcouplingmodes = :none, krylovrecycle::Integer = 0,
    krylovharvest::Integer = 8, krylovdeflationform::Symbol = :adef1,
    krylovkwargs::NamedTuple = (;),
    krylovdeflationkwargs::NamedTuple = (;),
    stagedkwargs::NamedTuple = (;),
    factorization = nothing, backend = CPU(),
    precision::Type{<:AbstractFloat} = Float64) where {N,M}

    # deprecation warning for maxpumpharmonics, whose role `Npumpharmonics`
    # took when the sampling grid became `Nevaluationharmonics`.
    if !isnothing(maxpumpharmonics)
        Base.depwarn(lazy"The `maxpumpharmonics` kwarg is deprecated and no longer used. `Npumpharmonics` is the retained set of pump modes and `Nevaluationharmonics` the grid on which the nonlinearity is sampled. Please remove it to avoid errors in future versions.", :hbsolve; force=true)
    end

    all(map(>=, Nevaluationharmonics, Npumpharmonics)) || throw(ArgumentError(
        lazy"`Nevaluationharmonics` = $(Nevaluationharmonics) must be at least `Npumpharmonics` = $(Npumpharmonics) in every pump."))

    # the pump modes: harmonics of each pump and their intermodulation
    # products, truncated, with the conjugate (negative frequency) modes
    # removed since the real signal determines them; the nonlinearity is
    # sampled on the larger `Nevaluationharmonics` grid
    freq = removeconjfreqs(
        truncfreqs(
            calcfreqsrdft(Nevaluationharmonics); dc = dc,
            odd = fourwavemixing, even = threewavemixing,
            maxintermodorder = maxpumpintermodorder,
            maxharmonics = Npumpharmonics, w = wp,
            frequencywindow = frequencywindow,
        )
    )

    indices = fourierindices(freq)

    Nmodes = length(freq.modes)

    # compile the circuit, build its graph, and assemble the matrices at
    # the pump mode count
    psc, cg, nm = preparecircuit(circuit, circuitdefs;
        sorting = sorting, Nmodes = Nmodes)


    # The nonlinear solve. `:staged` runs the source continuation driver,
    # which rebuilds its own truncations from the same arguments and
    # manages its own warm starts; every other method solves the system
    # assembled above.
    nonlinear = if method === :staged
        stagedhbnlsolve(wp, Npumpharmonics, sources, circuit, circuitdefs;
            iterations = iterations, ftol = ftol,
            Nevaluationharmonics = Nevaluationharmonics,
            maxintermodorder = maxpumpintermodorder,
            dc = dc, odd = fourwavemixing, even = threewavemixing,
            symfreqvar = symfreqvar, sorting = sorting,
            keyedarrays = keyedarrays, sensitivitynames = sensitivitynames,
            returnoperatingpoint = sensitivityoperatingpoint &&
                returnSsensitivity && !isempty(nm.Ljb.nzind),
            krylovcouplingmodes = krylovcouplingmodes,
            krylovrecycle = krylovrecycle, krylovharvest = krylovharvest,
            krylovdeflationform = krylovdeflationform,
            krylovdeflationkwargs = krylovdeflationkwargs,
            krylovkwargs = krylovkwargs, factorization = factorization,
            backend = backend, precision = precision, stagedkwargs...)
    else
        hbnlsolve(wp, sources, freq, indices, psc, cg, nm;
            iterations = iterations, x0 = x0, ftol = ftol,
            switchofflinesearchtol = switchofflinesearchtol,
            alphamin = alphamin,
            method = method, andersondepth = andersondepth,
            krylovcouplingmodes = krylovcouplingmodes,
            krylovrecycle = krylovrecycle, krylovharvest = krylovharvest,
            krylovdeflationform = krylovdeflationform,
            krylovdeflationkwargs = krylovdeflationkwargs,
            krylovkwargs = krylovkwargs,
                symfreqvar = symfreqvar, keyedarrays = keyedarrays,
            sensitivitynames = sensitivitynames,
            returnoperatingpoint = sensitivityoperatingpoint &&
                returnSsensitivity && !isempty(nm.Ljb.nzind),
            factorization = factorization, backend = backend,
            precision = precision)
    end

    # The derivative of the harmonic balance residual with respect to each
    # sensitivity component, needed for sensitivities which include the
    # shift of the operating point. The matrices `nm` were built at the
    # pump mode count, which is what `calcresidualsensitivity` expects. The
    # operating point derivatives dx/dr are computed in `hblinsolve` only
    # when it selects the forward contraction order. Without Josephson
    # junctions the linearized system does not depend on the operating
    # point at all, so the contribution is zero and is skipped.
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
            # one column per (component, parameter) pair with the direction
            # alpha folded in; hblinsolve merges them per parameter
            calcresidualsensitivity(nonlinear.operatingpoint, psc, cg, nm,
                [psc.componentnamedict[String(t[1])]
                 for t in sensitivitypairs],
                [Complex{Float64}(t[3]) for t in sensitivitypairs])
        end
        # the scattering block pairs add their columns after the lumped
        # ones, in the order hblinsolve expects
        if !isempty(sensitivityblockpairs)
            blockcols = calcblockresidualsensitivity(
                nonlinear.operatingpoint, psc, sensitivityblockpairs)
            cols = isnothing(cols) ? blockcols : hcat(cols, blockcols)
        end
        cols
    else
        nothing
    end

    # the signal modes: the signal and the idlers offset from it by pump
    # harmonics
    signalfreq =truncfreqs(
        calcfreqsdft(Nmodulationharmonics); dc = true, odd = threewavemixing,
        even = fourwavemixing, maxintermodorder = maxmodulationintermodorder,
        maxharmonics = maxmodulationharmonics,
    )

    # The linearized solve. The component values were resolved for the
    # nonlinear solve, so the signal side matrices are built from `nm.vvn`
    # rather than resolving the values again at the signal mode count.
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
        # the host factorization is the default here even on a device
        # backend, because the fallback path of the sweep is a host solve
        # of complex matrices; an explicit factorization is honored
        factorization = isnothing(factorization) ? KLUfactorization() :
            factorization,
        backend = backend)

    return HB(nonlinear, linearized)
end

# A fully numeric circuit needs no component definitions.
function hbsolve(ws, wp::NTuple{N,Number}, sources::Vector,
    Nmodulationharmonics::NTuple{M,Int}, Npumpharmonics::NTuple{N,Int},
    circuit; kwargs...) where {N,M}
    return hbsolve(ws, wp, sources, Nmodulationharmonics, Npumpharmonics,
        circuit, Dict{Any,Any}(); kwargs...)
end

"""
    hbsolve(ws, wp, sources, Nmodulationharmonics, Npumpharmonics,
        circuit::Circuit, circuitdefs = Dict{Symbol,Number}();
        sorting = :name, kwargs...)

Solve a typed [`Circuit`](@ref). The circuit is elaborated and lowered
with [`compile`](@ref) and every keyword of the general method applies.
`circuitdefs` is needed only when component values are symbolic. The
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
