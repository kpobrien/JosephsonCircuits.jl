
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
end

# backwards compatible constructor without the solver diagnostics
function NonlinearHB(w, frequencies, nodeflux, Rbnm, Ljb, Lb, Ljbm, Nmodes,
    Nbranches, nodes, ports, modes, S)
    return NonlinearHB(w, frequencies, nodeflux, Rbnm, Ljb, Lb, Ljbm,
        Nmodes, Nbranches, nodes, ports, modes, S,
        SolverInfo(IterationInfo[], NaN, NaN, true, NaN))
end

"""
    LinearizedHB(S, Snoise, QE, QEideal, CM, nodeflux, voltage, Nmodes, Nnodes,
        Nbranches, signalindex, w)

A simple structure to hold the linearized harmonic balance solutions.

# Fields
- `w`: the signal frequencies.
- `modes`: tuple of the signal mode indices where (0,) is the signal.
- `S`: the scattering matrix relating inputs and outputs for each combination
    of port and frequency.
- `Snoise`: the scattering matrix relating inputs at the noise ports.
    (lossy devices) and outputs at the physical ports for each combination of
    port and frequency.
- `Ssensitivity`:
- `Z`:
- `Zadjoint`: 
- `Zsensitivity`: 
- `Zsensitivityadjoint`: 
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
    Ssensitivity
    Z
    Zadjoint
    Zsensitivity
    Zsensitivityadjoint
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
- `stages`: a vector of [`IterationInfo`](@ref), one for each invocation
    of the nonlinear solver, in the order they ran. The labels record
    which solver stages ran and, for `method = :sourcestepping`, the
    parameters record the source scale of each continuation point.
- `initialresidual`: the norm of the residual at the initial value.
- `finalresidual`: the norm of the residual at the returned solution.
- `converged`: whether the solver reported convergence.
- `sourcefold`: the source scale at which the source stepping detected a
    fold in the branch of solutions, or NaN if no fold was detected. A
    fold is recorded even when the subsequent extrapolated jump over it
    succeeds; when the solve also failed, a finite value suggests that
    no solution exists on the branch at the full source amplitudes, for
    example because the sources exceed a self-oscillation threshold. The
    labels of `stages` record which continuation stage stalled.
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
        returnSsensitivity = false, returnZ = false, returnZadjoint = false,
        returnZsensitivity = false, returnZsensitivityadjoint = false,
        factorization = KLUfactorization()) where {N,M,K}

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
- `ftol = 1e-8`: the function tolerance defined we considered converged,
    defined as norm(F)/norm(x) < ftol or norm(F,Inf) <= ftol.
- `method = :quasinewton`: the nonlinear solution method. `:quasinewton`
    iterates with the complex holomorphic Jacobian `Jx` only, an
    approximation to the exact Jacobian. `:newton` solves the equivalent real
    system with the exact real Jacobian, which restores quadratic convergence
    near the solution, including for multi-tone problems: its mode coupling
    indices include the couplings which alias back onto the sampled grid,
    matching the cyclic Fourier transforms of the residual, and the
    assembled Jacobian agrees with the matrix-free Jacobian-vector product
    of `HBSystem` to machine precision. Both Jacobians are assembled
    directly from the Fourier coefficients of `cos(phi(t))` with precomputed
    plans (see [`plancomplexjacobian`](@ref) and
    [`planrealjacobian`](@ref)).
- `andersondepth::Integer = method == :newton ? 0 : 5`: the depth of the
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
- `sensitivitynames::Vector{String} = String[]`: the component names for which
    to return the sensitivities (in progress).
- `returnSsensitivity = false`: return the scattering parameter sensitivity
    matrix from the linearized simulations (in progress).
- `returnZ = false`: return the impedance matrix from the linearized
    simulations.
- `returnZadjoint = false`: return the impedance matrix from the linearized
    adjoint simulations.
- `returnZsensitivity = false`: return the Z parameter sensitivity
    matrix from the linearized simulations (in progress).
- `returnZsensitivityadjoint = false`: return the Z parameter sensitivity
    matrix from the linearized adjoint simulations (in progress).
- `factorization = KLUfactorization()`: the type of factorization to use for
    the nonlinear and the linearized simulations.

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
    alphamin = nothing, method = :quasinewton,
    andersondepth::Integer = method == :newton ? 0 : 5,
    x0 = nothing, symfreqvar = nothing, nbatches = Base.Threads.nthreads(),
    sorting = :number, returnS::Bool = true, returnSnoise::Bool = false,
    returnQE::Bool = true, returnCM::Bool = true, returnnodeflux::Bool = false,
    returnvoltage::Bool = false, returnnodefluxadjoint::Bool = false,
    returnvoltageadjoint::Bool = false, keyedarrays::Bool = true,
    sensitivitynames::Vector{String} = String[],
    returnSsensitivity::Bool = false, returnZ::Bool = false,
    returnZadjoint::Bool = false, returnZsensitivity::Bool = false,
    returnZsensitivityadjoint::Bool = false,
    factorization = KLUfactorization()) where {N,M}

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

    # solve the nonlinear problem
    nonlinear = hbnlsolve(wp, sources, freq, indices, psc, cg, nm;
        iterations = iterations, x0 = x0, ftol = ftol,
        switchofflinesearchtol = switchofflinesearchtol, alphamin = alphamin,
        method = method, andersondepth = andersondepth,
        symfreqvar = symfreqvar, keyedarrays = keyedarrays,
        sensitivitynames = sensitivitynames, factorization = factorization)

    # generate the signal modes
    signalfreq =truncfreqs(
        calcfreqsdft(Nmodulationharmonics); dc = true, odd = threewavemixing,
        even = fourwavemixing, maxintermodorder = maxmodulationintermodorder,
        maxharmonics = maxmodulationharmonics,
    )

    # solve the linearized problem
    linearized = hblinsolve(ws, psc, cg, circuitdefs, signalfreq;
        nonlinear = nonlinear, symfreqvar = symfreqvar, nbatches = nbatches,
        returnS = returnS, returnSnoise = returnSnoise, returnQE = returnQE,
        returnCM = returnCM, returnnodeflux = returnnodeflux,
        returnnodefluxadjoint = returnnodefluxadjoint,
        returnvoltage = returnvoltage,
        returnvoltageadjoint = returnvoltageadjoint, 
        keyedarrays = keyedarrays, sensitivitynames = sensitivitynames,
        returnSsensitivity = returnSsensitivity,
        returnZ = returnZ, returnZadjoint = returnZadjoint,
        returnZsensitivity = returnZsensitivity,
        returnZsensitivityadjoint = returnZsensitivityadjoint,
        factorization = factorization)

    return HB(nonlinear, linearized)
end

"""
    hblinsolve(w, circuit,circuitdefs; Nmodulationharmonics = (0,),
        nonlinear=nothing, symfreqvar=nothing, threewavemixing=false,
        fourwavemixing=true, maxintermodorder=Inf,
        nbatches::Integer = Base.Threads.nthreads(), returnS = true,
        returnSnoise = false, returnQE = true, returnCM = true,
        returnnodeflux = false, returnnodefluxadjoint = false,
        returnvoltage = false,
        )

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
- `sensitivitynames::Vector{String} = String[]`: the component names for which
    to return the sensitivities (in progress).
- `returnSsensitivity = false`: return the scattering parameter sensitivity
    matrix from the linearized simulations (in progress).
- `returnZ = false`: return the impedance matrix from the linearized
    simulations.
- `returnZadjoint = false`: return the impedance matrix from the linearized
    adjoint simulations.
- `returnZsensitivity = false`: return the Z parameter sensitivity
    matrix from the linearized simulations (in progress).
- `returnZsensitivityadjoint = false`: return the Z parameter sensitivity
    matrix from the linearized adjoint simulations (in progress).
- `factorization = KLUfactorization()`: the type of factorization to use for
    the nonlinear and the linearized simulations.

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
    sensitivitynames::Vector{String} = String[],
    returnSsensitivity::Bool = false, returnZ::Bool = false,
    returnZadjoint::Bool = false, returnZsensitivity::Bool = false,
    returnZsensitivityadjoint::Bool = false,
    factorization = KLUfactorization())

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
        keyedarrays = keyedarrays, sensitivitynames = sensitivitynames,
        returnSsensitivity = returnSsensitivity,
        returnZ = returnZ, returnZadjoint = returnZadjoint,
        returnZsensitivity = returnZsensitivity,
        returnZsensitivityadjoint = returnZsensitivityadjoint,
        factorization = factorization)
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
    sensitivitynames::Vector{String} = String[],
    returnSsensitivity::Bool = false, returnZ::Bool = false,
    returnZadjoint::Bool = false, returnZsensitivity::Bool = false,
    returnZsensitivityadjoint::Bool = false,
    factorization = KLUfactorization())

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
        phimatrixtd, irfftplan, rfftplan = plan_applynl(phimatrix)

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
    wpumpmodes = if isnothing(nonlinear)
        calcmodefreqs((0.0,),signalfreq.modes)
    else
        calcmodefreqs(nonlinear.w,signalfreq.modes)
    end


    # the flux basis represents voltages as v = im*w*phi, which is
    # degenerate at exactly zero total frequency; scattering parameters at
    # DC are the w -> 0 limit. detect signal plus pump-mode totals which
    # are zero or cancel to within the accumulated roundoff of the
    # contributing terms.
    all(isfinite, w) || throw(ArgumentError("All signal frequencies must be finite."))
    let wpumptuple = isnothing(nonlinear) ? (0.0,) : nonlinear.w
        for wi in w
            for (mi, wm) in enumerate(wpumpmodes)
                mode = signalfreq.modes[mi]
                terms = vcat(float(real(wi)),
                    [float(real(mode[j]*wpumptuple[j])) for j in eachindex(wpumptuple)])
                if isnumericallyzero(wi + wm, terms)
                    throw(ArgumentError("hblinsolve cannot evaluate a mode at (numerically) zero total frequency (signal frequency plus pump mode frequency, here signal $(wi/(2*pi)) Hz with mode $(mode)) because the node flux basis represents voltages as v = im*w*phi. Zero-frequency small-signal analysis is not supported; to estimate a DC limit, evaluate a sequence of decreasing nonzero frequencies and verify that the requested network parameters converge. For frequency independent resistive networks the result at any nonzero frequency equals the DC limit."))
                end
            end
        end
    end

    # this is the first signal frequency. we will use it for various setup tasks
    wmodes = w[1] .+ wpumpmodes
    wmodesm = Diagonal(repeat(wmodes,outer=psc.Nnodes-1));
    wmodes2m = Diagonal(repeat(wmodes.^2,outer=psc.Nnodes-1));

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
    # current variables to the resistors with constant real values and to the
    # mutually coupled inductor branches (whose inverse inductance entries
    # would otherwise diverge as the coupling coefficient approaches one).
    # eliminating the auxiliary variables recovers the nodal equations. no
    # gauge fixing equations are needed because modes at (numerically) zero
    # total frequency are not permitted.
    checkstaticstiffnessvalues(psc.componenttypes, signalnm.vvn)
    mnaindices = mnaresistorindices(psc.componenttypes, signalnm.vvn)
    coupledbranches = mnacoupledbranches(signalnm.Mb)
    Nauxmnar = length(mnaindices)*Nsignalmodes
    Nauxmna = Nauxmnar + length(coupledbranches)*Nsignalmodes
    Nnodalmna = (psc.Nnodes-1)*Nsignalmodes
    Amna0, AmnaG = calcAmnasplit(mnaindices, psc.nodeindices, signalnm.vvn,
        Nsignalmodes, psc.Nnodes)
    Amna0 = mnapad(Amna0, length(coupledbranches)*Nsignalmodes)
    AmnaG = mnapad(AmnaG, length(coupledbranches)*Nsignalmodes)
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
        wmodesm = Diagonal(repeat(wmodes,
            outer = (Nnodalmna + Nauxmna) ÷ Nsignalmodes))
        wmodes2m = Diagonal(repeat(wmodes.^2,
            outer = (Nnodalmna + Nauxmna) ÷ Nsignalmodes))
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
    lsys = HBLinearizedSystem(Amatrixindices, signalnm.Ljb, Rbnmmna,
        Nsignalmodes, cg.Nbranches, phimatrix, invLnmcopy, Gnmcopy, Cnmcopy,
        invLnm, Gnm, Cnm, invLnmfreqsubstindices, Gnmfreqsubstindices,
        Cnmfreqsubstindices, Amna0, AmnaG, symfreqvar, wpumpmodes, Nnodes)
    Asparse = lsys.Asparse

    portimpedances = [vvn[i] for i in portimpedanceindices]
    noiseportimpedances = [vvn[i] for i in noiseportimpedanceindices]

    # assemble Asparse once at the first frequency so we have something
    # reasonable to factorize.
    assemblesystemmatrix!(Asparse, lsys, wmodesm, wmodes2m)


    # make arrays for the voltages, node fluxes, scattering parameters,
    # quantum efficiency, and commutatio relations. if we aren't returning a
    # matrix, set it to be an empty array.
    S = if returnS || returnQE 
        zeros(Complex{Float64}, Nports*Nsignalmodes, Nports*Nsignalmodes,
            length(w))
    else
        zeros(Complex{Float64},0,0,0)
    end

    Z = if returnZ
        zeros(Complex{Float64}, Nports*Nsignalmodes, Nports*Nsignalmodes,
            length(w))
    else
        zeros(Complex{Float64},0,0,0)
    end

    Zadjoint = if returnZadjoint
        zeros(Complex{Float64}, Nports*Nsignalmodes,
            Nports*Nsignalmodes, length(w))
    else
        zeros(Complex{Float64},0,0,0)
    end

    Snoise = if returnSnoise
        zeros(Complex{Float64},
            length(noiseportimpedanceindices)*Nsignalmodes,
            Nports*Nsignalmodes, length(w))
    else
        zeros(Complex{Float64},0,0,0)
    end

    Ssensitivity = if returnSsensitivity
        zeros(Complex{Float64},
            length(sensitivitynames)*Nsignalmodes,
            Nports*Nsignalmodes, length(w))
    else
        zeros(Complex{Float64},0,0,0)
    end

    Zsensitivity = if returnZsensitivity
        zeros(Complex{Float64},
            length(sensitivitynames)*Nsignalmodes,
            Nports*Nsignalmodes, length(w))
    else
        zeros(Complex{Float64},0,0,0)
    end

    Zsensitivityadjoint = if returnZsensitivityadjoint
        zeros(Complex{Float64},
            length(sensitivitynames)*Nsignalmodes,
            Nports*Nsignalmodes, length(w))
    else
        zeros(Complex{Float64},0,0,0)
    end

    QE = if returnQE
        zeros(Float64,Nports*Nsignalmodes,Nports*Nsignalmodes,length(w))
    else
        zeros(Float64,0,0,0)
    end

    CM = if returnCM
        zeros(Float64,Nports*Nsignalmodes,length(w))
    else
        zeros(Float64,0,0)
    end

    nodeflux = if returnnodeflux
        zeros(Complex{Float64},Nsignalmodes*(Nnodes-1),
            Nsignalmodes*Nports, length(w))
    else
        Vector{Complex{Float64}}(undef,0)
    end

    nodefluxadjoint = if returnnodefluxadjoint
        zeros(Complex{Float64}, Nsignalmodes*(Nnodes-1),
            Nsignalmodes*Nports, length(w))
    else
        Vector{Complex{Float64}}(undef,0)
    end

    voltage = if returnvoltage
        zeros(Complex{Float64}, Nsignalmodes*(Nnodes-1),
            Nsignalmodes*Nports, length(w))
    else
        Vector{Complex{Float64}}(undef,0)
    end

    voltageadjoint = if returnvoltageadjoint
        zeros(Complex{Float64}, Nsignalmodes*(Nnodes-1),
            Nsignalmodes*Nports, length(w))
    else
        Vector{Complex{Float64}}(undef,0)
    end

    # generate the mode indices and find the signal index
    signalindex = 1

    # solve the linear system for the specified frequencies. the response for
    # each frequency is independent so it can be done in parallel; however
    # we want to reuse the factorization object and other input arrays. 
    # perform array allocations and factorization "nbatches" times.
    # parallelize using tasks
    batches = Base.Iterators.partition(1:length(w),1+(length(w)-1)÷nbatches)
    Threads.@sync for batch in batches
        Base.Threads.@spawn hblinsolve_inner!(S, Snoise, Ssensitivity, Z, Zadjoint, Zsensitivity, Zsensitivityadjoint,
            QE, CM, nodeflux, nodefluxadjoint, voltage, voltageadjoint,
            lsys, bnm,
            portindices, portimpedanceindices, noiseportimpedanceindices, sensitivityindices,
            portimpedances, noiseportimpedances, nodeindices, componenttypes,
            w, wpumpmodes, Nsignalmodes, Nnodes, symfreqvar, batch, factorization)
    end

    # calculate the `ideal` quantum efficiency based on the gain assuming an
    # ideal two mode amplifier
    QEideal = if returnQE
        # calculate the ideal quantum efficiency.
        zeros(Float64,size(S))
    else
        zeros(Float64,0,0,0)
    end
    if returnQE
        # calculate the ideal quantum efficiency.
        calcqeideal!(QEideal,S)
    end

    # turn all of the array outputs into keyed arrays if the
    # keyedarrays = true

    # if keyword argument keyedarrays = true then generate keyed arrays
    # for the scattering parameters
    Sout = if returnS && keyedarrays
        Stokeyed(S, modes, portnumbers, modes, portnumbers, w)
    elseif returnS
        S
    else
	zeros(Complex{Float64},0,0,0)
    end

    Zout = if returnZ && keyedarrays
        Stokeyed(Z, modes, portnumbers, modes, portnumbers, w)
    else
        Z
    end

    Zadjointout = if returnZadjoint && keyedarrays
        Stokeyed(Zadjoint, modes, portnumbers, modes, portnumbers, w)
    else
        Zadjoint
    end

    # if keyword argument keyedarrays = true then generate keyed arrays
    # for the noise scattering parameters
    Snoiseout = if returnSnoise && keyedarrays
        Snoisetokeyed(Snoise, modes,
            componentnames[noiseportimpedanceindices], modes, portnumbers, w)
    else
        Snoise
    end

    # if keyword argument keyedarrays = true then generate keyed arrays
    # for Ssensitivity
    Ssensitivityout = if returnSsensitivity && keyedarrays
        Snoisetokeyed(Ssensitivity, modes,
            sensitivitynames, modes, portnumbers, w)
    else
        Ssensitivity
    end

    Zsensitivityout = if returnZsensitivity && keyedarrays
        Snoisetokeyed(Zsensitivity, modes,
            sensitivitynames, modes, portnumbers, w)
    else
        Zsensitivity
    end

    Zsensitivityadjointout = if returnZsensitivityadjoint && keyedarrays
        Snoisetokeyed(Zsensitivityadjoint, modes,
            sensitivitynames, modes, portnumbers, w)
    else
        Zsensitivityadjoint
    end

    # if keyword argument keyedarrays = true then generate keyed arrays
    # for the quantum efficiency
    QEout = if returnQE && keyedarrays
        Stokeyed(QE, modes, portnumbers, modes, portnumbers, w)
    else
        QE
    end

    QEidealout = if returnQE && keyedarrays
        Stokeyed(QEideal, modes, portnumbers, modes, portnumbers, w)
    else
        QEideal
    end

    # if keyword argument keyedarrays = true then generate keyed arrays
    # for the commutation relations
    CMout = if returnCM && keyedarrays
        CMtokeyed(CM, modes, portnumbers, w)
    else
        CM
    end

    # if keyword argument keyedarrays = true then generate keyed arrays
    nodefluxout = if returnnodeflux && keyedarrays
        nodevariabletokeyed(nodeflux, modes, nodenames, modes,
            portnumbers, w)
    else
        nodeflux
    end

    # if keyword argument keyedarrays = true then generate keyed arrays
    nodefluxadjointout = if returnnodefluxadjoint && keyedarrays
        nodevariabletokeyed(nodefluxadjoint, modes,
            nodenames, modes, portnumbers, w)
    else
        nodefluxadjoint
    end

    # if keyword argument keyedarrays = true then generate keyed arrays
    voltageout = if returnvoltage && keyedarrays
        nodevariabletokeyed(voltage, modes,
            nodenames, modes, portnumbers, w)
    else
        voltage
    end

    # if keyword argument keyedarrays = true then generate keyed arrays
    voltageadjointout = if returnvoltageadjoint && keyedarrays
        nodevariabletokeyed(voltageadjoint, modes,
            nodenames, modes, portnumbers, w)
    else
        voltageadjoint
    end

    return LinearizedHB(w, modes, Sout, Snoiseout, Ssensitivityout, Zout,
        Zadjointout, Zsensitivityout, Zsensitivityadjointout, QEout,
        QEidealout, CMout, nodefluxout, nodefluxadjointout, voltageout,
        voltageadjointout, nodenames, nodeindices, componentnames,
        componenttypes, componentnamedict, mutualinductorbranchnames,
        portnumbers, portindices, portimpedanceindices,
        noiseportimpedanceindices, sensitivitynames, sensitivityindices,
        Nsignalmodes, Nnodes, Nbranches, Nports, signalindex)
end

"""
    hblinsolve_inner!(S, Snoise, QE, CM, nodeflux, voltage,
        lsys, bnm,
        portindices, portimpedanceindices, noiseportimpedanceindices,
        portimpedances, noiseportimpedances, nodeindices, componenttypes,
        w, indices, wp, Nmodes, Nnodes, symfreqvar, wi, factorization)

Solve the linearized harmonic balance problem for a subset of the frequencies
given by `wi`, assembling the per-frequency system matrices from the
[`HBLinearizedSystem`](@ref) `lsys` with [`assemblesystemmatrix!`](@ref).
This function is thread safe in that different frequencies can be computed
in parallel on separate threads; `lsys` is only read.
"""
function hblinsolve_inner!(S, Snoise, Ssensitivity, Z, Zadjoint, Zsensitivity,
    Zsensitivityadjoint, QE, CM, nodeflux, nodefluxadjoint, voltage,
    voltageadjoint, lsys, bnm,
    portindices, portimpedanceindices, noiseportimpedanceindices,
    sensitivityindices, portimpedances, noiseportimpedances, nodeindices,
    componenttypes, w, wpumpmodes, Nmodes, Nnodes, symfreqvar, wi, factorization)

    Nports = length(portindices)
    phin = zeros(Complex{Float64}, size(lsys.Asparse,1), Nmodes*Nports)
    # inputwave = Diagonal(zeros(Complex{Float64}, Nports*Nmodes))
    inputwave = zeros(Complex{Float64}, Nports*Nmodes, Nports*Nmodes)
    outputwave = zeros(Complex{Float64}, Nports*Nmodes, Nports*Nmodes)

    Nnoiseports = length(noiseportimpedanceindices)
    noiseoutputwave = zeros(Complex{Float64}, Nnoiseports*Nmodes,
        Nports*Nmodes)
    sensitivityoutputvoltage = zeros(Complex{Float64},
        length(sensitivityindices)*Nmodes, Nports*Nmodes)

    # operate on a copy of the system matrix because it is modified per
    # frequency, potentially by multiple threads at the same time.
    Asparsecopy = copy(lsys.Asparse)

    # if using the KLU factorization and sparse solver then make a 
    # factorization for the sparsity pattern.
    cache = FactorizationCache()

    # if the scattering matrix is empty define a new working matrix
    if isempty(S)
        Sview = zeros(Complex{Float64}, Nports*Nmodes, Nports*Nmodes)
    end

    if isempty(Z)
        Zview = zeros(Complex{Float64}, Nports*Nmodes, Nports*Nmodes)
    end

    if isempty(Zadjoint)
        Zadjointview = zeros(Complex{Float64}, Nports*Nmodes, Nports*Nmodes)
    end

    # if the noise scattering matrix is empty define a new working matrix
    if isempty(Snoise)
        Snoiseview = zeros(Complex{Float64}, Nnoiseports*Nmodes, Nports*Nmodes)
    end

    # if the noise scattering matrix is empty define a new working matrix
    if isempty(Zsensitivity)
        Zsensitivityview = zeros(Complex{Float64},
            length(sensitivityindices)*Nmodes, Nports*Nmodes)
    end

    # if the noise scattering matrix is empty define a new working matrix
    if isempty(Zsensitivityadjoint)
        Zsensitivityadjointview = zeros(Complex{Float64},
            length(sensitivityindices)*Nmodes, Nports*Nmodes)
    end

    # loop over the frequencies
    for i in wi

        # if the scattering matrix is not empty define a view
        if !isempty(S)
            Sview = view(S, :, :, i)
        end

        if !isempty(Z)
            Zview = view(Z, :, :, i)
        end

        if !isempty(Zadjoint)
            Zadjointview = view(Zadjoint, :, :, i)
        end

        # if the noise scattering matrix is not empty define a view
        if !isempty(Snoise)
            Snoiseview = view(Snoise, :, :, i)
        end

        # if the noise scattering matrix is not empty define a view
        if !isempty(Zsensitivity)
            Zsensitivityview = view(Zsensitivity, :, :, i)
        end

        # if the noise scattering matrix is not empty define a view
        if !isempty(Zsensitivityadjoint)
            Zsensitivityadjointview = view(Zsensitivityadjoint, :, :, i)
        end

        # calculate the frequency matrices
        ws = w[i]
        # wmodes = calcw(ws,indices,wp);
        wmodes = ws .+ wpumpmodes
        # the repeat count covers the auxiliary variables of the modified
        # nodal analysis augmentation as well as the node fluxes.
        wouter = size(lsys.Asparse,1) ÷ length(wmodes)
        wmodesm = Diagonal(repeat(wmodes, outer = wouter));
        wmodes2m = Diagonal(repeat(wmodes.^2, outer = wouter));

        # assemble the linearized system matrix at this frequency,
        # Asparsecopy = (AoLjnm + invLnm + im.*Gnm*wmodesm - Cnm*wmodes2m),
        # in a way that doesn't allocate significant memory, taking the
        # complex conjugates of the negative frequency mode entries of the
        # linear term matrices and substituting any symbolic frequency
        # variables.
        assemblesystemmatrix!(Asparsecopy, lsys, wmodesm, wmodes2m)

        # factor the sparse matrix
        # factorklu!(cache, Asparsecopy)
        tryfactorize!(cache, factorization, Asparsecopy)

        # solve the linear system
        trysolve!(phin, cache.factorization, bnm)

        # convert to node voltages. node flux is defined as the time integral
        # of node voltage so node voltage is derivative of node flux which can
        # be accomplished in the frequency domain by multiplying by j*w.
        if !isempty(voltage)
            @views voltage[:,:,i] .= im .* wmodesm.diag[1:size(voltage,1)] .* phin[1:size(voltage,1),:]
        end

        # copy the nodeflux for output. the auxiliary variables of the
        # modified nodal analysis augmentation are internal.
        if !isempty(nodeflux)
            copy!(view(nodeflux,:,:,i), view(phin, 1:size(nodeflux,1), :))
        end

        # calculate the scattering parameters
        if !isempty(S) || !isempty(QE) || !isempty(CM)
            calcinputoutput!(inputwave, outputwave, phin, bnm,
                portimpedanceindices, portimpedanceindices, portimpedances,
                portimpedances, nodeindices, componenttypes, wmodes, symfreqvar)
            calcscatteringmatrix!(Sview, inputwave, outputwave)
        end

        if !isempty(Z)
            calcinputcurrentoutputvoltage!(inputwave, outputwave, phin, bnm,
                portimpedanceindices, portimpedanceindices, nodeindices, wmodes)
            calcscatteringmatrix!(Zview, inputwave, outputwave)
        end

        if !isempty(Zsensitivity)
            calcinputcurrentoutputvoltage!(inputwave, sensitivityoutputvoltage,
                phin, bnm, portimpedanceindices, sensitivityindices,
                nodeindices, wmodes)
            calcscatteringmatrix!(Zsensitivityview, inputwave,
                sensitivityoutputvoltage)
        end

        if (Nnoiseports > 0 || !isempty(nodefluxadjoint) || !isempty(voltageadjoint) || !isempty(Zsensitivityadjoint) || !isempty(Zadjoint)) && (!isempty(Snoise) || !isempty(QE) || !isempty(CM) || !isempty(nodefluxadjoint) || !isempty(voltageadjoint) || !isempty(Zsensitivityadjoint) || !isempty(Zadjoint))

            # solve the linear system with the complex conjugate of the
            # pump modulation matrix
            assemblesystemmatrix!(Asparsecopy, lsys, wmodesm, wmodes2m;
                conjugatepump = true)

            # factor the sparse matrix
            tryfactorize!(cache, factorization, Asparsecopy)

            # solve the linear system
            trysolve!(phin, cache.factorization, bnm)

            # copy the nodeflux adjoint for output
            if !isempty(nodefluxadjoint)
                copy!(view(nodefluxadjoint,:,:,i), view(phin, 1:size(nodefluxadjoint,1), :))
            end

            if !isempty(voltageadjoint)
                @views voltageadjoint[:,:,i] .= im .* wmodesm.diag[1:size(voltageadjoint,1)] .* phin[1:size(voltageadjoint,1),:]
            end

            # calculate the noise scattering parameters
            if !isempty(Snoise)  || !isempty(QE) || !isempty(CM)
                calcinputoutputnoise!(inputwave, noiseoutputwave, phin, bnm,
                    portimpedanceindices, noiseportimpedanceindices,
                    portimpedances, noiseportimpedances, nodeindices,
                    componenttypes, wmodes, symfreqvar)
                calcscatteringmatrix!(Snoiseview, inputwave, noiseoutputwave)
            end

            if !isempty(Zadjoint)
                calcinputcurrentoutputvoltage!(inputwave, outputwave,
                    phin, bnm, portimpedanceindices, portimpedanceindices,
                    nodeindices, wmodes)
                calcscatteringmatrix!(Zadjointview, inputwave, outputwave)
            end

            if !isempty(Zsensitivityadjoint)
                calcinputcurrentoutputvoltage!(inputwave,
                    sensitivityoutputvoltage, phin, bnm, portimpedanceindices,
                    sensitivityindices, nodeindices, wmodes)
                calcscatteringmatrix!(Zsensitivityadjointview, inputwave,
                    sensitivityoutputvoltage)
            end

            # calculate the quantum efficiency
            if !isempty(QE)
                calcqe!(view(QE,:,:,i), Sview, transpose(Snoiseview))
            end

            # calculate the commutation relations (Manley-Rowe relations)
            if !isempty(CM)
                calccm!(view(CM,:,i), Sview, transpose(Snoiseview), wmodes)
            end
        else
            # calculate the quantum efficiency
            if !isempty(QE)
                calcqe!(view(QE,:,:,i), Sview)
            end

            # calculate the commutation relations (Manley-Rowe relations)
            if !isempty(CM)
                calccm!(view(CM,:,i), Sview, wmodes)
            end
        end
    end
    return nothing
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
- `ftol = 1e-8`: the function tolerance defined we considered converged,
    defined as norm(F)/norm(x) < ftol or norm(F,Inf) <= ftol.
- `method = :quasinewton`: the nonlinear solution method. `:quasinewton`
    iterates with the complex holomorphic Jacobian `Jx` only, an
    approximation to the exact Jacobian. `:newton` solves the equivalent real
    system with the exact real Jacobian, which restores quadratic convergence
    near the solution, including for multi-tone problems: its mode coupling
    indices include the couplings which alias back onto the sampled grid,
    matching the cyclic Fourier transforms of the residual, and the
    assembled Jacobian agrees with the matrix-free Jacobian-vector product
    of `HBSystem` to machine precision. Both Jacobians are assembled
    directly from the Fourier coefficients of `cos(phi(t))` with precomputed
    plans (see [`plancomplexjacobian`](@ref) and
    [`planrealjacobian`](@ref)).
- `andersondepth::Integer = method == :newton ? 0 : 5`: the depth of the
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
    method = :quasinewton, andersondepth::Integer = method == :newton ? 0 : 5,
    symfreqvar = nothing, sorting = :number, keyedarrays::Bool = true,
    sensitivitynames::Vector{String} = String[],
    factorization = KLUfactorization(), debugJacobian = false,
    ) where {N}

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
        sensitivitynames = sensitivitynames, factorization = factorization,
        debugJacobian = debugJacobian,
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
    method = :quasinewton, andersondepth::Integer = method == :newton ? 0 : 5,
    symfreqvar = nothing, keyedarrays::Bool = true,
    sensitivitynames::Vector{String} = String[],
    factorization = KLUfactorization(), debugJacobian = false,
    )

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
    phimatrix = zeros(Complex{Float64}, Nwtuple)

    # create an array to hold the time domain data for the RFFT. also generate
    # the plans.
    phimatrixtd, irfftplan, rfftplan = plan_applynl(phimatrix)

    # the number of frequency entries per Josephson junction in phimatrix
    Nfreq = prod(Nwtuple[1:end-1])

    x = if isnothing(x0)
        zeros(Complex{Float64}, (Nnodes-1)*Nmodes)
    else
        copy(x0)
    end
    F = zeros(Complex{Float64}, (Nnodes-1)*Nmodes)
    AoLjbmvector = zeros(Complex{Float64}, Nbranches*Nmodes)

    # make a sparse transpose (improves multiplication speed slightly)
    Rbnmt = sparse(transpose(Rbnm))

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
    # auxiliary branch current variables to the resistors, keeping their
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
    mnaindices = mnaresistorindices(componenttypes, vvn)
    Nauxr = length(mnaindices)*Nmodes
    Naux = Nauxr + length(coupledbranches)*Nmodes
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
    Amna = mnapad(Amna, length(coupledbranches)*Nmodes)
    AmnaL = calcAmnaind(coupledbranches, Lb, Mb, cg.Rbn, Nmodes,
        Nnodal + Nauxr, Nnodal + Naux, Lmean)
    Amna = spaddkeepzeros(Amna, AmnaL)
    # pad the system matrices and vectors with the auxiliary variables.
    # the incidence matrix gains empty columns so the branch fluxes and
    # the Josephson junction terms are unaffected, and the constant MNA
    # equations are folded into the static linear term, which is added
    # with unit coefficient and is its own contribution to the Jacobian.
    Rbnm = hcat(Rbnm, spzeros(eltype(Rbnm), size(Rbnm,1), Naux))
    Rbnmt = sparse(transpose(Rbnm))
    invLnm = spaddkeepzeros(mnapad(invLnm, Naux), Amna)
    Cnm = mnapad(Cnm, Naux)
    Gnm = mnapad(Gnm, Naux)
    bnm = vcat(bnm, zeros(eltype(bnm), Naux))
    wmodesm = Diagonal(repeat(wmodes,
        outer = Nnodes-1+length(mnaindices)+length(coupledbranches)))
    wmodes2m = Diagonal(repeat(wmodes.^2,
        outer = Nnodes-1+length(mnaindices)+length(coupledbranches)))
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
    Jx, complexjacobianplan = if method == :quasinewton || debugJacobian
        plancomplexjacobian(Amatrixindices, Ljb, Lmean, Rbnm, Nmodes,
            Nbranches, Nfreq, invLnm, Gnm, Cnm)
    else
        nothing, nothing
    end

    Jr, realjacobianplan = if method == :newton || debugJacobian
        planrealjacobian(Amatrixindicesaliased, Amatrixconjindices, Ljb,
            Lmean,
            Rbnm, Nmodes, Nbranches, Nfreq, invLnm, Gnm, Cnm, modelayout,
            modelayout)
    else
        nothing, nothing
    end

    # the unified evaluation object for the nonlinear system: the residual,
    # the matrix-free Jacobian-vector and Hessian-vector products, and the
    # assembled Jacobians are all evaluated through it, in the complex or
    # the equivalent real representation.
    sys = HBSystem(Rbnm, Rbnmt, invLnm, Gnm, Cnm, wmodesm, wmodes2m, bnm,
        Ljb, Ljbm, Lmean, freqindexmap, conjsourceindices,
        conjtargetindices, phimatrix, phimatrixtd, irfftplan, rfftplan,
        modelayout, realjacobianplan, complexjacobianplan)

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
        # write back the canonical real representation of the point
        complex_to_real!(xr, sys.x, modelayout.isreal)
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
        return (F=F, x=x, Fr=Fr, xr=xr, Jx=Jx, Jr=Jr, fj=fj!, fjreal=fjreal!,
            sys=sys, Nnodal=Nnodal, mnaindices=mnaindices,
            gaugeindices=gaugeindices, floatingcomponents=floatingcomponents,
            coupledbranches=coupledbranches,
            complexjacobianplan=complexjacobianplan,
            realjacobianplan=realjacobianplan, phimatrix=phimatrix,
            cosphimatrix=(x -> (setpoint!(sys, x);
                JosephsonCircuits._updatecosphimatrix!(sys);
                sys.phimatrix)), modelayout=modelayout,
            Amatrixindices=Amatrixindices,
            Amatrixindicesaliased=Amatrixindicesaliased,
            Amatrixconjindices=Amatrixconjindices, Ljb=Ljb, Ljbm=Ljbm,
            Lmean=Lmean, Rbnm=Rbnm, invLnm=invLnm, Gnm=Gnm, Cnm=Cnm,
            wmodesm=wmodesm, wmodes2m=wmodes2m, Nmodes=Nmodes,
            Nbranches=Nbranches, Nfreq=Nfreq)
    end

    # solve the nonlinear system
    info = if method == :quasinewton

        nlsolve!(fj!, F, Jx,x; iterations = iterations, ftol = ftol,
            andersondepth = andersondepth, factorization = factorization)

    elseif method == :newton

        # solve the equivalent real nonlinear system with the exact real
        # Jacobian then convert back to complex
        info = nlsolve!(fjreal!, Fr, Jr, xr; iterations = iterations,
            ftol = ftol, andersondepth = andersondepth,
            factorization = factorization)
        real_to_complex!(x,xr,modelayout.isreal)
        real_to_complex!(F,Fr,modelayout.isreal)
        info
    else
        throw(ArgumentError("Method $(method) is not defined."))
    end


    push!(solverstages, IterationInfo(info.label, 1.0,
        info.regularization, info.converged, info.iterations,
        info.normresidual, info.alpha, info.backtracks,
        info.andersonaccepted))

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
        residual!(F, sys)
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

    return NonlinearHB(w, frequencies, nodefluxout, Rbnmout, Ljb, Lb, Ljbm,
        Nmodes, Nbranches, nodenames, portnumbers, modes, Sout, solverinfo)

end

"""
    calcAoLjbmindices(Amatrixindices, Ljb::SparseVector, Nmodes, Nbranches,
        Nfreq)

Return the sparse matrix containing the indices from the frequency domain
RFFT data as well as the indices of the sparse matrix to conjugate.

# Examples
```jldoctest
Amatrixindices = [1 -2 -3 -4; 2 1 -2 -3; 3 2 1 -2; 4 3 2 1]
Ljb = JosephsonCircuits.SparseArrays.sparsevec([1,2],[1.0,1.0])
Nmodes = 4
Nbranches = length(Ljb)
Nfreq = 4
AoLjbmindices, conjindicessorted, nentries = JosephsonCircuits.calcAoLjbmindices(
    Amatrixindices,
    Ljb,
    Nmodes,
    Nbranches,
    Nfreq);
AoLjbmindices

# output
8×8 SparseArrays.SparseMatrixCSC{Int64, Int64} with 32 stored entries:
 1  2  3  4  ⋅  ⋅  ⋅  ⋅
 2  1  2  3  ⋅  ⋅  ⋅  ⋅
 3  2  1  2  ⋅  ⋅  ⋅  ⋅
 4  3  2  1  ⋅  ⋅  ⋅  ⋅
 ⋅  ⋅  ⋅  ⋅  5  6  7  8
 ⋅  ⋅  ⋅  ⋅  6  5  6  7
 ⋅  ⋅  ⋅  ⋅  7  6  5  6
 ⋅  ⋅  ⋅  ⋅  8  7  6  5
```
```jldoctest
Amatrixindices = [1 -2 -3 0; 2 1 -2 -3; 3 2 1 -2; 0 3 2 1]
Ljb = JosephsonCircuits.SparseArrays.sparsevec([1,2],[1.0,1.0])
Nmodes = 4
Nbranches = length(Ljb)
Nfreq = 4
AoLjbmindices, conjindicessorted, nentries = JosephsonCircuits.calcAoLjbmindices(
    Amatrixindices,
    Ljb,
    Nmodes,
    Nbranches,
    Nfreq);
AoLjbmindices

# output
8×8 SparseArrays.SparseMatrixCSC{Int64, Int64} with 28 stored entries:
 1  2  3  ⋅  ⋅  ⋅  ⋅  ⋅
 2  1  2  3  ⋅  ⋅  ⋅  ⋅
 3  2  1  2  ⋅  ⋅  ⋅  ⋅
 ⋅  3  2  1  ⋅  ⋅  ⋅  ⋅
 ⋅  ⋅  ⋅  ⋅  5  6  7  ⋅
 ⋅  ⋅  ⋅  ⋅  6  5  6  7
 ⋅  ⋅  ⋅  ⋅  7  6  5  6
 ⋅  ⋅  ⋅  ⋅  ⋅  7  6  5
```
```jldoctest
Amatrixindices = [1 -2 -3 -4; 2 1 -2 -3; 3 2 1 -2; 4 3 2 1]
Ljb = JosephsonCircuits.SparseArrays.sparsevec([1,3],[1.0,1.0])
Nmodes = 4
Nbranches = length(Ljb)
Nfreq = 4
AoLjbmindices, conjindicessorted, nentries = JosephsonCircuits.calcAoLjbmindices(
    Amatrixindices,
    Ljb,
    Nmodes,
    Nbranches,
    Nfreq);
for c in conjindicessorted;AoLjbmindices.nzval[c] = -AoLjbmindices.nzval[c];end;AoLjbmindices

# output
12×12 SparseArrays.SparseMatrixCSC{Int64, Int64} with 32 stored entries:
 1  -2  -3  -4  ⋅  ⋅  ⋅  ⋅  ⋅   ⋅   ⋅   ⋅
 2   1  -2  -3  ⋅  ⋅  ⋅  ⋅  ⋅   ⋅   ⋅   ⋅
 3   2   1  -2  ⋅  ⋅  ⋅  ⋅  ⋅   ⋅   ⋅   ⋅
 4   3   2   1  ⋅  ⋅  ⋅  ⋅  ⋅   ⋅   ⋅   ⋅
 ⋅   ⋅   ⋅   ⋅  ⋅  ⋅  ⋅  ⋅  ⋅   ⋅   ⋅   ⋅
 ⋅   ⋅   ⋅   ⋅  ⋅  ⋅  ⋅  ⋅  ⋅   ⋅   ⋅   ⋅
 ⋅   ⋅   ⋅   ⋅  ⋅  ⋅  ⋅  ⋅  ⋅   ⋅   ⋅   ⋅
 ⋅   ⋅   ⋅   ⋅  ⋅  ⋅  ⋅  ⋅  ⋅   ⋅   ⋅   ⋅
 ⋅   ⋅   ⋅   ⋅  ⋅  ⋅  ⋅  ⋅  5  -6  -7  -8
 ⋅   ⋅   ⋅   ⋅  ⋅  ⋅  ⋅  ⋅  6   5  -6  -7
 ⋅   ⋅   ⋅   ⋅  ⋅  ⋅  ⋅  ⋅  7   6   5  -6
 ⋅   ⋅   ⋅   ⋅  ⋅  ⋅  ⋅  ⋅  8   7   6   5
```
"""
function calcAoLjbmindices(Amatrixindices::Matrix, Ljb::SparseVector, Nmodes,
    Nbranches, Nfreq)

    # evaluate Amatrixindices to find the number of entries of each type
    nentries = 0
    nconjentries = 0
    nzeros = 0
    for j in 1:Nmodes
        for k in 1:Nmodes
            if Amatrixindices[j,k] == 0
                nzeros += 1
            else
                nentries += 1
                if Amatrixindices[j,k] < 0
                    nconjentries += 1
                end
            end
        end
    end

    # make these into a sparse matrix. skip any zeros
    conjindices = Vector{Int}(undef, nconjentries * nnz(Ljb))
    I = Vector{Int}(undef, nentries * nnz(Ljb))
    J = Vector{Int}(undef, nentries * nnz(Ljb))
    V = Vector{Int}(undef, nentries * nnz(Ljb))
    Vsort = Vector{Int}(undef, nentries * nnz(Ljb))

    # generate the contents of the sparse matrix 
    n = 1
    nconj = 1
    for i in 1:nnz(Ljb)
        for j in 1:Nmodes
            for k in 1:Nmodes
                if Amatrixindices[j,k] != 0
                    I[n] = j + (Ljb.nzind[i] - 1) * Nmodes
                    J[n] = k + (Ljb.nzind[i] - 1) * Nmodes
                    Vsort[n] = n
                    index = abs(Amatrixindices[j,k]) + Nfreq * (i - 1)
                    V[n] = index
                    if Amatrixindices[j,k] < 0
                        conjindices[nconj] = n
                        nconj += 1
                    end
                    n += 1
                end
            end
        end
    end

    # create the sparse matrix
    AoLjbmindices = sparse(I, J, Vsort, Nbranches * Nmodes, Nbranches * Nmodes)

    # find the sorting of nzvals in the sparse matrix and apply that same
    # sorting to 
    Vsort2 = copy(AoLjbmindices.nzval)
    conjindicessorted = Vsort2[conjindices]

    AoLjbmindices.nzval .= V[Vsort2]

    return AoLjbmindices, conjindicessorted, nentries

end

"""
    calcAoLjbm2(Am::Array, Amatrixindices::Matrix, Ljb::SparseVector, Lmean,
        Nmodes, Nbranches, Nfreq)

Return the harmonic balance matrix divided by the Josephson inductance.

# Examples
```jldoctest
Amatrix = ComplexF64[1.0 + 1.0im 1.0 + 1.0im; 1.0 + 1.0im 1.0 + 1.0im; 1.0 + 1.0im 1.0 + 1.0im]
Amatrixindices = [1 -2 -3; 2 1 -2; 3 2 1]
Ljb = JosephsonCircuits.SparseArrays.sparsevec([1,2],[1.0,2.0])
Lmean = 1
Nmodes = 3
Nbranches = 2
JosephsonCircuits.calcAoLjbm2(Amatrix, Amatrixindices, Ljb, Lmean, Nmodes, Nbranches)

# output
6×6 SparseArrays.SparseMatrixCSC{ComplexF64, Int64} with 18 stored entries:
 1.0+1.0im  1.0-1.0im  1.0-1.0im      ⋅          ⋅          ⋅    
 1.0+1.0im  1.0+1.0im  1.0-1.0im      ⋅          ⋅          ⋅    
 1.0+1.0im  1.0+1.0im  1.0+1.0im      ⋅          ⋅          ⋅    
     ⋅          ⋅          ⋅      0.5+0.5im  0.5-0.5im  0.5-0.5im
     ⋅          ⋅          ⋅      0.5+0.5im  0.5+0.5im  0.5-0.5im
     ⋅          ⋅          ⋅      0.5+0.5im  0.5+0.5im  0.5+0.5im
```
```jldoctest
Amatrix = ComplexF64[1.0 + 1.0im 1.0 + 1.0im; 1.0 + 1.0im 1.0 + 1.0im; 1.0 + 1.0im 1.0 + 1.0im]
Amatrixindices = [1 -2 0; 2 1 -2; 0 2 1]
Ljb = JosephsonCircuits.SparseArrays.sparsevec([1,2],[1.0,2.0])
Lmean = 1
Nmodes = 3
Nbranches = 2
JosephsonCircuits.calcAoLjbm2(Amatrix, Amatrixindices, Ljb, Lmean, Nmodes, Nbranches)

# output
6×6 SparseArrays.SparseMatrixCSC{ComplexF64, Int64} with 14 stored entries:
 1.0+1.0im  1.0-1.0im      ⋅          ⋅          ⋅          ⋅    
 1.0+1.0im  1.0+1.0im  1.0-1.0im      ⋅          ⋅          ⋅    
     ⋅      1.0+1.0im  1.0+1.0im      ⋅          ⋅          ⋅    
     ⋅          ⋅          ⋅      0.5+0.5im  0.5-0.5im      ⋅    
     ⋅          ⋅          ⋅      0.5+0.5im  0.5+0.5im  0.5-0.5im
     ⋅          ⋅          ⋅          ⋅      0.5+0.5im  0.5+0.5im
```
```jldoctest
@variables A11 A12 A21 A22 A31 A32 Lj1 Lj2
Amatrix = [A11 A12;A21 A22;A31 A32]
Amatrixindices = [1 -2 -3; 2 1 -2; 3 2 1]
Ljb = JosephsonCircuits.SparseArrays.sparsevec([1,2],[Lj1,Lj2])
Lmean = 1
Nmodes = 3
Nbranches = 2
JosephsonCircuits.calcAoLjbm2(Amatrix, Amatrixindices, Ljb, Lmean, Nmodes, Nbranches)

# output
6×6 SparseArrays.SparseMatrixCSC{Num, Int64} with 18 stored entries:
 A11 / Lj1  A21 / Lj1  A31 / Lj1          ⋅          ⋅          ⋅
 A21 / Lj1  A11 / Lj1  A21 / Lj1          ⋅          ⋅          ⋅
 A31 / Lj1  A21 / Lj1  A11 / Lj1          ⋅          ⋅          ⋅
         ⋅          ⋅          ⋅  A12 / Lj2  A22 / Lj2  A32 / Lj2
         ⋅          ⋅          ⋅  A22 / Lj2  A12 / Lj2  A22 / Lj2
         ⋅          ⋅          ⋅  A32 / Lj2  A22 / Lj2  A12 / Lj2
```
"""
function calcAoLjbm2(Am::Array, Amatrixindices::Matrix, Ljb::SparseVector,
    Lmean, Nmodes, Nbranches)

    Nfreq = prod(size(Am)[1:end-1])


    # calculate the sparse matrix filled with the indices of Am
    AoLjbmindices, conjindicessorted, Nfreq = calcAoLjbmindices(Amatrixindices,
        Ljb, Nmodes, Nbranches, Nfreq)

    # determine the type to use for AoLjbm
    type = promote_type(eltype(Am), eltype(1 ./Ljb.nzval))

    nzval = Vector{type}(undef, nnz(AoLjbmindices))

    AoLjbm = SparseMatrixCSC(AoLjbmindices.m, AoLjbmindices.n,
        AoLjbmindices.colptr, AoLjbmindices.rowval, nzval)

    updateAoLjbm2!(AoLjbm, Am, AoLjbmindices, conjindicessorted, Ljb, Lmean)

    return AoLjbm

end

"""
    updateAoLjbm2!(AoLjbm::SparseMatrixCSC, Am::Array, AoLjbmindices,
        conjindicessorted, Ljb::SparseVector, Lmean)

Update the values in the sparse AoLjbm matrix in place.

"""
function updateAoLjbm2!(AoLjbm::SparseMatrixCSC, Am::Array, AoLjbmindices,
    conjindicessorted, Ljb::SparseVector, Lmean)

    if nnz(AoLjbm) == 0
        nentries = 0
    else
        nentries = nnz(AoLjbm) ÷ nnz(Ljb)
    end

    # # does this run into problems if the inductance vector isn't sorted?
    # # can i guarantee that Ljb is always sorted? look into this
    # # copy over the values and scale by the inductance
    for i in eachindex(Ljb.nzval)
        for j in 1:nentries
            k = (i-1)*nentries+j
            AoLjbm.nzval[k] = Am[AoLjbmindices.nzval[k]] * (Lmean / Ljb.nzval[i])
        end
    end

    # for i in eachindex(AoLjbm.nzval)
    #     j = (i-1) ÷ (nentries) + 1
    #     AoLjbm.nzval[i] = Am[AoLjbmindices.nzval[i]] * (Lmean / Ljb.nzval[j])
    # end

    # take the complex conjugates
    for i in conjindicessorted
        AoLjbm.nzval[i] = conj(AoLjbm.nzval[i])
    end

    return nothing
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
