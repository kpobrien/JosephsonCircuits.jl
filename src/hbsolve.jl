# Global debug buffer
const DEBUG_MESSAGES = String[]

function debug_log(msg::String)
    push!(DEBUG_MESSAGES, msg)
end

function get_debug_log()
    return copy(DEBUG_MESSAGES)
end

function clear_debug_log()
    empty!(DEBUG_MESSAGES)
end


# Warning buffer (separate from debug)
const WARNING_MESSAGES = String[]

function warning_log(msg::String)
    push!(WARNING_MESSAGES, msg)
end

function get_warning_log()
    return copy(WARNING_MESSAGES)
end

function clear_warning_log()
    empty!(WARNING_MESSAGES)
end

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
    nonlinear_elements
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


# Add near the top of hbsolve.jl, after the struct definitions

"""
    NonlinearElement

Represents a nonlinear element with its type and parameters.
"""
struct NonlinearElement
    type::Symbol  # :josephson, :taylor
    indices::Vector{Int}  # Branch indices
    params::Dict{Symbol,Any}  # Parameters for the nonlinearity
end


#     hbsolve(ws, wp, Ip, Nsignalmodes::Int, Npumpmodes::Int, circuit,
#         circuitdefs; pumpports = [1], iterations = 1000, ftol = 1e-8,
#         switchofflinesearchtol = 1e-5, alphamin = 1e-4,
#         symfreqvar = nothing, nbatches = Base.Threads.nthreads(), sorting = :number,
#         returnS = true, returnSnoise = false, returnQE = true, returnCM = true,
#         returnnodeflux = false, returnvoltage = false, returnnodefluxadjoint = false,
#         returnvoltageadjoint = false, keyedarrays::Val{K} = Val(false),
#         sensitivitynames::Vector{String} = String[], returnSsensitivity = false,
#         returnZ = false, returnZadjoint = false,
#         returnZsensitivity = false, returnZsensitivityadjoint = false,
#         factorization = KLUfactorization())

# Calls the new harmonic balance solvers, [`hbnlsolve`](@ref) and
# [`hblinsolve`](@ref), which work for an arbitrary number of modes and ports),
# using an identical syntax to [`hbsolveold`](@ref), which only supports four
# wave mixing processes involving single strong tone and an arbitrary number of
# tone in the linearized solver. This function is primarily for testing the new
# solvers and is now deprecated.

# This function attempts to mimic [`hbsolveold`](@ref), but with the difference:
# The outputs of the linearized harmonic balance solver [`hblinsolve`](@ref) may
# not have the same ordering of signal modes as in [`hblinsolveold`](@ref). In
# [`hblinsolve`](@ref) the signal mode is always at index 1 and the location of
# the other modes can be found by inspecting the contents of `modes`.
function hbsolve(ws, wp, Ip, Nsignalmodes::Int, Npumpmodes::Int, circuit,
    circuitdefs; pumpports = [1], iterations = 1000, ftol = 1e-8,
    switchofflinesearchtol = 1e-5, alphamin = 1e-4,
    x0 = nothing,
    symfreqvar = nothing, nbatches = Base.Threads.nthreads(), sorting = :number,
    returnS = true, returnSnoise = false, returnQE = true, returnCM = true,
    returnnodeflux = false, returnvoltage = false, returnnodefluxadjoint = false,
    returnvoltageadjoint = false, keyedarrays::Val{K} = Val(false),
    sensitivitynames::Vector{String} = String[], returnSsensitivity = false,
    returnZ = false, returnZadjoint = false,
    returnZsensitivity = false, returnZsensitivityadjoint = false,
    factorization = KLUfactorization()) where K

    # Base.depwarn("""
    # Calls the new harmonic balance solvers, [`hbnlsolve`](@ref) and
    # [`hblinsolve`](@ref), which work for an arbitrary number of modes and ports),
    # using an identical syntax to [`hbsolveold`](@ref), which only supports four
    # wave mixing processes involving single strong tone and an arbitrary number of
    # tone in the linearized solver. This function is primarily for testing the new
    # solvers and is now deprecated. Please switch to the new syntax.
    #     """, :hbsolve)

    # solve the nonlinear system using the old syntax externally and the new
    # syntax internally
    w = (wp,)
    Nharmonics = (2*Npumpmodes,)

    # create the sources vector
    sources = [(mode = (1,), port = pumpports[1], current = Ip[1])]
    # if there are multiple currents and pumpports then add them as sources
    # they all have to have the same frequency when using the old interface
    @assert length(pumpports) == length(Ip)
    if length(pumpports) > 1
        for i in 2:length(pumpports)
            push!(sources,(mode = (1,), port = pumpports[i], current = Ip[i]))
        end
    end

    # calculate the frequency struct
    freq = removeconjfreqs(
        truncfreqs(
            calcfreqsrdft(Nharmonics),
            dc=false, odd=true, even=false, maxintermodorder=Inf,
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

    # identify nonlinear elements
    nonlinear_elements = identify_nonlinear_elements(
        psc.componenttypes,
        psc.componentvalues,
        psc.nodeindices,
        cg.edge2indexdict,
        circuitdefs 
    )

    
    nonlinear_elements::Dict{Int, NonlinearElement} = nonlinear_elements

    # solve the nonlinear problem
    nonlinear = hbnlsolve(w, sources, freq, indices, psc, cg, nm, nonlinear_elements;
        iterations = iterations, x0 = x0, ftol = ftol,
        switchofflinesearchtol = switchofflinesearchtol, alphamin = alphamin,
        symfreqvar = symfreqvar, keyedarrays = keyedarrays,
        factorization = factorization)

    # generate the signal modes
    signalfreq = truncfreqs(
        calcfreqsdft((Nsignalmodes,)),
        dc=true,odd=false,even=true,maxintermodorder=Inf,
    )

    # remove one of the signal modes if Nsignalmodes is even for compatibility
    # with old harmonic balance solver
    if mod(Nsignalmodes,2) == 0 && Nsignalmodes > 0
        signalfreq = JosephsonCircuits.removefreqs(
            signalfreq,
            [(Nsignalmodes,)],
        )
    end

    # solve the linearized problem
    # i should make this a tuple
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
    hbsolve(ws, wp::NTuple{N,Number}, sources::Vector,
        Nmodulationharmonics::NTuple{M,Int}, Npumpharmonics::NTuple{N,Int},
        circuit, circuitdefs;dc = false, threewavemixing = false,
        fourwavemixing = true, maxintermodorder=Inf, iterations = 1000,
        ftol = 1e-8, switchofflinesearchtol = 1e-5, alphamin = 1e-4,
        symfreqvar = nothing, nbatches = Base.Threads.nthreads(),
        sorting = :number, returnS = true, returnSnoise = false, returnQE = true,
        returnCM = true, returnnodeflux = false, returnvoltage = false,
        returnnodefluxadjoint = false, returnvoltageadjoint = false,
        keyedarrays::Val{K} = Val(true), sensitivitynames::Vector{String} = String[],
        returnSsensitivity = false, returnZ = false, returnZadjoint = false,
        returnZsensitivity = false, returnZsensitivityadjoint = false,
        factorization = KLUfactorization()) where {N,M,K}

Calls the harmonic balance solvers, [`hbnlsolve`](@ref) and
[`hblinsolve`](@ref), which work for an arbitrary number of modes and ports,
and for both three and four wave mixing processes. See also [`hbnlsolve`](@ref)
and [`hblinsolve`](@ref).

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
- `switchofflinesearchtol = 1e-5`: the function tolerance at which we switch
    from Newton with linesearch to only Newton. For easily converging
    functions, setting this to zero can speed up simulations.
- `alphamin = 1e-4`: the minimum step size relative to 1 for the linesearch.
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
- `keyedarrays::Val{K} = Val(true)`: when Val(true) return the output matrices
    and vectors as keyed arrays for more intuitive indexing. When Val(false)
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
    circuit, circuitdefs;dc = false, threewavemixing = false,
    fourwavemixing = true, maxintermodorder=Inf, iterations = 1000,
    ftol = 1e-8, switchofflinesearchtol = 1e-5, alphamin = 1e-4,
    x0 = nothing,
    symfreqvar = nothing, nbatches = Base.Threads.nthreads(),
    sorting = :number, returnS = true, returnSnoise = false, returnQE = true,
    returnCM = true, returnnodeflux = false, returnvoltage = false,
    returnnodefluxadjoint = false, returnvoltageadjoint = false,
    keyedarrays::Val{K} = Val(true), sensitivitynames::Vector{String} = String[],
    returnSsensitivity = false, returnZ = false, returnZadjoint = false,
    returnZsensitivity = false, returnZsensitivityadjoint = false,
    factorization = KLUfactorization()) where {N,M,K}    

    # calculate the Frequencies struct
    freq = removeconjfreqs(
        truncfreqs(
            calcfreqsrdft(Npumpharmonics),
            dc=dc, odd=fourwavemixing, even=threewavemixing,
            maxintermodorder=maxintermodorder,
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

    # identify nonlinear elements
    nonlinear_elements = identify_nonlinear_elements(
        psc.componenttypes,
        psc.componentvalues,
        psc.nodeindices,
        cg.edge2indexdict,
        circuitdefs 
    )
    
    nonlinear_elements::Dict{Int, NonlinearElement} = nonlinear_elements


    # solve the nonlinear problem
    nonlinear = hbnlsolve(wp, sources, freq, indices, psc, cg, nm, nonlinear_elements;
        iterations = iterations, x0 = x0, ftol = ftol,
        switchofflinesearchtol = switchofflinesearchtol, alphamin = alphamin,
        symfreqvar = symfreqvar, keyedarrays = keyedarrays,
        factorization = factorization)

    # generate the signal modes
    signalfreq = truncfreqs(
        calcfreqsdft(Nmodulationharmonics),
        dc=true, odd=threewavemixing, even=fourwavemixing,
        maxintermodorder=maxintermodorder,
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
- `keyedarrays::Val{K} = Val(true)`: when Val(true) return the output matrices
    and vectors as keyed arrays for more intuitive indexing. When Val(false)
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
    fourwavemixing=true, returnnodeflux=true, keyedarrays = Val(false))
isapprox(linearized.nodeflux,
    ComplexF64[9.901008591291e-12 - 6.40587007644028e-14im 2.164688307719963e-14 - 2.90852607344097e-16im 6.671563044645655e-14 - 8.585524364135119e-16im; 2.1633104519765224e-14 - 8.251861334047893e-16im 1.0099063486905209e-11 - 1.948847859339803e-13im -8.532003011745068e-15 + 3.234788465760295e-16im; 6.671648606599472e-14 + 7.892709980649199e-16im -8.53757633177974e-15 - 9.748395563374129e-17im 9.856580758892428e-12 + 5.859984004390703e-14im; 1.5888896262186103e-11 - 1.0303480614499543e-13im -2.557126237504446e-12 + 1.759201163407723e-14im -8.475819811683215e-12 + 5.3531443609574795e-14im; -2.5781681021577177e-13 + 4.757590640631487e-15im 2.36818731889176e-12 - 4.569646499606389e-14im 1.116372367616482e-13 - 2.039935997276492e-15im; -1.0210743447568219e-11 - 5.905490368441375e-14im 1.3377918536056493e-12 + 7.190105205618706e-15im 2.5392856657302323e-11 + 1.5143842454586225e-13im; 2.4781693042536835e-11 - 1.6057018472176702e-13im -2.5342360504077476e-12 + 1.7306764301173096e-14im -8.40554044664581e-12 + 5.269404591748149e-14im; -2.348528974341763e-13 + 3.949450668269274e-15im 1.1449271118157543e-11 - 2.2093702114766968e-13im 1.0261871618968225e-13 - 1.7240213938923877e-15im; -1.0140560031409567e-11 - 5.828587508192886e-14im 1.3288225860409326e-12 + 7.0954601524623594e-15im 3.423954321087654e-11 + 2.0403371894291513e-13im],
    atol = 1e-6)

# output
true
```
"""
function hblinsolve(w, circuit,circuitdefs; Nmodulationharmonics = (0,),
    nonlinear=nothing, symfreqvar=nothing, threewavemixing=false,
    fourwavemixing=true, maxintermodorder=Inf,
    nbatches::Integer = Base.Threads.nthreads(), sorting = :number, returnS = true,
    returnSnoise = false, returnQE = true, returnCM = true,
    returnnodeflux = false, returnnodefluxadjoint = false,
    returnvoltage = false, returnvoltageadjoint = false,
    keyedarrays::Val{K} = Val(true), sensitivitynames::Vector{String} = String[],
    returnSsensitivity = false, returnZ = false, returnZadjoint = false,
    returnZsensitivity = false, returnZsensitivityadjoint = false,
    factorization = KLUfactorization()) where K

    # parse and sort the circuit
    psc = parsesortcircuit(circuit, sorting = sorting)

    # calculate the circuit graph
    cg = calccircuitgraph(psc)

    # generate the signal modes
    signalfreq = truncfreqs(
        calcfreqsdft(Nmodulationharmonics),
        dc=true, odd=threewavemixing, even=fourwavemixing, 
            maxintermodorder=maxintermodorder,
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
signalfreq = JosephsonCircuits.truncfreqs(
    JosephsonCircuits.calcfreqsdft(Nmodulationharmonics),
    dc = true, odd = threewavemixing, even = fourwavemixing,
    maxintermodorder = Inf,
)
linearized = JosephsonCircuits.hblinsolve(ws, psc, cg, circuitdefs,
    signalfreq;nonlinear = nonlinear, returnnodeflux=true, keyedarrays = Val(false))
isapprox(linearized.nodeflux,
    ComplexF64[9.901008591291e-12 - 6.40587007644028e-14im 2.164688307719963e-14 - 2.90852607344097e-16im 6.671563044645655e-14 - 8.585524364135119e-16im; 2.1633104519765224e-14 - 8.251861334047893e-16im 1.0099063486905209e-11 - 1.948847859339803e-13im -8.532003011745068e-15 + 3.234788465760295e-16im; 6.671648606599472e-14 + 7.892709980649199e-16im -8.53757633177974e-15 - 9.748395563374129e-17im 9.856580758892428e-12 + 5.859984004390703e-14im; 1.5888896262186103e-11 - 1.0303480614499543e-13im -2.557126237504446e-12 + 1.759201163407723e-14im -8.475819811683215e-12 + 5.3531443609574795e-14im; -2.5781681021577177e-13 + 4.757590640631487e-15im 2.36818731889176e-12 - 4.569646499606389e-14im 1.116372367616482e-13 - 2.039935997276492e-15im; -1.0210743447568219e-11 - 5.905490368441375e-14im 1.3377918536056493e-12 + 7.190105205618706e-15im 2.5392856657302323e-11 + 1.5143842454586225e-13im; 2.4781693042536835e-11 - 1.6057018472176702e-13im -2.5342360504077476e-12 + 1.7306764301173096e-14im -8.40554044664581e-12 + 5.269404591748149e-14im; -2.348528974341763e-13 + 3.949450668269274e-15im 1.1449271118157543e-11 - 2.2093702114766968e-13im 1.0261871618968225e-13 - 1.7240213938923877e-15im; -1.0140560031409567e-11 - 5.828587508192886e-14im 1.3288225860409326e-12 + 7.0954601524623594e-15im 3.423954321087654e-11 + 2.0403371894291513e-13im],
    atol = 1e-6)

# output
true
"""
function hblinsolve(w, psc::ParsedSortedCircuit,
    cg::CircuitGraph, circuitdefs, signalfreq::Frequencies{N}; nonlinear = nothing,
    symfreqvar = nothing, nbatches::Integer = Base.Threads.nthreads(),
    returnS = true, returnSnoise = false, returnQE = true, returnCM = true,
    returnnodeflux = false, returnnodefluxadjoint = false, returnvoltage = false,
    returnvoltageadjoint = false, keyedarrays::Val{K} = Val(true),
    sensitivitynames::Vector{String} = String[], returnSsensitivity = false,
    returnZ = false, returnZadjoint = false,
    returnZsensitivity = false, returnZsensitivityadjoint = false,
    factorization = KLUfactorization()) where {N,K}

    Nsignalmodes = length(signalfreq.modes)
    # calculate the numeric matrices
    signalnm = numericmatrices(psc, cg, circuitdefs, Nmodes = Nsignalmodes)
    

    if isnothing(nonlinear)

        allpumpfreq = calcfreqsrdft((0,))
        Amatrixmodes, Amatrixindices = hbmatind(allpumpfreq, signalfreq)
        Nwtuple = NTuple{length(allpumpfreq.Nw)+1,Int}((allpumpfreq.Nw..., length(signalnm.Ljb.nzval)))
        phimatrix = ones(Complex{Float64}, Nwtuple)
        wpumpmodes = calcmodefreqs((0.0,),signalfreq.modes)

    else

        pumpfreq = nonlinear.frequencies

        Npumpmodes = length(pumpfreq.modes)
        
        # pumpfreq = JosephsonCircuits.calcfreqsrdft((2*Npumpmodes,))
        allpumpfreq = calcfreqsrdft(pumpfreq.Nharmonics)
        
        pumpindices = fourierindices(pumpfreq)        

        Amatrixmodes, Amatrixindices = hbmatind(allpumpfreq, signalfreq)    

        # calculate the dimensions of the array which holds the frequency
        # domain information for the fourier transform
        # Count all nonlinear elements
        num_nl_elements = length(nonlinear.nonlinear_elements)
        Nwtuple = NTuple{length(pumpfreq.Nw)+1,Int}((pumpfreq.Nw..., num_nl_elements))
        
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

        # Recreate all_nl_branches from nonlinear_elements
        all_nl_branches = SparseVector(cg.Nbranches, Int[], Float64[])
        for (branch, elem) in sort(collect(nonlinear.nonlinear_elements), by=x->x[1])
            push!(all_nl_branches.nzind, branch)
            push!(all_nl_branches.nzval, 1.0)
        end 
               

        # Create the repeated version
        all_nl_branches_m = SparseVector(cg.Nbranches * Nsignalmodes, Int[], Float64[])
        for idx in all_nl_branches.nzind
            for m in 1:Nsignalmodes
                push!(all_nl_branches_m.nzind, (idx-1)*Nsignalmodes + m)
                push!(all_nl_branches_m.nzval, 1.0)
            end
        end        

        # Create pump mode indices for extraction        
        pump_nl_indices = Int[]
        for idx in all_nl_branches.nzind
            for m in 1:Npumpmodes
                push!(pump_nl_indices, (idx-1)*Npumpmodes + m)
            end
        end
        

        # Use all_nl_branches_m instead of Ljbm for indices
        phivectortomatrix!(
            branchflux[pump_nl_indices], phimatrix,
            pumpindices.vectomatmap,
            pumpindices.conjsourceindices,
            pumpindices.conjtargetindices,
            length(all_nl_branches.nzval)
        )

        # Apply the linearization of nonlinearities
        # For Josephson: cos(φ), for Taylor: derivative of polynomial
        apply_nonlinearities!(phimatrix, phimatrixtd, nonlinear.nonlinear_elements,
                            :jacobian, irfftplan, rfftplan) 
                                   
                            

        wpumpmodes = calcmodefreqs(nonlinear.w,signalfreq.modes)
        
        
    end    


    # this is the first signal frequency. we will use it for various setup tasks
    wmodes = w[1] .+ wpumpmodes
    wmodesm = Diagonal(repeat(wmodes,outer=psc.Nnodes-1));
    wmodes2m = Diagonal(repeat(wmodes.^2,outer=psc.Nnodes-1));

    Nfreq = prod(signalfreq.Nw)

    # Modified to handle all nonlinear elements, not just Josephson
    if isnothing(nonlinear)
        # Use signalnm.Ljb for the no-nonlinear case
        AoLjbmindices, conjindicessorted = calcAoLjbmindices(
            Amatrixindices, signalnm.Ljb, Nsignalmodes, cg.Nbranches, Nfreq
        )
        
        AoLjbm = calcAoLjbm2(phimatrix, 
            Amatrixindices, signalnm.Ljb, 1, Nsignalmodes, cg.Nbranches
        )
    else        
                
        all_nl_branches_for_calc = SparseVector(cg.Nbranches, Int[], Float64[])

        # First, add all Josephson junctions in the order they appear in signalnm.Ljb
        for i in 1:length(signalnm.Ljb.nzind)
            branch = signalnm.Ljb.nzind[i]
            if haskey(nonlinear.nonlinear_elements, branch)
                push!(all_nl_branches_for_calc.nzind, branch)
                push!(all_nl_branches_for_calc.nzval, signalnm.Ljb.nzval[i])
            end
        end

        # Then, add any Taylor elements that aren't already included
        for (branch, elem) in sort(collect(nonlinear.nonlinear_elements), by=x->x[1])
            if elem.type == :taylor && !(branch in all_nl_branches_for_calc.nzind)
                push!(all_nl_branches_for_calc.nzind, branch)
                if haskey(elem.params, :inductance)
                    push!(all_nl_branches_for_calc.nzval, elem.params[:inductance])
                else
                    push!(all_nl_branches_for_calc.nzval, 1.0)  # Fallback
                end
            end
        end                        
        
        
        AoLjbmindices, conjindicessorted = calcAoLjbmindices(
            Amatrixindices, all_nl_branches_for_calc, Nsignalmodes, cg.Nbranches, Nfreq
        )                                
                
        AoLjbm = calcAoLjbm2(phimatrix, 
            Amatrixindices, all_nl_branches_for_calc, 1, Nsignalmodes, cg.Nbranches
        )                      
                

    end

    AoLjnm = signalnm.Rbnm'*AoLjbm*signalnm.Rbnm  


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
    # return bnm
    # if there is a symbolic frequency variable, then we need to redo the noise
    # port calculation because calcnoiseportimpedanceindices() can't tell if a
    # symbolic expression is complex. 
    if isnothing(symfreqvar)
        noiseportimpedanceindices = signalnm.noiseportimpedanceindices
    else
        noiseportimpedanceindices = calcnoiseportimpedanceindices(
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
    # Asparse = (AoLjnm + invLnmcopy + Gnmcopy + Cnmcopy)
    Asparse = spaddkeepzeros(
        spaddkeepzeros(
            spaddkeepzeros(AoLjnm,invLnmcopy),Gnmcopy),Cnmcopy)

    # make the index maps so we can efficiently add the sparse matrices
    # together without allocations or changing the sparsity structure. 
    Cnmindexmap = sparseaddmap(Asparse,Cnmcopy)
    Gnmindexmap = sparseaddmap(Asparse,Gnmcopy)
    invLnmindexmap = sparseaddmap(Asparse,invLnmcopy)
    AoLjnmindexmap = sparseaddmap(Asparse,AoLjnm)

    portimpedances = [vvn[i] for i in portimpedanceindices]
    noiseportimpedances = [vvn[i] for i in noiseportimpedanceindices]

    # solve for Asparse once so we have something reasonable to
    # factorize.
    fill!(Asparse.nzval, zero(eltype(Asparse.nzval)))
    sparseadd!(Asparse, 1, AoLjnm, AoLjnmindexmap)

    # take the complex conjugate of the negative frequency terms in
    # the capacitance and conductance matrices. substitute in the symbolic
    # frequency variable if present.
    sparseaddconjsubst!(Asparse, -1, Cnm, wmodes2m, Cnmindexmap,
        real.(wmodesm) .< 0, wmodesm, Cnmfreqsubstindices, symfreqvar)
    sparseaddconjsubst!(Asparse, im, Gnm, wmodesm, Gnmindexmap,
        real.(wmodesm) .< 0, wmodesm, Gnmfreqsubstindices, symfreqvar)
    sparseaddconjsubst!(Asparse, 1, invLnm,
        Diagonal(ones(size(invLnm,1))), invLnmindexmap, real.(wmodesm) .< 0,
        wmodesm, invLnmfreqsubstindices, symfreqvar)

    # make arrays for the voltages, node fluxes, scattering parameters,
    # quantum efficiency, and commutatio relations. if we aren't returning a
    # matrix, set it to be an empty array.
    if returnS
        S = zeros(Complex{Float64}, Nports*Nsignalmodes, Nports*Nsignalmodes,
            length(w))
    else
        S = zeros(Complex{Float64},0,0,0)
    end

    if returnZ
        Z = zeros(Complex{Float64}, Nports*Nsignalmodes, Nports*Nsignalmodes,
            length(w))
    else
        Z = zeros(Complex{Float64},0,0,0)
    end

    if returnZadjoint
        Zadjoint = zeros(Complex{Float64}, Nports*Nsignalmodes,
            Nports*Nsignalmodes, length(w))
    else
        Zadjoint = zeros(Complex{Float64},0,0,0)
    end

    if returnSnoise
        Snoise = zeros(Complex{Float64},
            length(noiseportimpedanceindices)*Nsignalmodes,
            Nports*Nsignalmodes, length(w))
    else
        Snoise = zeros(Complex{Float64},0,0,0)
    end

    if returnSsensitivity
        Ssensitivity = zeros(Complex{Float64},
            length(sensitivitynames)*Nsignalmodes,
            Nports*Nsignalmodes, length(w))
    else
        Ssensitivity = zeros(Complex{Float64},0,0,0)
    end

    if returnZsensitivity
        Zsensitivity = zeros(Complex{Float64},
            length(sensitivitynames)*Nsignalmodes,
            Nports*Nsignalmodes, length(w))
    else
        Zsensitivity = zeros(Complex{Float64},0,0,0)
    end

    if returnZsensitivityadjoint
        Zsensitivityadjoint = zeros(Complex{Float64},
            length(sensitivitynames)*Nsignalmodes,
            Nports*Nsignalmodes, length(w))
    else
        Zsensitivityadjoint = zeros(Complex{Float64},0,0,0)
    end

    if returnQE
        QE = zeros(Float64,Nports*Nsignalmodes,Nports*Nsignalmodes,length(w))
    else
        QE = zeros(Float64,0,0,0)
    end

    if returnCM
        CM = zeros(Float64,Nports*Nsignalmodes,length(w))
    else
        CM = zeros(Float64,0,0)
    end

    if returnnodeflux
        nodeflux = zeros(Complex{Float64},Nsignalmodes*(Nnodes-1),
            Nsignalmodes*Nports, length(w))
    else
        nodeflux = Vector{Complex{Float64}}(undef,0)
    end

    if returnnodefluxadjoint
        nodefluxadjoint = zeros(Complex{Float64}, Nsignalmodes*(Nnodes-1),
            Nsignalmodes*Nports, length(w))
    else
        nodefluxadjoint = Vector{Complex{Float64}}(undef,0)
    end

    if returnvoltage
        voltage = zeros(Complex{Float64}, Nsignalmodes*(Nnodes-1),
            Nsignalmodes*Nports, length(w))
    else
        voltage = Vector{Complex{Float64}}(undef,0)
    end

    if returnvoltageadjoint
        voltageadjoint = zeros(Complex{Float64}, Nsignalmodes*(Nnodes-1),
            Nsignalmodes*Nports, length(w))
    else
        voltageadjoint = Vector{Complex{Float64}}(undef,0)
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
            QE, CM, nodeflux, nodefluxadjoint, voltage, voltageadjoint, Asparse,
            AoLjnm, invLnm, Cnm, Gnm, bnm,
            AoLjnmindexmap, invLnmindexmap, Cnmindexmap, Gnmindexmap,
            Cnmfreqsubstindices, Gnmfreqsubstindices, invLnmfreqsubstindices,
            portindices, portimpedanceindices, noiseportimpedanceindices, sensitivityindices,
            portimpedances, noiseportimpedances, nodeindices, componenttypes,
            w, wpumpmodes, Nsignalmodes, Nnodes, symfreqvar, batch, factorization)
    end

    # calculate the `ideal` quantum efficiency based on the gain assuming an
    # ideal two mode amplifier
    if returnQE
        # calculate the ideal quantum efficiency.
        QEideal = zeros(Float64,size(S))
        calcqeideal!(QEideal,S)
    else
        QEideal = zeros(Float64,0,0,0)
    end

    # turn all of the array outputs into keyed arrays if the
    # keyedarrays = Val(true)

    # if keyword argument keyedarrays = Val(true) then generate keyed arrays
    # for the scattering parameters
    if returnS && K
        Sout = Stokeyed(S, modes, portnumbers, modes, portnumbers, w)
    else
        Sout = S
    end

    if returnZ && K
        Zout = Stokeyed(Z, modes, portnumbers, modes, portnumbers, w)
    else
        Zout = Z
    end

    if returnZadjoint && K
        Zadjointout = Stokeyed(Zadjoint, modes, portnumbers, modes, portnumbers, w)
    else
        Zadjointout = Zadjoint
    end

    # if keyword argument keyedarrays = Val(true) then generate keyed arrays
    # for the noise scattering parameters
    if returnSnoise && K
        Snoiseout =  Snoisetokeyed(Snoise, modes,
            componentnames[noiseportimpedanceindices], modes, portnumbers, w)
    else
        Snoiseout = Snoise
    end

    # if keyword argument keyedarrays = Val(true) then generate keyed arrays
    # for Ssensitivity
    if returnSsensitivity && K
        Ssensitivityout =  Snoisetokeyed(Ssensitivity, modes,
            sensitivitynames, modes, portnumbers, w)
    else
        Ssensitivityout = Ssensitivity
    end

    if returnZsensitivity && K
        Zsensitivityout =  Snoisetokeyed(Zsensitivity, modes,
            sensitivitynames, modes, portnumbers, w)
    else
        Zsensitivityout = Zsensitivity
    end

    if returnZsensitivityadjoint && K
        Zsensitivityadjointout =  Snoisetokeyed(Zsensitivityadjoint, modes,
            sensitivitynames, modes, portnumbers, w)
    else
        Zsensitivityadjointout = Zsensitivityadjoint
    end

    # if keyword argument keyedarrays = Val(true) then generate keyed arrays
    # for the quantum efficiency
    if returnQE && K
        QEout = Stokeyed(QE, modes, portnumbers, modes, portnumbers, w)
        QEidealout = Stokeyed(QEideal, modes, portnumbers, modes, portnumbers, w)
    else
        QEout = QE
        QEidealout = QEideal
    end

    # if keyword argument keyedarrays = Val(true) then generate keyed arrays
    # for the commutation relations
    if returnCM && K
        CMout = CMtokeyed(CM, modes, portnumbers, w)
    else
        CMout = CM
    end

    # if keyword argument keyedarrays = Val(true) then generate keyed arrays
    if returnnodeflux && K
        nodefluxout = nodevariabletokeyed(nodeflux, modes, nodenames, modes,
            portnumbers, w)
    else
        nodefluxout = nodeflux
    end

    # if keyword argument keyedarrays = Val(true) then generate keyed arrays
    if returnnodefluxadjoint && K
        nodefluxadjointout = nodevariabletokeyed(nodefluxadjoint, modes,
            nodenames, modes, portnumbers, w)
    else
        nodefluxadjointout = nodefluxadjoint
    end

    # if keyword argument keyedarrays = Val(true) then generate keyed arrays
    if returnvoltage && K
        voltageout = nodevariabletokeyed(voltage, modes,
            nodenames, modes, portnumbers, w)
    else
        voltageout = voltage
    end

    # if keyword argument keyedarrays = Val(true) then generate keyed arrays
    if returnvoltageadjoint && K
        voltageadjointout = nodevariabletokeyed(voltageadjoint, modes,
            nodenames, modes, portnumbers, w)
    else
        voltageadjointout = voltageadjoint
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
    hblinsolve_inner!(S, Snoise, QE, CM, nodeflux, voltage, Asparse,
        AoLjnm, invLnm, Cnm, Gnm, bnm,
        AoLjnmindexmap, invLnmindexmap, Cnmindexmap, Gnmindexmap,
        Cnmfreqsubstindices, Gnmfreqsubstindices, invLnmfreqsubstindices,
        portindices, portimpedanceindices, noiseportimpedanceindices,
        portimpedances, noiseportimpedances, nodeindices, componenttypes,
        w, indices, wp, Nmodes, Nnodes, symfreqvar, wi, factorization)

Solve the linearized harmonic balance problem for a subset of the frequencies
given by `wi`. This function is thread safe in that different frequencies can
be computed in parallel on separate threads.
"""
function hblinsolve_inner!(S, Snoise, Ssensitivity, Z, Zadjoint, Zsensitivity,
    Zsensitivityadjoint, QE, CM, nodeflux, nodefluxadjoint, voltage,
    voltageadjoint, Asparse, AoLjnm, invLnm, Cnm, Gnm, bnm,
    AoLjnmindexmap, invLnmindexmap, Cnmindexmap, Gnmindexmap,
    Cnmfreqsubstindices, Gnmfreqsubstindices, invLnmfreqsubstindices,
    portindices, portimpedanceindices, noiseportimpedanceindices,
    sensitivityindices, portimpedances, noiseportimpedances, nodeindices,
    componenttypes, w, wpumpmodes, Nmodes, Nnodes, symfreqvar, wi, factorization)
    

    Nports = length(portindices)
    phin = zeros(Complex{Float64}, Nmodes*(Nnodes-1), Nmodes*Nports)
    # inputwave = Diagonal(zeros(Complex{Float64}, Nports*Nmodes))
    inputwave = zeros(Complex{Float64}, Nports*Nmodes, Nports*Nmodes)
    outputwave = zeros(Complex{Float64}, Nports*Nmodes, Nports*Nmodes)

    Nnoiseports = length(noiseportimpedanceindices)
    noiseoutputwave = zeros(Complex{Float64}, Nnoiseports*Nmodes,
        Nports*Nmodes)
    sensitivityoutputvoltage = zeros(Complex{Float64},
        length(sensitivityindices)*Nmodes, Nports*Nmodes)

    # operate on a copy of Asparse because it may be modified by multiple
    # threads at the same time.
    Asparsecopy = copy(Asparse)

    # calculate the conjugate of AoLjnm
    AoLjnmconj = copy(AoLjnm)
    conj!(AoLjnmconj.nzval)

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
        wmodesm = Diagonal(repeat(wmodes, outer = Nnodes-1));
        wmodes2m = Diagonal(repeat(wmodes.^2, outer = Nnodes-1));
        

        # perform the operation below in a way that doesn't allocate
        # significant memory, plus take the conjugates mentioned below.
        # Asparsecopy = (AoLjnm + invLnm - im.*Gnm*wmodesm - Cnm*wmodes2m)

        fill!(Asparsecopy.nzval, zero(eltype(Asparsecopy.nzval)))
        sparseadd!(Asparsecopy, 1, AoLjnm, AoLjnmindexmap)
        

        # take the complex conjugate of the negative frequency terms in
        # the capacitance and conductance matrices. substitute in the symbolic
        # frequency variable if present.
        sparseaddconjsubst!(Asparsecopy, -1, Cnm, wmodes2m, Cnmindexmap,
            real.(wmodesm) .< 0, wmodesm, Cnmfreqsubstindices, symfreqvar)
        sparseaddconjsubst!(Asparsecopy, im, Gnm, wmodesm, Gnmindexmap,
            real.(wmodesm) .< 0, wmodesm, Gnmfreqsubstindices, symfreqvar)
        sparseaddconjsubst!(Asparsecopy, 1, invLnm,
            Diagonal(ones(size(invLnm,1))), invLnmindexmap, real.(wmodesm) .< 0,
            wmodesm, invLnmfreqsubstindices, symfreqvar)        
                            
        

        # factor the sparse matrix
        # factorklu!(cache, Asparsecopy)
        tryfactorize!(cache, factorization, Asparsecopy)        

        # solve the linear system
        trysolve!(phin, cache.factorization, bnm)        

        # convert to node voltages. node flux is defined as the time integral
        # of node voltage so node voltage is derivative of node flux which can
        # be accomplished in the frequency domain by multiplying by j*w.
        if !isempty(voltage)
            voltage[:,:,i] .= im*wmodesm*phin
        end

        # copy the nodeflux for output
        if !isempty(nodeflux)
            copy!(view(nodeflux,:,:,i),phin)
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

            # solve the nonlinear system with the complex conjugate of the pump
            # modulation matrix
            fill!(Asparsecopy.nzval, zero(eltype(Asparsecopy.nzval)))
            sparseadd!(Asparsecopy, 1, AoLjnmconj, AoLjnmindexmap)

            # take the complex conjugate of the negative frequency terms in
            # the capacitance and conductance matrices. substitute in the symbolic
            # frequency variable if present. 
            sparseaddconjsubst!(Asparsecopy, -1, Cnm, wmodes2m, Cnmindexmap,
                real.(wmodesm) .< 0, wmodesm, Cnmfreqsubstindices, symfreqvar)
            sparseaddconjsubst!(Asparsecopy, im, Gnm, wmodesm, Gnmindexmap,
                real.(wmodesm) .< 0, wmodesm, Gnmfreqsubstindices, symfreqvar)
            sparseaddconjsubst!(Asparsecopy, 1, invLnm,
                Diagonal(ones(size(invLnm,1))),invLnmindexmap, real.(wmodesm) .< 0,
                wmodesm, invLnmfreqsubstindices, symfreqvar)

            # factor the sparse matrix
            tryfactorize!(cache, factorization, Asparsecopy)

            # solve the linear system
            trysolve!(phin, cache.factorization, bnm)

            # copy the nodeflux adjoint for output
            if !isempty(nodefluxadjoint)
                copy!(view(nodefluxadjoint,:,:,i), phin)
            end

            if !isempty(voltageadjoint)
                copy!(view(voltageadjoint,:,:,i), im*wmodesm*phin)
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
        x0 = nothing, ftol = 1e-8, switchofflinesearchtol = 1e-5,
        alphamin = 1e-4, symfreqvar = nothing, sorting= :number)

Harmonic balance solver supporting an arbitrary number of large signals
(strong tones or pumps) and arbitrary numbers of ports, sources, and drives
including direct current (zero frequency) or flux pumping using a current
source and a mutual inductor. Use `hblinsolve` to linearize the system of
equations about the operating point found with `hbnlsolve`.

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
- `switchofflinesearchtol = 1e-5`:
- `alphamin = 1e-4`: the function tolerance at which we switch
    from Newton with linesearch to only Newton. For easily converging
    functions, setting this to zero can speed up simulations.
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
    maxintermodorder = Inf, dc = false, odd = true, even = false,
    x0 = nothing, ftol = 1e-8, switchofflinesearchtol = 1e-5, alphamin = 1e-4,
    symfreqvar = nothing, sorting= :number, keyedarrays::Val{K} = Val(true),
    sensitivitynames::Vector{String} = String[],
    factorization = KLUfactorization()) where {N,K}    

    # calculate the frequency struct
    freq = removeconjfreqs(
        truncfreqs(
            calcfreqsrdft(Nharmonics),
            dc=dc, odd=odd, even=even, maxintermodorder=maxintermodorder,
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

    # identify nonlinear elements
    nonlinear_elements = identify_nonlinear_elements(
        psc.componenttypes,
        psc.componentvalues,
        psc.nodeindices,
        cg.edge2indexdict,
        circuitdefs 
    )
    nonlinear_elements::Dict{Int, NonlinearElement} = nonlinear_elements

    return hbnlsolve(w, sources, freq, indices, psc, cg, nm, nonlinear_elements;
        iterations = iterations, x0 = x0, ftol = ftol,
        switchofflinesearchtol = switchofflinesearchtol, alphamin = alphamin,
        symfreqvar = symfreqvar, keyedarrays = keyedarrays,
        sensitivitynames = sensitivitynames,
        factorization = factorization)
end

"""
    hbnlsolve(w::NTuple{N,Number}, sources, frequencies::Frequencies{N},
        indices::FourierIndices{N}, psc::ParsedSortedCircuit, cg::CircuitGraph,
        nm::CircuitMatrices; iterations = 1000, x0 = nothing,
        ftol = 1e-8, switchofflinesearchtol = 1e-5, alphamin = 1e-4,
        symfreqvar = nothing)

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
nonlinear_elements = JosephsonCircuits.identify_nonlinear_elements(psc.componenttypes, psc.componentvalues, cg.nodeindices, cg.edge2indexdict, circuitdefs)

out=hbnlsolve(
    (wp,),
    [
        (mode=(0,),port=1,current=Idc),
        (mode=(1,),port=1,current=Ip),
    ],
    frequencies, fi, psc, cg, nm, nonlinear_elements)
isapprox(out.nodeflux[:],
    ComplexF64[15.190314040027522 - 8.56492651167657e-24im, 2.991103820177504e-6 - 1.8501001011477133e-8im, -6.835392148510984 - 1.0356102442254259e-14im, 7.396422335315908e-6 - 4.5749403967992827e-8im, 6.835392148539885 - 1.0356102451770844e-14im, 1.008026285172782e-5 - 6.23498762664213e-8im],
    atol = 1e-6)

# output
true
```
"""
function hbnlsolve(w::NTuple{N,Number}, sources, frequencies::Frequencies{N},
    indices::FourierIndices{N}, psc::ParsedSortedCircuit, cg::CircuitGraph,
    nm::CircuitMatrices, nonlinear_elements::Dict{Int, NonlinearElement}; iterations = 1000, x0 = nothing,
    ftol = 1e-8, switchofflinesearchtol = 1e-5, alphamin = 1e-4,
    symfreqvar = nothing, keyedarrays::Val{K} = Val(true),
    sensitivitynames::Vector{String} = String[],
    factorization = KLUfactorization()) where {N,K}
    

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

    # generate the frequencies of the modes
    Nmodes = length(modes)
    wmodes = calcmodefreqs(w,modes)

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
    Cnm = nm.Cnm
    Gnm = nm.Gnm
    invLnm = nm.invLnm
    portindices = nm.portindices
    portnumbers = nm.portnumbers
    portimpedanceindices = nm.portimpedanceindices
    vvn = nm.vvn
    Lmean = nm.Lmean
    # if there are no inductors, then Lmean will be zero so set it to be one
    if iszero(Lmean)
        Lmean = one(eltype(Lmean))
    end
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

    # Create a sparse vector that includes ALL nonlinear elements
    all_nl_branches = SparseVector(Nbranches, Int[], Float64[])
    for (branch, elem) in sort(collect(nonlinear_elements), by=x->x[1])
        push!(all_nl_branches.nzind, branch)
        
        # Get the actual inductance value
        if elem.type == :josephson
            # Find the inductance in Ljb
            idx = findfirst(x -> x == branch, Ljb.nzind)
            if !isnothing(idx)                
                push!(all_nl_branches.nzval, Ljb.nzval[idx])
            else
                error("Josephson junction at branch $branch not found in Ljb")
            end
        elseif elem.type == :taylor
            # Use the inductance from params
            if haskey(elem.params, :inductance)
                push!(all_nl_branches.nzval, elem.params[:inductance])
            else
                error("Taylor element at branch $branch missing inductance parameter")
            end
        else
            error("Unknown nonlinear element type: $(elem.type)")
        end
    end

    # Recalculate Lmean based on ALL nonlinear elements (not just Josephson)
    if iszero(nm.Lmean) && nnz(all_nl_branches) > 0
        # Calculate mean inductance from all nonlinear inductances
        Lmean = sum(all_nl_branches.nzval) / nnz(all_nl_branches)
    else
        # Keep the existing Lmean
        # debug_log("Using original Lmean: $Lmean")
    end

    # Create the repeated version for all nonlinear branches
    all_nl_branches_m = SparseVector(Nbranches * Nmodes, Int[], Float64[])
    for idx in all_nl_branches.nzind
        for m in 1:Nmodes
            push!(all_nl_branches_m.nzind, (idx-1)*Nmodes + m)            
            push!(all_nl_branches_m.nzval, all_nl_branches.nzval[findfirst(==(idx), all_nl_branches.nzind)])
        end
    end

    # calculate the dimensions of the array which holds the frequency
    # domain information for the fourier transform
    # Nwtuple = NTuple{length(Nw)+1,Int}((Nw..., length(Ljb.nzval)))
    # Calculate the dimensions based on ALL nonlinear elements
    Nwtuple = NTuple{length(Nw)+1,Int}((Nw..., length(all_nl_branches.nzval)))

    # create an array to hold the frequency domain data for the
    # fourier transform
    phimatrix = zeros(Complex{Float64}, Nwtuple)

    # create an array to hold the time domain data for the RFFT. also generate
    # the plans.
    phimatrixtd, irfftplan, rfftplan = plan_applynl(phimatrix)

    # initializing this with zeros seems to cause problems
    # ideally i should initialize the vector of ones then convert to the
    # matrix.
    Amatrix = rand(Complex{Float64}, Nwtuple)
    Nfreq = prod(size(Amatrix)[1:end-1])

    
    # AoLjbmindices, conjindicessorted = calcAoLjbmindices(Amatrixindices, Ljb, Nmodes, Nbranches, Nfreq)

    AoLjbmindices, conjindicessorted = calcAoLjbmindices(Amatrixindices, all_nl_branches, Nmodes, Nbranches, Nfreq)

    # right now i redo the calculation of AoLjbmindices, conjindicessorted in
    # calcAoLjbm2
    AoLjbm = calcAoLjbm2(Amatrix, Amatrixindices, all_nl_branches, Lmean, Nmodes, Nbranches)
    AoLjbmcopy = calcAoLjbm2(Amatrix, Amatrixindices, all_nl_branches, Lmean, Nmodes, Nbranches)   
    

    # convert to a sparse node matrix. Note: I was having problems with type 
    # instability when i used AoLjbm here instead of AoLjbmcopy. 
    AoLjnm = transpose(Rbnm)*AoLjbmcopy*Rbnm    
   

    if isnothing(x0)
        x = zeros(Complex{Float64}, (Nnodes-1)*Nmodes)   
    else
        x = x0    
    end    



    F = zeros(Complex{Float64}, (Nnodes-1)*Nmodes)
    AoLjbmvector = zeros(Complex{Float64}, Nbranches*Nmodes)

    # make a sparse transpose (improves multiplication speed slightly)
    Rbnmt = sparse(transpose(Rbnm))

    # substitute in the mode frequencies for components which have frequency
    # defined symbolically.
    Cnm = freqsubst(Cnm, wmodes, symfreqvar)
    Gnm = freqsubst(Gnm, wmodes, symfreqvar)
    invLnm = freqsubst(invLnm, wmodes, symfreqvar)

    # scale the matrices for numerical reasons
    Cnm *= Lmean
    Gnm *= Lmean
    invLnm *= Lmean

    # Calculate an initial Jacobian in order to create the factorization object.
    # This need to have the same sparsity structure as the actual Jacobian. If
    # the numerical values are vastly different from the actual Jacobian this
    # can cause a singular value error in klu! when we attempt to reuse the
    # symbolic factorization. We perform the sparse matrix addition keeping
    # numerical zeros (the usual sparse matrix addition converts these to
    # structural zeros which would change the sparsity structure).

    # Always use the same Jacobian pattern for consistency
    # When there are no Josephson junctions, AoLjnm will just contain zeros
    J = spaddkeepzeros(spaddkeepzeros(spaddkeepzeros(AoLjnm, invLnm), im*Gnm), -Cnm)

    # make the arrays and datastructures we need for
    # the non-allocating sparse matrix multiplication.
    AoLjbmRbnm = AoLjbmcopy*Rbnm
    xbAoLjbmRbnm = fill(false, size(AoLjbmcopy, 1))
    AoLjnm = Rbnmt*AoLjbmRbnm
    xbAoLjnm = fill(false, size(Rbnmt, 1))

    # make the index maps so we can add the sparse matrices together without
    # memory allocations. 
    AoLjnmindexmap = sparseaddmap(J, AoLjnm)
    invLnmindexmap = sparseaddmap(J, invLnm)
    Gnmindexmap = sparseaddmap(J, Gnm)
    Cnmindexmap = sparseaddmap(J, Cnm)

    # build the function and Jacobian for solving the nonlinear system
    # build the function and Jacobian for solving the nonlinear system
    function fj!(F, J, x)
        calcfj2!(F, J, x, wmodesm, wmodes2m, Rbnm, Rbnmt, invLnm,
            Cnm, Gnm, bnm, Ljb, all_nl_branches_m, Nmodes,
            Nbranches, Lmean, AoLjbmvector, AoLjbm,
            AoLjnmindexmap, invLnmindexmap, Gnmindexmap, Cnmindexmap,
            AoLjbmindices, conjindicessorted,
            freqindexmap, conjsourceindices, conjtargetindices, phimatrix,
            AoLjnm, xbAoLjnm, AoLjbmRbnm, xbAoLjbmRbnm,
            phimatrixtd, irfftplan, rfftplan,
            nonlinear_elements  # important
        )
        return nothing
    end

    # # use this for debugging purposes to return the function and the
    # # Jacobian
    # fj!(F,J,x)
    # return (F,J,x)

    # Debug: Check dimensions before calling fj!
    # debug_log("Before fj! call:")
    # debug_log("  x length = $(length(x))")
    # debug_log("  F length = $(length(F))")
    # debug_log("  J size = $(size(J))")
    # debug_log("  Nnodes = $Nnodes, Nmodes = $Nmodes")
    # debug_log("  Nbranches = $Nbranches")

    # solve the nonlinear system
    nlsolve!(fj!, F, J, x; iterations = iterations, ftol = ftol,
        switchofflinesearchtol = switchofflinesearchtol,
        alphamin = alphamin, factorization = factorization)


    nodeflux = x    

    # calculate the scattering parameters for the pump
    Nports = length(portindices)        


    S = zeros(Complex{Float64}, Nports*Nmodes, Nports*Nmodes)
    inputwave = zeros(Complex{Float64}, Nports*Nmodes)
    outputwave = zeros(Complex{Float64}, Nports*Nmodes)
    portimpedances = [vvn[i] for i in portimpedanceindices]
    if !isempty(S)
        calcinputoutput!(inputwave, outputwave, nodeflux, bnm/Lmean,
            portimpedanceindices, portimpedanceindices, portimpedances,
            portimpedances, nodeindices, componenttypes, wmodes, symfreqvar)
        calcscatteringmatrix!(S, inputwave, outputwave)
    end

    if K
        nodefluxout = nodevariabletokeyed(nodeflux, modes, nodenames)
    else
        nodefluxout = nodeflux
    end

    #
    if !isempty(S) && K
        Sout = Stokeyed(S, modes, portnumbers, modes, portnumbers)
    else
        Sout = S
    end

    return NonlinearHB(w, frequencies, nodefluxout, Rbnm, Ljb, Lb, Ljbm,
        Nmodes, Nbranches, nodenames, portnumbers, modes, Sout, nonlinear_elements)

end

"""
    calcfj2!(F,J,phin,wmodesm,wmodes2m,Rbnm,invLnm,Cnm,Gnm,bm,Ljb,Ljbindices,
        Ljbindicesm,Nmodes,Lmean,AoLjbm)
        
Calculate the residual and the Jacobian. These are calculated with one
function in order to reuse as much as possible.

Leave off the type signatures on F and J because the solver will pass
`nothing` if it only wants to calculate F or J.
"""
function calcfj2!(F,
        J,
        nodeflux::AbstractVector,
        wmodesm::AbstractMatrix,
        wmodes2m::AbstractMatrix,
        Rbnm::AbstractMatrix,
        Rbnmt::AbstractMatrix,
        invLnm::AbstractMatrix,
        Cnm::AbstractMatrix,
        Gnm::AbstractMatrix,
        bnm::AbstractVector,
        Ljb::SparseVector,
        Ljbm::SparseVector,
        Nmodes::Int,
        Nbranches::Int,
        Lmean,
        AoLjbmvector::AbstractVector,
        AoLjbm, AoLjnmindexmap, invLnmindexmap, Gnmindexmap, Cnmindexmap,
        AoLjbmindices, conjindicessorted,
        freqindexmap, conjsourceindices, conjtargetindices, phimatrix,
        AoLjnm, xbAoLjnm, AoLjbmRbnm, xbAoLjbmRbnm,
        phimatrixtd, irfftplan, rfftplan,
        nonlinear_elements::Dict{Int, NonlinearElement},
        )

    # Convert from node flux to branch flux
    phib = Rbnm*nodeflux

    # Extract ALL nonlinear branches (not just Josephson)
    all_nl_branches = Int[]
    branch_to_nl_idx = Dict{Int,Int}()
    nl_idx = 0
    
    for (branch, elem) in sort(collect(nonlinear_elements), by=x->x[1])
        nl_idx += 1
        branch_to_nl_idx[branch] = nl_idx
        push!(all_nl_branches, branch)
    end
    
    # Build indices for ALL nonlinear elements
    all_nl_indices = Int[]
    for branch in all_nl_branches
        for m in 1:Nmodes
            push!(all_nl_indices, (branch-1)*Nmodes + m)
        end
    end

    # Create a sparse vector with inductances for ALL nonlinear elements
    # Determine the appropriate type for the inductances array
    all_nl_inductances = Float64[]
    all_nl_branches_for_sparse = Int[]

    for (branch, elem) in sort(collect(nonlinear_elements), by=x->x[1])
        push!(all_nl_branches_for_sparse, branch)
        if elem.type == :josephson
            # Use stored inductance or find in Ljb
            if haskey(elem.params, :inductance)                
                push!(all_nl_inductances, elem.params[:inductance])                
            else
                idx = findfirst(x -> x == branch, Ljb.nzind)
                if !isnothing(idx)            
                    push!(all_nl_inductances, Ljb.nzval[idx])
                else                    
                    push!(all_nl_inductances, Lmean)
                end
            end
        elseif elem.type == :taylor
            # Use the inductance stored in params
            if haskey(elem.params, :inductance)
                L0_value = elem.params[:inductance]
            else
                # Fallback to Lmean if inductance not found
                L0_value = Lmean
                # debug_log("  Taylor branch $branch: WARNING - no inductance found, using Lmean = $L0_value")
            end
            push!(all_nl_inductances, L0_value)
        else
            push!(all_nl_inductances, Lmean)
        end
    end
    
    # Create the sparse vector
    all_nl_Lb = SparseVector(Nbranches, all_nl_branches_for_sparse, all_nl_inductances)
    

    if !isnothing(F)
        if !isempty(all_nl_indices)
            # Extract flux values for ALL nonlinear elements
            phivector_all = phib[all_nl_indices]
            
            # Map to phimatrix - sized for ALL nonlinear elements
            phivectortomatrix!(phivector_all, phimatrix, freqindexmap,
                conjsourceindices, conjtargetindices, length(nonlinear_elements))
            
            # Apply nonlinearities through FFT (now handles mixed types!)
            apply_nonlinearities!(phimatrix, phimatrixtd, nonlinear_elements, 
                                :function, irfftplan, rfftplan)
            
            # Convert back to vector
            fill!(AoLjbmvector, 0)
            AoLjbmvectorview = view(AoLjbmvector, all_nl_indices)
            phimatrixtovector!(AoLjbmvectorview, phimatrix, freqindexmap,
                conjsourceindices, conjtargetindices, length(nonlinear_elements))
                
            # Scale by inductance for each element
            for (idx, nl_idx) in enumerate(all_nl_indices)
                branch = (nl_idx - 1) ÷ Nmodes + 1
                elem = nonlinear_elements[branch]
        

                # Get the appropriate inductance value for scaling
                if elem.type == :josephson
                    # Find this branch in Ljbm
                    ljbm_idx = findfirst(x -> x == nl_idx, Ljbm.nzind)
                    if !isnothing(ljbm_idx)
                        AoLjbmvectorview[idx] *= (Lmean/Ljbm.nzval[ljbm_idx])
                    end
                elseif elem.type == :taylor
                    # For Taylor, use the L0 value from params
                    L0 = get(elem.params, :inductance, Lmean)
                    AoLjbmvectorview[idx] *= (Lmean/L0)                    
                end
            end
        end
        
        # Calculate F
        F .= Rbnmt*AoLjbmvector .+ invLnm*nodeflux .+ im*Gnm*wmodesm*nodeflux .- Cnm*wmodes2m*nodeflux .- bnm
    end

    # Calculate the Jacobian
    if !isnothing(J)
        if !isempty(all_nl_indices)
            # Extract flux values for ALL nonlinear elements
            phivector_all = phib[all_nl_indices]
            
            # Map to phimatrix
            phivectortomatrix!(phivector_all, phimatrix, freqindexmap,
                conjsourceindices, conjtargetindices, length(nonlinear_elements))
            
            # Apply derivative of nonlinearities for Jacobian
            apply_nonlinearities!(phimatrix, phimatrixtd, nonlinear_elements,
                                :jacobian, irfftplan, rfftplan)                        
        
            
            # Update AoLjbm with ALL nonlinear contributions
            updateAoLjbm2!(AoLjbm, phimatrix, AoLjbmindices, conjindicessorted,
                all_nl_Lb, Lmean)            
        end
        
        # Convert to sparse node matrix        
        spmatmul!(AoLjbmRbnm, AoLjbm, Rbnm, xbAoLjbmRbnm)
        spmatmul!(AoLjnm, Rbnmt, AoLjbmRbnm, xbAoLjnm)

        # Build Jacobian
        fill!(J, 0)
        sparseadd!(J, AoLjnm, AoLjnmindexmap)
        sparseadd!(J, invLnm, invLnmindexmap)
        sparseadd!(J, im, Gnm, wmodesm, Gnmindexmap)
        sparseadd!(J, -1, Cnm, wmodes2m, Cnmindexmap)
    end
    
    return nothing
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
    # it looks like there was a bug here.
    # AoLjbmindices, conjindicessorted, Nfreq = calcAoLjbmindices(Amatrixindices,
    #     Ljb, Nmodes, Nbranches, Nfreq)
    AoLjbmindices, conjindicessorted, nentries = calcAoLjbmindices(Amatrixindices,
        Ljb, Nmodes, Nbranches, Nfreq)    


    # determine the type to use for AoLjbm
    type = promote_type(eltype(Am), eltype(1 ./Ljb.nzval))

    # Check if we have symbolic values
    if type <: Symbolic
        type = Any
    end

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
    

    # Copy over the values and scale by the inductance
    for i in eachindex(Ljb.nzval)
        for j in 1:nentries
            k = (i-1)*nentries+j            
            AoLjbm.nzval[k] = Am[AoLjbmindices.nzval[k]] * (Lmean / Ljb.nzval[i])
        end
    end

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
            
            # new code 20250710
            # portindex is a component index, we need to find the actual branch
            # port_nodes = nodeindices[:, portindex]
            # key = (port_nodes[1], port_nodes[2])
            # end new

            # check if the mode is in the dictionary of modes
            if haskey(modedict, mode)
                # if we find the mode and the port, set that branch in bbm
                # equal to the current scaled by the mean inductance and the
                # flux quantum.
                modeindex = modedict[mode]                
                # old code 20250710:
                key = (nodeindices[1, portindex], nodeindices[2, portindex])                                    
                
                # old code 20250710: branch-major ordering
                bbm[(edge2indexdict[key]-1)*Nmodes+modeindex] += Lmean*current/phi0    
                # new code 20250710: mode-major ordering 
                # bbm[(modeindex-1)*Nbranches+edge2indexdict[key]] += Lmean*current/phi0           
            else
                throw(ArgumentError("Source mode $(mode) not found."))
            end
        else
            throw(ArgumentError("Source port $(port) not found."))
        end
    end     

    return nothing
end



function apply_nonlinearities!(phimatrix, phimatrixtd, nonlinear_elements, 
                               mode::Symbol, irfftplan, rfftplan)
    
    # Create a mapping of which nonlinearity to apply to each column
    nl_functions = Vector{Function}(undef, length(nonlinear_elements))

    # Sort nonlinear_elements by branch index to ensure consistent ordering
    sorted_elements = sort(collect(nonlinear_elements), by=x->x[1])
    
    for (i, (branch, elem)) in enumerate(sorted_elements)
        if elem.type == :josephson
            nl_functions[i] = (mode == :function) ? sin : cos
        elseif elem.type == :taylor
            coeffs = elem.params[:coeffs]
            powers = elem.params[:powers]
            
            # Get the inductance for normalization from params
            # L0 = elem.params[:inductance]  # This should be 3.29e-10
            # Ic_equiv = phi0 / L0  # Equivalent critical current for normalization
            
            # Handle empty coefficients/powers
            if isempty(coeffs) || isempty(powers)
                nl_functions[i] = φ -> φ
                # debug_log("Warning: Taylor branch $branch has empty coefficients, using linear")
            else  # This is the Taylor expansion case
                # CHANGE THIS VALUE TO ADJUST WRAPPING
                wrap_boundary = Inf #π  # e.g., 2.0, π/2, 1000*π for no wrapping, etc.
                
                # Define the wrapping function once
                wrap_phase = function(x)
                    if isinf(wrap_boundary)
                        # No wrapping
                        return x
                    else
                        # Proper wrapping to [-wrap_boundary, wrap_boundary]
                        # Using mod to handle both positive and negative values correctly
                        return mod(x + wrap_boundary, 2*wrap_boundary) - wrap_boundary
                    end
                end
                
                if mode == :function
                    # already normalized: no need to divide by Ic_equiv
                    nl_functions[i] = φ -> begin
                        # Wrap φ
                        if isa(φ, Complex)
                            φ_wrapped = complex(wrap_phase(real(φ)), wrap_phase(imag(φ)))
                        else
                            φ_wrapped = wrap_phase(φ)
                        end
                        
                        result = zero(typeof(φ))
                        for (c, p) in zip(coeffs, powers)
                            result += (c) * φ_wrapped^p
                        end
                        return result
                    end
                else # :jacobian
                    # Derivative of polynomial with wrapping
                    nl_functions[i] = φ -> begin
                        # Wrap φ (same as above)
                        if isa(φ, Complex)
                            φ_wrapped = complex(wrap_phase(real(φ)), wrap_phase(imag(φ)))
                        else
                            φ_wrapped = wrap_phase(φ)
                        end
                        
                        # Compute derivative
                        result = zero(typeof(φ))
                        for (c, p) in zip(coeffs, powers)
                            if p > 0
                                result += p * (c) * φ_wrapped^(p-1)
                            end
                        end
                        return result
                    end
                end
            end            
        else
            # Default to linear
            nl_functions[i] = φ -> φ
            debug_log("Warning: Unknown nonlinearity type $(elem.type) for branch $branch")
        end
    end
    

    # Apply mixed nonlinearities
    applynl_mixed!(phimatrix, phimatrixtd, nl_functions, irfftplan, rfftplan)
end