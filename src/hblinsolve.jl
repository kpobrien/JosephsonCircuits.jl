# =========================================================================
# The linearized harmonic balance frequency sweep.
# =========================================================================

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

# A fully numeric circuit (such as the output of a circuit builder) needs no
# component definitions, so `circuitdefs` is optional.
function hblinsolve(w, circuit; kwargs...)
    return hblinsolve(w, circuit, Dict{Any,Any}(); kwargs...)
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
    sensitivitypairs::AbstractVector =
        Tuple{String,Int,Complex{Float64}}[],
    sensitivityblockpairs::AbstractVector = Tuple{String,Int,Any}[],
    nsensitivityparameters::Integer = 0,
    sensitivitylabels::Union{Nothing,Vector{String}} = nothing,
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
    # no promoted components; see the matching comment in hbnlsolve
    mnaindices = Int[]
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

    # Design parameter sensitivities: the caller names physical parameters
    # rather than components, through the pair list (componentindex,
    # parameterindex, alpha) with alpha = (dc/dtheta)/c the direction of the
    # component value under the parameter. Every component whose value
    # depends on a parameter contributes to that parameter's derivative;
    # the stamps of one parameter merge into one contraction.
    Nlumpedpairs = length(sensitivitypairs)
    if !isempty(sensitivitypairs) || !isempty(sensitivityblockpairs)
        isempty(sensitivitynames) || throw(ArgumentError(
            "pass either sensitivitynames or sensitivitypairs, not both"))
        nsensitivityparameters > 0 || throw(ArgumentError(
            "sensitivitypairs requires nsensitivityparameters"))
        # pairs are keyed by component name because the parse sorts the
        # circuit, so positional indices from the caller would not survive.
        # A scattering block pair is keyed by its first port's component.
        sensitivitynames = vcat(
            String[String(t[1]) for t in sensitivitypairs],
            String[String(t[1]) for t in sensitivityblockpairs])
    end

    # find the indices associated with the components for which we will
    # calculate sensitivities
    sensitivityindices = zeros(Int,length(sensitivitynames))
    for i in eachindex(sensitivitynames)
        sensitivityindices[i] = componentnamedict[sensitivitynames[i]]
        # entries past the lumped pairs are the block pairs; in the legacy
        # sensitivitynames path there are none
        blocktail = !isempty(sensitivityblockpairs) && i > Nlumpedpairs
        if componenttypes[sensitivityindices[i]] == :S && !blocktail
            throw(ArgumentError(lazy"Sensitivities with respect to scattering block components require a block derivative; got $(sensitivitynames[i]) as a plain component."))
        end
        if blocktail && componenttypes[sensitivityindices[i]] != :S
            throw(ArgumentError(lazy"The block pair component $(sensitivitynames[i]) is not a scattering block port."))
        end
    end

    # the grouping of the pairs into merged stamps, and the design parameter
    # slot each group accumulates into. Computed here, before the residual
    # derivative columns and the operating point solves, because both must
    # be merged with the same grouping as the stamps. Each scattering block
    # pair is its own group: its stamp values are rebuilt per frequency and
    # cannot concatenate with anything.
    stampgrouping, stampslots = if isempty(sensitivitypairs) &&
            isempty(sensitivityblockpairs)
        Vector{Int}[], Int[]
    else
        g, sl = parametergrouping(sensitivitypairs,
            view(sensitivityindices, 1:Nlumpedpairs), componenttypes,
            Dict(idx => p
                for (p, idx) in enumerate(signalnm.portimpedanceindices)))
        for (bi, bp) in enumerate(sensitivityblockpairs)
            push!(g, [Nlumpedpairs + bi])
            push!(sl, Int(bp[2]))
        end
        g, sl
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
            substitutefreq.(vvn, symfreqvar, wmodes[1]))
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
    sensitivitystamps, sensitivityblockentries = if returnSsensitivity
        st = calcsensitivitystamps(
            sensitivityindices[1:(isempty(sensitivityblockpairs) ?
                end : Nlumpedpairs)],
            psc, cg, signalnm, lsys, phimatrix, mnaindices,
            coupledbranches, Nnodalmna, Nsignalmodes, Nnodes)
        if !isempty(sensitivitypairs)
            st = [reparameterize(st[k], sensitivitypairs[k][3],
                sensitivitypairs[k][2]) for k in eachindex(st)]
        end
        entries = Tuple{Int,Any,Any}[]
        if !isempty(sensitivityblockpairs)
            isnothing(ssys) && throw(ArgumentError(
                "the circuit has no scattering blocks to take a block sensitivity of"))
            for (bi, bp) in enumerate(sensitivityblockpairs)
                targetdef = psc.componentvalues[
                    sensitivityindices[Nlumpedpairs + bi]].block
                dsys, zsys = derivativestampsystems(ssys, targetdef, bp[3])
                push!(st, blocksensitivitystamp(ssys, Int(bp[2])))
                push!(entries, (0, dsys, zsys))
            end
        end
        if !isempty(stampgrouping)
            st = mergestamps(st, stampgrouping)
            # locate each block pair's stamp in the merged vector: block
            # groups are the singletons pointing past the lumped pairs
            for (gi, g) in enumerate(stampgrouping)
                if length(g) == 1 && g[1] > Nlumpedpairs
                    bi = g[1] - Nlumpedpairs
                    entries[bi] = (gi, entries[bi][2], entries[bi][3])
                end
            end
        end
        st, entries
    else
        SensitivityStamp[], Tuple{Int,Any,Any}[]
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
        if !isempty(sensitivitypairs) && !isnothing(sensitivitynodeflux)
            throw(ArgumentError("sensitivitynodeflux is per component and cannot be combined with sensitivitypairs; supply sensitivityresidual (or neither) instead."))
        end
    end
    # The residual derivative columns are per pair; the stamps of one group
    # merge into one contraction, so the columns must merge with the same
    # grouping before anything is sized or solved from them. Summing is
    # right because the residual is linear in each component's contribution.
    if !isnothing(sensitivityresidual) && !isempty(stampgrouping) &&
            length(stampgrouping) != size(sensitivityresidual, 2)
        sensitivityresidual = hcat(
            [sum(sensitivityresidual[:, g], dims = 2)
             for g in stampgrouping]...)
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
        isempty(sensitivityblockpairs) || throw(ArgumentError("The reverse mode sensitivity contraction does not support scattering block parameters; use sensitivitymode = :forward."))
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
        !isnothing(sensitivityresidual) &&
            isempty(sensitivityblockpairs) && Ncomponents >
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
            sensitivityresidual,
            isempty(stampslots) ? collect(1:size(sensitivityresidual, 2)) :
                stampslots)
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
        Ncomponents = isempty(sensitivitypairs) ?
            length(sensitivitynames) : nsensitivityparameters,
        Nnodes = Nnodes,
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
        reverse = sensitivityreverse,
        blockentries = sensitivityblockentries,
        patternindex = isempty(sensitivityblockentries) ? Int[] :
            ssys.patternindex)
    # the block sensitivity stamps are rebuilt per frequency on the host
    usedevice = !(backend isa CPU) && cansweepondevice(lsys) &&
        isempty(sensitivityblockentries)
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
            @warn "A scattering block of this circuit has scattering parameters which cannot be evaluated on the backend, so its noise channels are formed on the host and the whole adjoint solution is copied back at every frequency, which is the largest transfer in the sweep. Give the block's data as a constant matrix, as tabulated data, or as a callable of one entry at a time (`ScatteringParameters(f; nports, form = :entry)`) to keep the noise on the backend." maxlog=1
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
        Ssensitivitytokeyed(outputarrays.Ssensitivity, modes, portnumbers,
            modes, portnumbers,
            isnothing(sensitivitylabels) ? sensitivitynames :
                sensitivitylabels, w)
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
    # this worker's private scattering block sensitivity state, or nothing
    blocksens::Any
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
        zeros(Float64, Nnoisechannels*Nmodes),
        isempty(sensitivity.blockentries) ? nothing :
            WorkerBlockSensitivity(sensitivity.stamps,
                sensitivity.blockentries, sensitivity.patternindex))
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

    # this worker's private scattering block sensitivity state; bound
    # before the loop because `ws` is rebound to the signal frequency there
    blocksens = ws.blocksens

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
                # a scattering block stamp depends on the signal frequency
                # through S itself, so each worker rebuilds its private
                # copy here rather than rescaling a constant
                stamps = if isnothing(blocksens)
                    sensitivity.stamps
                else
                    refreshblockstamps!(blocksens, wmodes)
                    blocksens.stamps
                end
                calcSsensitivity!(view(arrays.Ssensitivity,:,:,:,i),
                    stamps, sensitivity.dAop, dAsparse, dAphin,
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

