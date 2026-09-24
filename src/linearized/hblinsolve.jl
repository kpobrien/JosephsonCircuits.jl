# =========================================================================
# The linearized harmonic balance frequency sweep.
# =========================================================================

"""
    hblinsolve(w, circuit, circuitdefs; Nmodulationharmonics = (0,),
        nonlinear = nothing, threewavemixing = false,
        fourwavemixing = true, maxharmonics = Nmodulationharmonics,
        maxintermodorder = Inf, nbatches = Base.Threads.nthreads(),
        returnS = true, returnSnoise = false,
        returnCnoise = false, returnQE = true, returnCM = true,
        returnnodeflux = false, returnnodefluxadjoint = false,
        returnvoltage = false, returnvoltageadjoint = false,
        keyedarrays = true, temperature = 0.0,
        sensitivitynames::AbstractVector = String[],
        sensitivitynodeflux = nothing, sensitivityresidual = nothing,
        sensitivitymode = :auto, returnSsensitivity = false,
        factorization = nothing, backend = CPU())

Sweep the weak signal frequencies `w` through the circuit linearized about
the operating point `nonlinear` found by [`hbnlsolve`](@ref), or through
the linear circuit when `nonlinear = nothing`. Any number of signal and
idler modes, ports and pumps is supported. The scattering parameters,
noise scattering parameters, quantum efficiency, commutation relations,
node fluxes and voltages, and sensitivities are computed on request.

The linearized system is solved in the same modified nodal analysis
formulation as the nonlinear one (see [`hbnlsolve`](@ref)), without gauge
fixing rows because a mode at (numerically) zero total frequency is
rejected with an `ArgumentError`; estimate a direct current limit from a
sequence of decreasing nonzero frequencies instead.

# Arguments
- `w`: the signal angular frequency or frequencies in radians per second,
    a real number or any iterable of them.
- `circuit`: a typed [`Circuit`](@ref) or a [`CompiledCircuit`](@ref). A
    `Circuit` is compiled with the default node ordering; for another one,
    pass `compile(circuit; sorting = ...)`.
- `circuitdefs`: a dictionary from the symbols or symbolic variables used
    as component values to their numerical values. Optional when every
    component value is numeric.

# Keywords
- `Nmodulationharmonics = (0,)`: how many harmonics of each pump to retain
    around the signal, which sets the signal and idler modes; `(0,)` is
    the signal alone.
- `nonlinear = nothing`: the [`NonlinearHB`](@ref) operating point to
    linearize about, or `nothing` for a linear circuit.
- `threewavemixing = false`: retain the odd pump harmonics around the
    signal, through which three wave mixing couples the modes.
- `fourwavemixing = true`: retain the even pump harmonics around the
    signal, through which four wave mixing couples the modes.
- `maxharmonics = Nmodulationharmonics`: an upper bound on the absolute
    harmonic index retained for each pump; see [`truncfreqs`](@ref).
- `maxintermodorder = Inf`: keep only the modes whose harmonic indices
    have an absolute sum of at most this order.
$(_DOC_NBATCHES)
$(_DOC_RETURNS)
$(_DOC_TEMPERATURE)
$(_DOC_SENSNAMES)
- `sensitivitynodeflux = nothing`: the derivatives of the pump operating
    point with respect to each component value, the columns of
    [`calcnodefluxsensitivity`](@ref), to include the shift of the
    operating point in the sensitivities.
- `sensitivityresidual = nothing`: the derivatives of the harmonic balance
    residual with respect to each component value, the columns of
    [`calcresidualsensitivity`](@ref), which serve the same purpose;
    required by the reverse contraction order and sufficient for either.
    [`hbsolve`](@ref) supplies these. When neither is given the pump
    operating point is held fixed.
$(_DOC_SENSMODE)
$(_DOC_SSENS)
- `factorization = nothing`: the factorization of the linearized system
    matrix at each frequency: [`KLUfactorization`](@ref),
    [`LUfactorization`](@ref), [`CUDSSFactorization`](@ref) on a device,
    or [`BlockFactorization`](@ref) for dense node blocks (`Nmodes`
    unknowns per node, [`SparseBlockFactorization`](@ref)). `nothing`
    chooses by the number of tones and the memory
    ([`linearizedfactorization`](@ref)): the sparse factorization for one
    tone, the block factorization in double for two or more when its
    factors fit in half the free memory of the backend. On a device the
    choice also picks the solver of the batch: a sparse factorization is
    solved by cuDSS, a `BlockFactorization` by the batched block
    factorization.
    The precision of the solutions is the factorization's:
    `BlockFactorization(precision = Float32)` refines single precision
    factors against the double residual to double accuracy, and
    `BlockFactorization(precision = Float32, refine = false)` solves each
    frequency entirely in single precision, for the cases where single
    precision scattering parameters are enough (see
    [`BlockFactorization`](@ref) for the measured accuracy); the outputs
    are returned in double either way.
$(_DOC_LINBACKEND)
- `symfreqvar = nothing`: deprecated, the parameter a frequency dependent
    value was written as an expression in. Write the value as a
    [`FrequencyDependent`](@ref) closure of the frequency instead.
- `returnZ`, `returnZadjoint`, `returnZsensitivity`,
    `returnZsensitivityadjoint`: removed; passing any of them warns.
    Compute impedances from the scattering parameters instead.

# Returns
- `LinearizedHB`: the linearized solution; see [`LinearizedHB`](@ref).

# Examples
```jldoctest
circuit = Circuit(
    [:p1 => Port(1; Z0 = :Rleft),
     :l1 => Inductor(:Lm),
     :l2 => Inductor(:Lm),
     :k1 => MutualInductor(:K1, :l1, :l2),
     :cc => Capacitor(:Cc),
     :jj3 => JosephsonJunction(:Lj),
     :jj4 => JosephsonJunction(:Lj),
     :cj => Capacitor(:Cj),
     :gnd => Ground()],
    [[(:p1, 1), (:l1, 1), (:cc, 1)],
     [(:cc, 2), (:l2, 1), (:jj4, 1), (:cj, 1)],
     [(:l2, 2), (:jj3, 1)],
     [(:p1, 2), (:l1, 2), (:jj3, 2), (:jj4, 2), (:cj, 2), (:gnd, 1)]])
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
    nonlinear = nonlinear, threewavemixing=false,
    fourwavemixing=true, returnnodeflux=true, keyedarrays = false)
isapprox(linearized.nodeflux,
    ComplexF64[9.901008591291e-12 - 6.40587007644028e-14im 2.164688307719963e-14 - 2.90852607344097e-16im 6.671563044645655e-14 - 8.585524364135119e-16im; 2.1633104519765224e-14 - 8.251861334047893e-16im 1.0099063486905209e-11 - 1.948847859339803e-13im -8.532003011745068e-15 + 3.234788465760295e-16im; 6.671648606599472e-14 + 7.892709980649199e-16im -8.53757633177974e-15 - 9.748395563374129e-17im 9.856580758892428e-12 + 5.859984004390703e-14im; 1.5888896262186103e-11 - 1.0303480614499543e-13im -2.557126237504446e-12 + 1.759201163407723e-14im -8.475819811683215e-12 + 5.3531443609574795e-14im; -2.5781681021577177e-13 + 4.757590640631487e-15im 2.36818731889176e-12 - 4.569646499606389e-14im 1.116372367616482e-13 - 2.039935997276492e-15im; -1.0210743447568219e-11 - 5.905490368441375e-14im 1.3377918536056493e-12 + 7.190105205618706e-15im 2.5392856657302323e-11 + 1.5143842454586225e-13im; 2.4781693042536835e-11 - 1.6057018472176702e-13im -2.5342360504077476e-12 + 1.7306764301173096e-14im -8.40554044664581e-12 + 5.269404591748149e-14im; -2.348528974341763e-13 + 3.949450668269274e-15im 1.1449271118157543e-11 - 2.2093702114766968e-13im 1.0261871618968225e-13 - 1.7240213938923877e-15im; -1.0140560031409567e-11 - 5.828587508192886e-14im 1.3288225860409326e-12 + 7.0954601524623594e-15im 3.423954321087654e-11 + 2.0403371894291513e-13im],
    atol = 1e-6)

# output
true
```
"""
function hblinsolve(w, circuit::CompilableCircuit,
    circuitdefs::AbstractDict = Dict{Symbol,Any}(); Nmodulationharmonics = (0,),
    threewavemixing::Bool = false, fourwavemixing::Bool = true,
    maxharmonics = Nmodulationharmonics, maxintermodorder = Inf,
    symfreqvar = nothing, kwargs...)

    # compile the circuit; the matrices are assembled by the method below
    # at the signal mode count
    psc = compile(circuit)
    # the deprecated symbolic frequency variable, in circuit/legacy.jl
    isnothing(symfreqvar) || (psc = frequencydependentcircuit(psc,
        circuitdefs, symfreqvar, :hblinsolve))

    signalfreq = signalfrequencies(Nmodulationharmonics; threewavemixing,
        fourwavemixing, maxintermodorder, maxharmonics)

    # every other keyword is the sweep's, with its defaults there
    return hblinsolve(w, psc, circuitdefs, signalfreq; kwargs...)
end


"""
    hblinsolve(w, psc::CompiledCircuit, circuitdefs,
        signalfreq::Frequencies; nonlinear = nothing, kwargs...)

The linearized sweep on an already compiled circuit `psc`, at the signal
mode set `signalfreq`. `circuitdefs` is the usual
dictionary; in its place the vector of resolved component values `nm.vvn` of
a [`CircuitMatrices`](@ref) already built for the circuit is accepted, which
is what [`hbsolve`](@ref) passes so that the values are resolved once; the sweep
`w` is a real number or any iterable of them ([`sweepfrequencies`](@ref)).
This is what the other methods call after building those; it takes every keyword
of the general method except the ones which describe the mode set
(`Nmodulationharmonics`, `threewavemixing`, `fourwavemixing`,
`maxharmonics`, `maxintermodorder`), plus the design parameter
sensitivity keywords described under [`hbsolve`](@ref) (`sensitivitypairs`,
`sensitivityblockpairs`, `nsensitivityparameters`, `sensitivitylabels`,
which only this method accepts) and `debuglsys = false`, which returns
the [`HBLinearizedSystem`](@ref) and its ingredients instead of solving,
for building reference implementations in tests.

# Examples
```jldoctest
circuit = Circuit(
    [:p1 => Port(1; Z0 = :Rleft),
     :l1 => Inductor(:Lm),
     :l2 => Inductor(:Lm),
     :k1 => MutualInductor(:K1, :l1, :l2),
     :cc => Capacitor(:Cc),
     :jj3 => JosephsonJunction(:Lj),
     :jj4 => JosephsonJunction(:Lj),
     :cj => Capacitor(:Cj),
     :gnd => Ground()],
    [[(:p1, 1), (:l1, 1), (:cc, 1)],
     [(:cc, 2), (:l2, 1), (:jj4, 1), (:cj, 1)],
     [(:l2, 2), (:jj3, 1)],
     [(:p1, 2), (:l1, 2), (:jj3, 2), (:jj4, 2), (:cj, 2), (:gnd, 1)]])
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
psc = JosephsonCircuits.compile(circuit)
nm = JosephsonCircuits.numericmatrices(psc, circuitdefs, Nmodes = Nmodes)
nonlinear = hbnlsolve(
    (wp,),
    [
        (mode=(0,),port=1,current=Idc),
        (mode=(1,),port=1,current=Ip),
    ],
    frequencies, fi, psc, nm)
signalfreq =JosephsonCircuits.truncfreqs(
    JosephsonCircuits.calcfreqsdft(Nmodulationharmonics),
    dc = true, odd = threewavemixing, even = fourwavemixing,
    maxintermodorder = Inf,
)
linearized = JosephsonCircuits.hblinsolve(ws, psc, circuitdefs,
    signalfreq;nonlinear = nonlinear, returnnodeflux=true, keyedarrays = false)
isapprox(linearized.nodeflux,
    ComplexF64[9.901008591291e-12 - 6.40587007644028e-14im 2.164688307719963e-14 - 2.90852607344097e-16im 6.671563044645655e-14 - 8.585524364135119e-16im; 2.1633104519765224e-14 - 8.251861334047893e-16im 1.0099063486905209e-11 - 1.948847859339803e-13im -8.532003011745068e-15 + 3.234788465760295e-16im; 6.671648606599472e-14 + 7.892709980649199e-16im -8.53757633177974e-15 - 9.748395563374129e-17im 9.856580758892428e-12 + 5.859984004390703e-14im; 1.5888896262186103e-11 - 1.0303480614499543e-13im -2.557126237504446e-12 + 1.759201163407723e-14im -8.475819811683215e-12 + 5.3531443609574795e-14im; -2.5781681021577177e-13 + 4.757590640631487e-15im 2.36818731889176e-12 - 4.569646499606389e-14im 1.116372367616482e-13 - 2.039935997276492e-15im; -1.0210743447568219e-11 - 5.905490368441375e-14im 1.3377918536056493e-12 + 7.190105205618706e-15im 2.5392856657302323e-11 + 1.5143842454586225e-13im; 2.4781693042536835e-11 - 1.6057018472176702e-13im -2.5342360504077476e-12 + 1.7306764301173096e-14im -8.40554044664581e-12 + 5.269404591748149e-14im; -2.348528974341763e-13 + 3.949450668269274e-15im 1.1449271118157543e-11 - 2.2093702114766968e-13im 1.0261871618968225e-13 - 1.7240213938923877e-15im; -1.0140560031409567e-11 - 5.828587508192886e-14im 1.3288225860409326e-12 + 7.0954601524623594e-15im 3.423954321087654e-11 + 2.0403371894291513e-13im],
    atol = 1e-6)

# output
true
```
"""
function hblinsolve(w, psc::CompiledCircuit,
        circuitdefs::Union{AbstractDict,AbstractVector}, signalfreq::Frequencies;
        sensitivitypairs = Tuple{String,Int,ComplexF64}[],
        sensitivityblockpairs = Tuple{String,Int,Any}[],
        symfreqvar = nothing, kwargs...)
    # the deprecated symbolic frequency variable, in circuit/legacy.jl
    if !isnothing(symfreqvar)
        circuitdefs isa AbstractDict || throw(ArgumentError(
            "symfreqvar needs the circuit definitions, not a resolved value table."))
        psc = frequencydependentcircuit(psc, circuitdefs, symfreqvar,
            :hblinsolve)
    end
    # the resolved value table, the sweep and the sensitivity pairs in their
    # canonical forms, so that the sweep below is compiled once for every
    # way of writing them
    return hblinsolve(sweepfrequencies(w), psc, resolvedvalues(psc, circuitdefs),
        signalfreq; sensitivitypairs = sensitivitypairtable(sensitivitypairs),
        sensitivityblockpairs = sensitivityblockpairtable(sensitivityblockpairs),
        kwargs...)
end

# the flat value table of a compiled circuit, resolved from the definitions
# or given already resolved
resolvedvalues(psc::CompiledCircuit, circuitdefs::AbstractDict) =
    componentvaluestonumber(psc.componentvalues, definitiontable(circuitdefs))
resolvedvalues(psc::CompiledCircuit, vvn::AbstractVector) = Vector{Any}(vvn)

function hblinsolve(w::Vector{Float64}, psc::CompiledCircuit,
    vvn::Vector{Any}, signalfreq::Frequencies;
    nonlinear = nothing,
    nbatches::Integer = Base.Threads.nthreads(),
    returnS::Bool = true, returnSnoise::Bool = false, returnQE::Bool = true,
    returnCM::Bool = true, returnnodeflux::Bool = false,
    returnnodefluxadjoint::Bool = false, returnvoltage::Bool = false,
    returnvoltageadjoint::Bool = false, keyedarrays::Bool = true,
    temperature = 0.0, returnCnoise::Bool = false,
    sensitivitynames::AbstractVector = String[],
    sensitivitypairs::Vector{Tuple{String,Int,ComplexF64}} =
        Tuple{String,Int,ComplexF64}[],
    sensitivityblockpairs::Vector{Tuple{String,Int,Any}} =
        Tuple{String,Int,Any}[],
    nsensitivityparameters::Integer = 0,
    sensitivitylabels::Union{Nothing,Vector{String}} = nothing,
    sensitivitynodeflux = nothing, sensitivityresidual = nothing,
    sensitivitymode::Symbol = :auto,
    returnSsensitivity::Bool = false, returnZ = nothing,
    returnZadjoint = nothing, returnZsensitivity = nothing,
    returnZsensitivityadjoint = nothing,
    factorization = nothing, backend = CPU(), debuglsys = false)


    # deprecation warnings for  `returnZ`, `returnZadjoint`,
    # `returnZsensitivity`, and `returnZsensitivityadjoint`
    if !isnothing(returnZ) || !isnothing(returnZadjoint) || 
        !isnothing(returnZsensitivity) || !isnothing(returnZsensitivityadjoint)
        Base.depwarn(lazy"The `returnZ`, `returnZadjoint`, `returnZsensitivity`, and `returnZsensitivityadjoint` kwargs have been removed. Please compute them from scattering parameters matrices.", :hblinsolve; force=true)
    end

    # the sweep in four stages, which read each other's fields by name:
    # the system and its noise channels, the sensitivity stamps, the sweep
    # over the frequencies, and the outputs. A stage is compiled again only
    # when what it depends on changes, and the operating point, the
    # sensitivity arrays and the requested outputs each reach one stage.
    wantsnoise = returnSnoise || returnQE || returnCM || returnCnoise
    s = linearizedsetup(w, psc, vvn, signalfreq, nonlinear,
        factorization, backend, temperature,
        String[String(n) for n in sensitivitynames],
        sensitivitypairs, sensitivityblockpairs; nbatches = nbatches,
        nsensitivityparameters = nsensitivityparameters,
        wantsnoise = wantsnoise)
    (; Nsignalmodes, signalnm, phimatrix, wpumpmodes, Nnodes,
        nodeindices, componenttypes, Nbranches,
        coupledbranches, Nauxmna, Nnodalmna, portindices, portnumbers,
        portimpedances, vvn, modes, Nlumpedpairs, sensitivitynames,
        sensitivityindices, stampgrouping, stampslots, Nports, bnm,
        noiseportimpedanceindices, ssys, lsys, factorization, refine,
        noiseplan, Nnoisechannels, channeltemperatures, channelsigns,
        noiseportimpedances, pumpfactorization, ondevice) = s
    sens = linearizedsensitivity(; psc, nonlinear, sensitivitynodeflux,
        sensitivityresidual, sensitivitypairs, sensitivityblockpairs,
        Nsignalmodes, signalnm, phimatrix, Nnodes, coupledbranches,
        Nnodalmna, Nlumpedpairs, sensitivitynames, sensitivityindices,
        stampgrouping, stampslots, ssys, lsys, pumpfactorization,
        sensitivitymode, returnSsensitivity)
    (; sensitivitystamps, sensitivityblockentries, sensitivitydAop,
        sensitivityreverse) = sens
    if debuglsys
        return (lsys=lsys, bnm=bnm, wpumpmodes=wpumpmodes,
            factorization=factorization,
            phimatrix=phimatrix, Ljb=signalnm.Ljb, Nnodes=Nnodes,
            Nmodes=Nsignalmodes, Nnodalmna=Nnodalmna, Nauxmna=Nauxmna,
            coupledbranches=coupledbranches, vvn=vvn,
            portindices=portindices, portimpedances=portimpedances,
            noiseportimpedanceindices=noiseportimpedanceindices,
            nodeindices=nodeindices, componenttypes=componenttypes,
            )
    end


    # The output arrays. An output which was not requested is a zero size
    # array, which signals that it is not to be computed. The sensitivities
    # are scaled by the input waves and depend on S itself, so S is computed
    # whenever they are requested even if it is not returned.
    outputarrays = linearizedsweep!(; w, backend, nbatches,
        nsensitivityparameters, sensitivitypairs, sensitivityblockpairs,
        Nsignalmodes, wpumpmodes, Nnodes, nodeindices, componenttypes,
        portindices, portimpedances, sensitivitynames, Nports, bnm,
        noiseportimpedanceindices, ssys, lsys, factorization, refine,
        noiseplan, Nnoisechannels, channeltemperatures, channelsigns,
        noiseportimpedances, sensitivitystamps, sensitivityblockentries,
        sensitivitydAop, sensitivityreverse, ondevice, returnS, returnSnoise,
        returnCnoise, returnSsensitivity, returnQE, returnCM, returnnodeflux,
        returnnodefluxadjoint, returnvoltage, returnvoltageadjoint)

    return linearizedoutputs(; psc, outputarrays, w, keyedarrays,
        sensitivitylabels, Nsignalmodes, Nnodes, nodeindices,
        componenttypes, Nbranches, portindices, portnumbers,
        portimpedances, modes, sensitivitynames, sensitivityindices, Nports,
        noiseportimpedanceindices, ssys, noiseplan, returnS, returnSnoise,
        returnCnoise, returnSsensitivity, returnQE, returnCM, returnnodeflux,
        returnnodefluxadjoint, returnvoltage, returnvoltageadjoint)
end

# the signal modes: the signal and the idlers offset from it by pump
# harmonics
signalfrequencies(Nmodulationharmonics; threewavemixing, fourwavemixing,
        maxintermodorder, maxharmonics) =
    truncfreqs(calcfreqsdft(Nmodulationharmonics); dc = true,
        odd = threewavemixing, even = fourwavemixing,
        maxintermodorder = maxintermodorder, maxharmonics = maxharmonics)

# the factorization of the pump Jacobian of the operating point
# sensitivities: the linearized solve's, except that a block factorization
# does not apply to it, since it has the real layout's blocks rather than
# the signal modes'
pumpjacobianfactorization(factorization) =
    factorization isa BlockFactorization ? KLUfactorization() : factorization

# the node voltages from the node fluxes: the voltage is the time
# derivative of the flux, which in the frequency domain is multiplication by
# im*w, for the rows `vv` has
function fluxtovoltage!(vv, phin, wmodes, Nmodes)
    @inbounds for t in axes(vv, 1)
        wm = im*wmodes[(t-1) % Nmodes + 1]
        for kc in axes(vv, 2)
            vv[t, kc] = wm*phin[t, kc]
        end
    end
    return vv
end

# the signal frequencies are finite and there is one at least, and the sweep
# has a batch to run in; checked before anything is built, and by `hbsolve`
# before the pump is solved
function checksweepinputs(w, nbatches)
    all(isfinite, w) || throw(ArgumentError("All signal frequencies must be finite."))
    isempty(w) && throw(ArgumentError("At least one signal frequency is required."))
    nbatches >= 1 || throw(ArgumentError(lazy"`nbatches` = $(nbatches) must be at least 1."))
    return nothing
end

"""
    linearizedsetup(w, psc, vvn, signalfreq, nonlinear,
        factorization, backend, temperature, sensitivitynames,
        sensitivitypairs, sensitivityblockpairs; nbatches,
        nsensitivityparameters, wantsnoise)

The first stage of [`hblinsolve`](@ref): the circuit matrices at the
signal mode count, the pump's cosine transform from the nonlinear solution
(or unity without one), the mode frequencies and the checks on them, the
modified nodal analysis padding, the sensitivity component indices and
their grouping, the port sources, the noise channels with their
temperatures, the [`HBLinearizedSystem`](@ref) with its system matrix
assembled, and the factorization the sweep uses. Returned as a named
tuple whose fields the later stages read by name.
"""
function linearizedsetup(w::Vector{Float64}, psc::CompiledCircuit,
    vvn::Vector{Any}, signalfreq::Frequencies,
    nonlinear, factorization, backend, temperature,
    sensitivitynames::Vector{String},
    sensitivitypairs::Vector{Tuple{String,Int,ComplexF64}},
    sensitivityblockpairs::Vector{Tuple{String,Int,Any}};
    nbatches::Integer, nsensitivityparameters::Integer, wantsnoise::Bool)
    checksweepinputs(w, nbatches)
    Nsignalmodes = length(signalfreq.modes)
    # the numeric matrices at the signal mode count, which differs from the
    # pump's
    topology = psc.topology
    signalnm = numericmatrices(psc, vvn; Nmodes = Nsignalmodes)

    if isnothing(nonlinear)

        allpumpfreq = calcfreqsrdft((0,))
        Amatrixmodes, Amatrixindices = hbmatind(allpumpfreq, signalfreq)
        Nwtuple = NTuple{length(allpumpfreq.Nw)+1,Int}((allpumpfreq.Nw..., length(signalnm.Ljb.nzval)))
        phimatrix = ones(Complex{Float64}, Nwtuple)

    else

        # the operating point has to be this circuit's: its junction
        # branches are the columns of the pump modulation read below
        (length(nonlinear.Ljb) == length(signalnm.Ljb) &&
            nonlinear.Ljb.nzind == signalnm.Ljb.nzind) || throw(ArgumentError(
            "the nonlinear solution is of another circuit: its branches or its junctions differ from this circuit's."))

        pumpfreq = nonlinear.frequencies
        # the signal modes are modulation harmonics of each pump tone
        length(pumpfreq.Nw) == length(signalfreq.Nw) || throw(ArgumentError(
            lazy"The signal modes have $(length(signalfreq.Nw)) tones but the nonlinear solution has $(length(pumpfreq.Nw)) pump frequencies; `Nmodulationharmonics` takes one count per pump tone."))

        allpumpfreq = calcfreqsrdft(pumpfreq.Nharmonics)
        # the maps between the pump's flux vector and its transform array,
        # which is all of its Fourier indices this reads
        vectomatmap, conjsourceindices, conjtargetindices =
            calcphiindices(pumpfreq, conjsym(pumpfreq))
        Npumpmodes = length(pumpfreq.modes)

        Amatrixmodes, Amatrixindices = hbmatind(allpumpfreq, signalfreq)

        # the frequency domain array of the Fourier transform, with one
        # column per junction
        Nwtuple = NTuple{length(pumpfreq.Nw)+1,Int}((pumpfreq.Nw..., length(nonlinear.Ljb.nzval)))

        phimatrix = zeros(Complex{Float64}, Nwtuple)

        # the time domain array and the transform plans
        phimatrixtd, irfftplan, rfftplan = plan_applynl(phimatrix, CPU())

        # arrange the pump branch fluxes for the inverse real transform,
        # conjugate modes included
        branchflux = nonlinear.Rbnm*nonlinear.nodeflux[:]
        phivectortomatrix!(
            branchflux[nonlinear.Ljbm.nzind], phimatrix,
            vectomatmap, conjsourceindices, conjtargetindices,
            length(nonlinear.Ljb.nzval)
        )

        # The Fourier coefficients of the derivative of the current-phase
        # relation at the pump, `cos(phi(t))` for the Josephson relation,
        # which is what modulates the linearized system.
        relations = calcjunctionrelations(psc.componenttypes, psc.nodeindices,
            psc.junctioncprs, topology.edge2indexdict, nonlinear.Ljb)
        if isnothing(relations)
            applynl!(
                phimatrix,
                phimatrixtd,
                cos,
                irfftplan,
                rfftplan,
            )
        else
            applyrelationnl!(phimatrix, phimatrixtd, similar(phimatrixtd),
                relations, relations.derivative, cos, irfftplan, rfftplan)
        end

    end

    # `wpumpmodes` is captured below; bind it to a concrete type
    wpumpmodes::Vector{Float64} = if isnothing(nonlinear)
        calcmodefreqs((0.0,),signalfreq.modes)
    else
        calcmodefreqs(nonlinear.w,signalfreq.modes)
    end

    # the fields of the nonlinear solution are untyped, so make the pump
    # frequencies a concrete tuple: everything below is per (frequency,
    # mode) and would box every value otherwise
    let wpumptuple = isnothing(nonlinear) ? (0.0,) :
            map(x -> Float64(real(x)), nonlinear.w)
        # the magnitudes of the pump terms of each mode do not depend on the
        # signal frequency, so their sum is precomputed per mode and only
        # abs(wi) joins per frequency
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

    # the first signal frequency, used for the setup below
    wmodes = w[1] .+ wpumpmodes

    # extract the elements we need
    Nnodes = psc.Nnodes
    nodeindices = psc.nodeindices
    componentnamedict = psc.componentnamedict
    componenttypes = psc.componenttypes
    Nbranches = topology.Nbranches
    edge2indexdict = topology.edge2indexdict
    Ljb = signalnm.Ljb
    Rbnm = signalnm.Rbnm
    Cnm = signalnm.Cnm
    Gnm = signalnm.Gnm
    invLnm = signalnm.invLnm

    # fail now if a component value still depends on a parameter which was
    # not defined (a frequency dependent value carries a closure and is fine)
    checkcomponentvaluesdefined(psc.componentnames, signalnm.vvn)

    # The modified nodal analysis augmentation: auxiliary branch current
    # variables for the mutually coupled inductor branches, whose inverse
    # inductance entries would otherwise diverge as the coupling coefficient
    # approaches one. No gauge fixing rows are needed, because a mode at
    # zero total frequency is not permitted.
    checkstaticstiffnessvalues(psc.componenttypes, signalnm.vvn)
    # the auxiliary unknowns: the branch currents of the mutually coupled
    # inductors and the port currents of the scattering blocks; resistors
    # are node conductances, as in hbnlsolve
    coupledbranches = mnacoupledbranches(signalnm.Mb)
    Nauxscattering = countscatteringports(psc)*Nsignalmodes
    Nauxmna = length(coupledbranches)*Nsignalmodes + Nauxscattering
    Nnodalmna = (psc.Nnodes-1)*Nsignalmodes
    # the frequency independent augmentation, filled below
    Amna0 = spzeros(Complex{Float64}, Nnodalmna + Nauxmna, Nnodalmna + Nauxmna)
    if !isempty(coupledbranches)
        # the coupled inductor branches, which `numericmatrices` excludes
        # from the inverse inductance matrix: their constitutive equations
        # and Kirchhoff current law couplings, with unscaled branch currents
        # as the auxiliary variables (Lscale = 1) to match the unscaled
        # matrices of this solver
        AmnaL = calcAmnaind(coupledbranches, signalnm.Lb, signalnm.Mb,
            topology.Rbn, Nsignalmodes, Nnodalmna, Nnodalmna + Nauxmna, 1)
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
    portimpedances = signalnm.portimpedances
    vvn = signalnm.vvn
    modes = signalfreq.modes

    # Design parameter sensitivities: the caller names physical parameters
    # rather than components, through pairs (componentname, parameterindex,
    # alpha) with alpha = (dv/dp)/v the direction of the component value
    # under the parameter. The stamps of one parameter merge into one
    # contraction.
    # A component is named by its identifier or given by its flat index; a
    # scattering block has no table entry and is named by its instance
    # path, so a block pair resolves to the block's ordinal instead, which
    # only the block stamps read. The names are kept beside the indices to
    # label the sensitivity axis.
    Nlumpedpairs = length(sensitivitypairs)
    sensitivityindices = Int[]
    if !isempty(sensitivitypairs) || !isempty(sensitivityblockpairs)
        isempty(sensitivitynames) || throw(ArgumentError(
            "pass either sensitivitynames or sensitivitypairs, not both"))
        nsensitivityparameters > 0 || throw(ArgumentError(
            "sensitivitypairs requires nsensitivityparameters"))
        for t in sensitivitypairs
            push!(sensitivityindices, componentindex(psc, t[1]))
        end
        for t in sensitivityblockpairs
            b = scatteringblockindex(psc, t[1])
            iszero(b) && throw(ArgumentError(lazy"The block pair component $(t[1]) is not a scattering block of this circuit."))
            # the derivative is stamped in the block's place, so it is of
            # the block's kind
            target = psc.scatteringblocks[b].definition
            (t[3] isa LinearizedScattering) == (target isa LinearizedScattering) ||
                throw(ArgumentError(lazy"The derivative in the block pair of $(t[1]) is a $(nameof(typeof(t[3]))) but the block is a $(nameof(typeof(target))); a derivative must be of its block's kind."))
            push!(sensitivityindices, b)
        end
        sensitivitynames = vcat(
            String[psc.componentnames[i]
                for i in view(sensitivityindices, 1:Nlumpedpairs)],
            String[psc.scatteringblocks[i].path
                for i in view(sensitivityindices,
                    (Nlumpedpairs+1):length(sensitivityindices))])
    else
        for name in sensitivitynames
            idx = get(componentnamedict, String(name), 0)
            if iszero(idx)
                iszero(scatteringblockindex(psc, name)) &&
                    throw(ArgumentError(lazy"The component $(name) is not in this circuit."))
                throw(ArgumentError(lazy"Sensitivities with respect to scattering blocks require a block derivative; got $(name) as a plain component."))
            end
            push!(sensitivityindices, idx)
        end
    end

    # The grouping of the pairs into merged stamps and the design parameter
    # slot each group accumulates into, computed before the residual
    # derivative columns and the operating point solves because both merge
    # with the same grouping. A scattering block pair is its own group: its
    # stamp is rebuilt per frequency and cannot be concatenated.
    stampgrouping, stampslots = if isempty(sensitivitypairs) &&
            isempty(sensitivityblockpairs)
        Vector{Int}[], Int[]
    else
        g, sl = parametergrouping(sensitivitypairs,
            view(sensitivityindices, 1:Nlumpedpairs), componenttypes,
            Dict(idx => p
                for (p, idx) in enumerate(signalnm.portenvironmentindices)
                if !iszero(idx)))
        # the block pairs of one parameter form one group: their stamps
        # reach the entries of different instances, so they concatenate
        blockgroup = Dict{Int,Int}()
        for (bi, bp) in enumerate(sensitivityblockpairs)
            parameter = Int(bp[2])
            k = get!(blockgroup, parameter) do
                push!(g, Int[]); push!(sl, parameter)
                length(g)
            end
            push!(g[k], Nlumpedpairs + bi)
        end
        g, sl
    end

    Nports = length(portindices)

    # the source terms in the branch basis: a unit current source at each
    # port and mode
    bbm = zeros(Complex{Float64},Nbranches*Nsignalmodes,Nsignalmodes*Nports)

    for (i,val) in enumerate(portindices)
        key = (nodeindices[1,val],nodeindices[2,val])
        for j = 1:Nsignalmodes
            bbm[(edge2indexdict[key]-1)*Nsignalmodes+j,(i-1)*Nsignalmodes+j] = 1
        end
    end

    # and in the node basis
    bnm = transpose(Rbnm)*bbm
    if Nauxmna > 0
        bnm = vcat(bnm, zeros(eltype(bnm), Nauxmna, size(bnm, 2)))
    end
    # A lossy capacitor or inductor is a noise channel when its value has a
    # nonzero imaginary part. A frequency dependent value can have one at
    # some modes of the sweep and not at others, so it is a channel when it
    # has one at any of them, and where it has none its channel carries no
    # noise (see `noisewavescale`). Which resistor is a port's termination
    # is a role, not a value, so that part does not depend on the values.
    noiseportimpedanceindices = if any(checkissymbolic, vvn)
        lossy(v) = !iszero(imag(v))
        lossyinsweep(v) = any(wi -> any(wm -> lossy(substitutefreq(v,
            wi + wm)), wpumpmodes), w)
        [i for i in noisecandidates(psc) if psc.componenttypes[i] === :R ||
            (checkissymbolic(vvn[i]) ? lossyinsweep(vvn[i]) :
                vvn[i] isa Complex && lossy(vvn[i]))]
    else
        signalnm.noiseportimpedanceindices
    end

    # whether any entry holds a symbolic value, in which case the stored
    # values themselves change with the frequency and the sweep stays on
    # the host
    symbolicvalues = !isempty(symbolicindices(Cnm)) ||
        !isempty(symbolicindices(Gnm)) || !isempty(symbolicindices(invLnm))


    Cnmcopy = freqsubst(Cnm,wmodes)
    Gnmcopy = freqsubst(Gnm,wmodes)
    invLnmcopy = freqsubst(invLnm,wmodes)

    # The linearized system object: the sparsity structure of the system
    # matrix Asparse = AoLjnm + invLnm + Gnm + Cnm with stored zeros, a
    # plan scattering the pump modulation term AoLjnm = Rbnm'*AoLjbm*Rbnm
    # into it from the Fourier coefficients of cos(phi(t)), index maps for
    # the frequency dependent linear terms, and the pump modulation term
    # with its conjugate. This shares the machinery of the nonlinear
    # Jacobians (see `HBLinearizedSystem`); the per frequency matrices are
    # assembled from it by `assemblesystemmatrix!`. The scattering blocks
    # contribute constant Kirchhoff current law couplings of their auxiliary
    # port currents, folded into the augmentation, and frequency dependent
    # constitutive equations assembled per frequency.
    ssys = scatteringstampsystem(psc.scatteringblocks, Nsignalmodes;
        auxoffset = Nnodalmna + Nauxmna - Nauxscattering,
        Ntotal = Nnodalmna + Nauxmna, scale = 1.0,
        modeoffsets = wpumpmodes,
        iscale = auxcurrentscale(calcsolverscale(w, psc.componenttypes,
            signalnm.vvn, signalnm.portimpedances, signalnm.Lmean)))
    if !isnothing(ssys)
        Amna0 = spaddkeepzeros(Amna0, ssys.kcl)
    end
    lsys = HBLinearizedSystem(Amatrixindices, signalnm.Ljb, Rbnmmna,
        Nsignalmodes, topology.Nbranches, phimatrix, invLnmcopy, Gnmcopy, Cnmcopy,
        invLnm, Gnm, Cnm, symbolicvalues, Amna0, wpumpmodes,
        Nnodes; scattering = ssys)
    Asparse = lsys.Asparse

    # The sweep runs on the backend when it is a device, no stored value
    # depends on the frequency and no block sensitivity is asked for, whose
    # stamps are rebuilt per frequency on the host; otherwise on host
    # threads. The factorization when none was given is chosen for where
    # the sweep runs, by the number of tones and the memory (see
    # `linearizedfactorization`).
    ondevice = !(backend isa CPU) && cansweepondevice(lsys) &&
        isempty(sensitivityblockpairs)
    if isnothing(factorization)
        factorization = linearizedfactorization(Asparse, Nsignalmodes,
            length(first(signalfreq.modes)), ondevice ? backend : CPU();
            nbatches = nbatches)
    elseif factorization isa CUDSSFactorization && !ondevice
        throw(ArgumentError("the sweep of this circuit runs on the host, because a component value depends on the frequency or a block sensitivity is asked for, and CUDSSFactorization factorizes on a device; leave `factorization` to its default or pass a host factorization such as KLUfactorization()."))
    end
    # whether single precision block factors refine to double: the
    # factorization's choice, and moot for any other factorization
    refine = !(factorization isa BlockFactorization) || factorization.refine

    # the vacuum noise channels of the dissipative scattering blocks, which
    # follow the lumped noise channels in the rows of the noise scattering
    # matrix; a block declared lossless has none. A pumped block's
    # declaration, its losslessness or its stated covariance, is checked
    # over the modes of every solve whatever outputs are asked for, so
    # that the block is accepted or refused on its declaration and not
    # on the outputs
    checklosslessblocks(ssys, w, wpumpmodes)
    checkpumpedblockmodels(ssys, w, wpumpmodes)
    # every output derived from the noise waves needs the channels planned,
    # the covariance included, and the noise a block states admitted at
    # every frequency they are formed at
    noiseplan = if wantsnoise
        checkblocknoisemodels(ssys, w, wpumpmodes)
        planscatteringnoise(ssys)
    else
        nothing
    end
    Nnoisechannels = length(noiseportimpedanceindices) +
        (isnothing(noiseplan) ? 0 : noiseplan.Nchannels)
    # the temperature of each noise channel: the analysis default, or the
    # one a component or block states for itself
    channeltemperatures = noisechanneltemperatures(psc,
        noiseportimpedanceindices, noiseplan, ssys, temperature)
    # the sign of each channel in the commutation relations, which only
    # the channels of the conjugate kind of a block which states its noise
    # reverse
    channelsigns = noisechannelsigns(noiseportimpedanceindices, noiseplan,
        ssys)

    noiseportimpedances = [vvn[i] for i in noiseportimpedanceindices]

    # assemble the system matrix at the first frequency for the symbolic
    # analysis of the factorization
    assemblesystemmatrix!(Asparse, lsys, wmodes)

    pumpfactorization = pumpjacobianfactorization(factorization)
    return (; Nsignalmodes, signalnm, phimatrix, wpumpmodes, Nnodes, nodeindices, componenttypes, Nbranches, coupledbranches, Nauxmna, Nnodalmna, portindices, portnumbers, portimpedances, vvn, modes, Nlumpedpairs, sensitivitynames, sensitivityindices, stampgrouping, stampslots, Nports, bnm, noiseportimpedanceindices, ssys, lsys, factorization, refine, noiseplan, Nnoisechannels, channeltemperatures, channelsigns, noiseportimpedances, pumpfactorization, ondevice)
end

"""
    linearizedsensitivity(; ...)

The second stage of [`hblinsolve`](@ref): the sensitivity stamps of the
requested components and scattering blocks, and the operating point
contribution in the contraction order chosen (`sensitivitymode`), either
as forward stamps or as a [`ReverseSensitivity`](@ref). Empty when no
sensitivity was asked for. The keywords are the fields of
[`linearizedsetup`](@ref) this stage reads, the operating point and the
sensitivity arrays.
"""
function linearizedsensitivity(;
        psc, nonlinear, sensitivitynodeflux, sensitivityresidual,
        sensitivitypairs, sensitivityblockpairs, Nsignalmodes, signalnm,
        phimatrix, Nnodes, coupledbranches, Nnodalmna, Nlumpedpairs,
        sensitivitynames, sensitivityindices, stampgrouping, stampslots,
        ssys, lsys, pumpfactorization, sensitivitymode,
        returnSsensitivity)
    sensitivitystamps, sensitivityblockentries = if returnSsensitivity
        st = calcsensitivitystamps(
            sensitivityindices[1:(isempty(sensitivityblockpairs) ?
                end : Nlumpedpairs)],
            psc, signalnm, lsys, phimatrix, coupledbranches, Nnodalmna,
            Nsignalmodes, Nnodes)
        if !isempty(sensitivitypairs)
            st = [reparameterize(st[k], sensitivitypairs[k][3],
                sensitivitypairs[k][2]) for k in eachindex(st)]
        end
        entries = Tuple{Int,Any,Any,Vector{Int}}[]
        if !isempty(sensitivityblockpairs)
            isnothing(ssys) && throw(ArgumentError(
                "the circuit has no scattering blocks to take a block sensitivity of"))
            for (bi, bp) in enumerate(sensitivityblockpairs)
                # the block ordinal held by the block tail of
                # `sensitivityindices`
                dsys, zsys, position, rows, cols = targetstampsystems(ssys,
                    sensitivityindices[Nlumpedpairs + bi], bp[3])
                push!(st, blocksensitivitystamp(rows, cols, Int(bp[2])))
                push!(entries, (0, dsys, zsys, position))
            end
        end
        if !isempty(stampgrouping)
            lengths = [length(s.vals) for s in st]
            st = mergestamps(st, stampgrouping)
            # each block pair's place in the merged vector: its group's
            # stamp, from the offset of its own entries in the concatenation
            for (gi, g) in enumerate(stampgrouping)
                offset = 0
                for i in g
                    if i > Nlumpedpairs
                        bi = i - Nlumpedpairs
                        entries[bi] = (gi, entries[bi][2], entries[bi][3],
                            entries[bi][4] .+ offset)
                    end
                    offset += lengths[i]
                end
            end
        end
        st, entries
    else
        SensitivityStamp[], Tuple{Int,Any,Any,Vector{Int}}[]
    end

    # The contribution of the shift of the pump operating point, when its
    # derivatives are supplied either directly (`sensitivitynodeflux`) or as
    # residual derivatives (`sensitivityresidual`), from which the former
    # are computed only if the forward contraction order needs them. The
    # forward order costs one product against the sparsity structure of the
    # system matrix per component and frequency; the reverse order costs
    # transposed solves per output port mode pair and a sparse inner product
    # per component, so it wins once there are more components than output
    # port mode pairs.
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
    # Validate the operating point inputs here: the contractions index
    # arrays sized from them inside @inbounds loops, so malformed inputs
    # must be rejected rather than discovered as memory corruption. The two
    # contraction orders read different inputs, so supplying both would let
    # inconsistent inputs give order dependent results; reject that too.
    if useoperatingpoint
        op = nonlinear.operatingpoint
        if !isnothing(sensitivitynodeflux) && !isnothing(sensitivityresidual)
            throw(ArgumentError("Supply either sensitivitynodeflux or sensitivityresidual, not both: the forward contraction order reads the former and the reverse order the latter, so inconsistent inputs would make the result depend on the contraction order."))
        end
        if !isnothing(sensitivitynodeflux) &&
                size(sensitivitynodeflux) != (length(op.x), length(sensitivitynames))
            throw(DimensionMismatch(lazy"sensitivitynodeflux must be (operating point state length) x (number of sensitivity components) = ($(length(op.x)), $(length(sensitivitynames))), got $(size(sensitivitynodeflux))."))
        end
        # the residual derivatives are sized by the system the implicit
        # function theorem is applied to: the canonical one when a direct
        # current block is active, the harmonic one otherwise
        if !isnothing(sensitivityresidual) &&
                size(sensitivityresidual) != (sensitivitydim(op), length(sensitivitynames))
            throw(DimensionMismatch(lazy"sensitivityresidual must be (rows of the Jacobian the sensitivity is taken through) x (number of sensitivity components) = ($(sensitivitydim(op)), $(length(sensitivitynames))), got $(size(sensitivityresidual))."))
        end
        if !isempty(sensitivitypairs) && !isnothing(sensitivitynodeflux)
            throw(ArgumentError("sensitivitynodeflux is per component and cannot be combined with sensitivitypairs; supply sensitivityresidual (or neither) instead."))
        end
    end
    # The residual derivative columns are per pair and merge with the same
    # grouping as the stamps before anything is sized or solved from them;
    # summing is right because the residual is linear in each component's
    # contribution.
    if !isnothing(sensitivityresidual) && !isempty(stampgrouping) &&
            length(stampgrouping) != size(sensitivityresidual, 2)
        sensitivityresidual = hcat(
            [sum(sensitivityresidual[:, g], dims = 2)
             for g in stampgrouping]...)
    end
    # Without Josephson junctions the system matrix does not depend on the
    # operating point, so its contribution is zero: drop the operating point
    # inputs rather than build junction shaped contractions on an empty
    # junction set.
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
        # The forward order costs one product per component; the reverse
        # order costs one transposed solve per output port mode pair whatever
        # the component count. A solve costs several times a product, so the
        # crossover is at a multiple of the pair count; the factor of eight
        # is an empirical estimate, and near the crossover either choice
        # costs about the same.
        Ncomponents = isnothing(sensitivityresidual) ?
            size(sensitivitynodeflux, 2) : size(sensitivityresidual, 2)
        !isnothing(sensitivityresidual) &&
            isempty(sensitivityblockpairs) && Ncomponents >
            8*(length(signalnm.portindices)*Nsignalmodes)^2
    end

    sensitivitydAop = if useoperatingpoint && !usereverse
        # the forward order contracts against the operating point shift
        # itself, so compute it from the residual derivatives if it was not
        # supplied; the reverse order works from the residual derivatives
        # directly and skips these per component solves
        dx = if isnothing(sensitivitynodeflux)
            calcnodefluxsensitivity(nonlinear.operatingpoint,
                sensitivityresidual; factorization = pumpfactorization)
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


    return (; sensitivitystamps, sensitivityblockentries, sensitivitydAop, sensitivityreverse)
end

# The number of BLAS threads is a setting of the whole process, so the
# sweeps which run with one are counted: the first to start saves the
# setting and the last to finish restores it, which is right however many
# run at once.
const BLASLIMITED = Ref(0)
const BLASSAVED = Ref(1)
const BLASLIMITLOCK = ReentrantLock()

function withoneblasthread(f)
    lock(BLASLIMITLOCK) do
        if BLASLIMITED[] == 0
            BLASSAVED[] = BLAS.get_num_threads()
            BLAS.set_num_threads(1)
        end
        BLASLIMITED[] += 1
    end
    try
        return f()
    finally
        lock(BLASLIMITLOCK) do
            BLASLIMITED[] -= 1
            BLASLIMITED[] == 0 && BLAS.set_num_threads(BLASSAVED[])
        end
    end
end

"""
    linearizedsweep!(; ...)

The third stage of [`hblinsolve`](@ref): the sweep over the signal
frequencies, on the backend in batches when the system allows it and on
host threads otherwise, each frequency solved by [`hblinsolve_inner!`](@ref)
into the [`LinearizedArrays`](@ref) returned. The keywords are the fields
of the two stages before which this one reads, and the requested outputs.
"""
function linearizedsweep!(;
        w, backend, nbatches, nsensitivityparameters,
        sensitivitypairs, sensitivityblockpairs, Nsignalmodes, wpumpmodes,
        Nnodes, nodeindices, componenttypes, portindices, portimpedances,
        sensitivitynames, Nports, bnm, noiseportimpedanceindices, ssys, lsys,
        factorization, refine, noiseplan, Nnoisechannels,
        channeltemperatures, channelsigns, noiseportimpedances,
        sensitivitystamps, sensitivityblockentries, sensitivitydAop,
        sensitivityreverse, ondevice, returnS, returnSnoise, returnCnoise,
        returnSsensitivity, returnQE, returnCM, returnnodeflux,
        returnnodefluxadjoint, returnvoltage, returnvoltageadjoint)
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
        # with pairs of either kind the output axis is the design
        # parameters, whose count the pairs' slots may exceed the pair count
        Ncomponents = isempty(sensitivitypairs) &&
            isempty(sensitivityblockpairs) ?
            length(sensitivitynames) : nsensitivityparameters,
        Nnodes = Nnodes,
        Nfrequencies = length(w))

    # Solve the linear system at each frequency. The frequencies are
    # independent, so they are split into `nbatches` batches solved by
    # tasks in parallel, each batch reusing one workspace and one symbolic
    # factorization. On a device backend the solve of a batch is done there
    # instead: the matrices of a batch share one sparsity pattern, so they
    # are assembled by one kernel and factorized and solved as a uniform
    # batch, while the frequency loop and every output it computes stay the
    # host ones. The device path is not used when a component value is
    # frequency dependent, since the assembly is then not a constant
    # quadratic in the frequency.
    sensitivitytuple = (stamps = sensitivitystamps, dAop = sensitivitydAop,
        reverse = sensitivityreverse,
        blockentries = sensitivityblockentries)
    if ondevice
        # what the solutions are read for decides whether the whole solution
        # comes back from the device or only some rows: the scattering
        # parameters read the port rows of the forward solution and the
        # noise port rows of the adjoint one, while the node flux, voltage
        # and sensitivity outputs read all of both
        fullforward = !isempty(outputarrays.nodeflux) ||
            !isempty(outputarrays.voltage) ||
            !isempty(outputarrays.Ssensitivity)
        fulladjoint = !isempty(outputarrays.nodefluxadjoint) ||
            !isempty(outputarrays.voltageadjoint) ||
            !isempty(outputarrays.Ssensitivity)
        # The noise scattering parameters are computed where the adjoint
        # solution is, so when they are its only reader the adjoint solution
        # is never copied back. The noise channels of the scattering blocks
        # are formed on the device when the blocks' scattering parameters
        # can be evaluated there (the same condition as for their stamps);
        # otherwise the whole adjoint solution comes back and every channel
        # is formed on the host.
        # a block which states its noise has channels of both kinds, which
        # the kernels do not form, so its circuit's channels stay on the host
        blocknoiseondevice = isnothing(noiseplan) ||
            (candeviceevaluate(lsys.scattering) && isnothing(channelsigns))
        wantsnoise = !isempty(outputarrays.Snoise) ||
            !isempty(outputarrays.QE) || !isempty(outputarrays.CM) ||
            !isempty(outputarrays.Cnoise)
        if wantsnoise && !blocknoiseondevice
            if isnothing(channelsigns)
                @warn "A scattering block of this circuit has scattering parameters which cannot be evaluated on the backend, so its noise channels are formed on the host and the whole adjoint solution is copied back at every frequency, which is the largest transfer in the sweep. Give the block's data as a constant matrix, as tabulated data, or as a callable of one entry at a time (`ScatteringParameters(f; nports, form = :entry)`) to keep the noise on the backend." maxlog=1
            else
                @warn "A scattering block of this circuit states its noise with a NoiseCovariance, whose channels are formed on the host, so the whole adjoint solution is copied back at every frequency, which is the largest transfer in the sweep." maxlog=1
            end
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
                portindices, Nsignalmodes)),
            adjointspec; factorization = factorization, refine = refine)
        # The device solves a batch of frequencies, then the host computes
        # their outputs. Each worker's workspace holds a solution buffer the
        # size of the whole circuit's state, so the number of workers is
        # chosen by how much host work a frequency carries: without the
        # sensitivities the outputs are cheap and one worker wins, while the
        # sensitivities cost a product per component and frequency and are
        # worth spreading across workers.
        nworkers = isempty(outputarrays.Ssensitivity) ? 1 :
            max(1, min(nbatches, Base.Threads.nthreads()))
        wss = [LinearizedWorkspace(outputarrays, sensitivitytuple, lsys,
            Nports, Nsignalmodes, Nnoisechannels,
            length(wpumpmodes), factorization; assembles = false)
            for _ in 1:nworkers]
        # the noise scattering parameters are reduced on the device against
        # scratch of their own, so each worker computing them needs its own
        noisecbs = isnothing(devicenoiseplan) ? nothing :
            [devicenoise(devicenoiseplan,
                (isnothing(deviceblocknoiseplan) || t == 1) ?
                    deviceblocknoiseplan : withfactors(deviceblocknoiseplan),
                solutions.providers, i -> adjointdevice(solutions, i),
                size(bnm, 2), wpumpmodes, w, !isempty(outputarrays.Snoise),
                channeltemperatures)
             for t in 1:nworkers]
        nb = solutions.nb
        inner!(t, batch) = hblinsolve_inner!(wss[t], outputarrays,
            sensitivitytuple, lsys, bnm,
            portindices, noiseportimpedanceindices,
            portimpedances, noiseportimpedances, nodeindices,
            componenttypes, w, wpumpmodes, Nsignalmodes, Nnodes,
            batch, factorization;
            presolved = (i, phin) -> forwardsolution!(phin, solutions, i),
            presolvedadjoint = isnothing(solutions.adj) ? nothing :
                (i, phin) -> adjointsolution!(phin, solutions, i),
            presolvednoise = isnothing(noisecbs) ? nothing : noisecbs[t],
            noiseplan = noiseplan,
            channeltemperatures = channeltemperatures,
            channelsigns = channelsigns)
        for lo in 1:nb:length(w)
            hi = min(lo + nb - 1, length(w))
            # only this touches the device, serially; the outputs of the
            # batch it staged are host work on disjoint frequencies
            solvebatch!(solutions, lo)
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
        # every worker factorizes the same pattern, and with the reverse
        # sensitivities the same pump Jacobian, so their fill reducing
        # orderings are chosen once here for every worker's first
        # factorization
        ordering = fillordering(factorization, lsys.Asparse)
        pumpordering = isnothing(sensitivityreverse) ? nothing :
            fillordering(pumpjacobianfactorization(factorization),
                sensitivityjacobian(sensitivityreverse.op))
        runbatches() = Threads.@sync for batch in batches
            Base.Threads.@spawn hblinsolve_inner!(
                LinearizedWorkspace(outputarrays, sensitivitytuple, lsys,
                    Nports, Nsignalmodes, Nnoisechannels,
                    length(wpumpmodes), factorization;
                    ordering = ordering, pumpordering = pumpordering),
                outputarrays, sensitivitytuple,
                lsys, bnm,
                portindices, noiseportimpedanceindices,
                portimpedances, noiseportimpedances, nodeindices, componenttypes,
                w, wpumpmodes, Nsignalmodes, Nnodes, batch,
                factorization; noiseplan = noiseplan,
                channeltemperatures = channeltemperatures,
                channelsigns = channelsigns, refine = refine)
        end
        # a block factorization's dense products would each spawn BLAS
        # threads; with several batches the batches are the parallelism, so
        # each task gets one BLAS thread while the sweep runs
        if factorization isa BlockFactorization && length(batches) > 1
            withoneblasthread(runbatches)
        else
            runbatches()
        end
    end

    return outputarrays
end

"""
    linearizedoutputs(; ...)

The last stage of [`hblinsolve`](@ref): the [`LinearizedHB`](@ref) of the
sweep, its arrays keyed by mode, port, node and frequency when asked. It
takes the compiled circuit for the names it records, which no stage
before it reads.
"""
function linearizedoutputs(; psc,
        outputarrays, w, keyedarrays, sensitivitylabels, Nsignalmodes,
        Nnodes, nodeindices,
        componenttypes, Nbranches, portindices,
        portnumbers, portimpedances, modes, sensitivitynames,
        sensitivityindices, Nports, noiseportimpedanceindices, ssys,
        noiseplan, returnS, returnSnoise, returnCnoise, returnSsensitivity,
        returnQE, returnCM, returnnodeflux, returnnodefluxadjoint,
        returnvoltage, returnvoltageadjoint)
    # the names the result records are the compiled circuit's own
    nodenames = psc.nodenames
    componentnames = psc.componentnames
    componentnamedict = psc.componentnamedict
    mutualinductorbranchnames = coupledinductornames(psc)
    signalindex = 1
    # the quantum efficiency of an ideal two mode amplifier with the same
    # gain
    QEideal = outputarrays.QEideal

    # the requested outputs as keyed arrays when `keyedarrays = true`, and
    # as they are otherwise: `keyed(f, requested, a)` applies the keying `f`
    keyed(f, requested, a) = requested && keyedarrays ? f(a) : a
    byports(a) = Stokeyed(a, modes, portnumbers, modes, portnumbers, w)
    bynodes(a) = nodevariabletokeyed(a, modes, nodenames, modes, portnumbers,
        w)
    # S is computed whenever the sensitivities are, which scale by it, and
    # is returned only when asked for
    Sout = returnS ? keyed(byports, true, outputarrays.S) :
        zeros(Complex{Float64}, 0, 0, 0)
    Snoiseout = keyed(returnSnoise, outputarrays.Snoise) do a
        Snoisetokeyed(a, modes, noisechannelnames(componentnames,
            noiseportimpedanceindices, noiseplan, ssys), modes, portnumbers, w)
    end
    # the added noise covariance is indexed by output port mode on both
    # sides, like the scattering matrix
    Cnoiseout = keyed(byports, returnCnoise, outputarrays.Cnoise)
    Ssensitivityout = keyed(returnSsensitivity, outputarrays.Ssensitivity) do a
        Ssensitivitytokeyed(a, modes, portnumbers, modes, portnumbers,
            isnothing(sensitivitylabels) ? sensitivitynames :
                sensitivitylabels, w)
    end
    QEout = keyed(byports, returnQE, outputarrays.QE)
    QEidealout = keyed(byports, returnQE, QEideal)
    CMout = keyed(a -> CMtokeyed(a, modes, portnumbers, w), returnCM,
        outputarrays.CM)
    nodefluxout = keyed(bynodes, returnnodeflux, outputarrays.nodeflux)
    nodefluxadjointout = keyed(bynodes, returnnodefluxadjoint,
        outputarrays.nodefluxadjoint)
    voltageout = keyed(bynodes, returnvoltage, outputarrays.voltage)
    voltageadjointout = keyed(bynodes, returnvoltageadjoint,
        outputarrays.voltageadjoint)

    return LinearizedHB(w, modes, Sout, Snoiseout, Cnoiseout, Ssensitivityout,
        QEout,
        QEidealout, CMout, nodefluxout, nodefluxadjointout, voltageout,
        voltageadjointout, nodenames, nodeindices, componentnames,
        componenttypes, componentnamedict, mutualinductorbranchnames,
        portnumbers, portindices, portimpedances,
        noiseportimpedanceindices, sensitivitynames, sensitivityindices,
        Nsignalmodes, Nnodes, Nbranches, Nports, signalindex)
end

"""
    LinearizedArrays(; requestS, requestSnoise, requestCnoise,
        requestSsensitivity, requestQE, requestCM, requestnodeflux,
        requestnodefluxadjoint, requestvoltage, requestvoltageadjoint,
        Nports, Nmodes,
        Nnoisechannels, Ncomponents, Nnodes, Nfrequencies)

The preallocated output arrays of one [`hblinsolve`](@ref) run, filled per
frequency by [`hblinsolve_inner!`](@ref). An output which was not requested
is a zero size array of the same dimensionality, which signals, through
`isempty`, that it is not to be computed. The output shapes and the request
conditions are defined here and nowhere else.
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

The largest buffers are the size of the solution and of the system matrix,
so a worker is given one workspace once and reuses it across every range of
frequencies it is handed. This is what lets the device path hand out one
batch of frequencies at a time while keeping the host work parallel.
"""
struct LinearizedWorkspace{TA,TC,TR}
    phin::Matrix{Complex{Float64}}
    inputwave::Vector{Complex{Float64}}
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
    # the scratch of the scattering block vacuum noise channels, and of
    # the block evaluation at each frequency of the assembly
    scatteringnoisework::ScatteringNoiseWorkspace
    scatteringwork::ScatteringWorkspace
    # the occupation of each noise channel mode, rebuilt per frequency
    occupation::Vector{Float64}
    # the reduction of the noise scattering matrix the quantum efficiency
    # and the commutation relations read, formed per frequency, and the
    # scratch of their compensated row sums
    noise::NoiseReduction{Vector{Float64}}
    rowsum::Vector{Float64}
    rowcomp::Vector{Float64}
    # this worker's private scattering block sensitivity state, or nothing
    blocksens::Any
end

"""
    LinearizedWorkspace(arrays::LinearizedArrays, sensitivity, lsys, Nports,
        Nmodes, Nnoisechannels, Nwpumpmodes, factorization;
        assembles::Bool = true, ordering = nothing, pumpordering = nothing)

Make the scratch of one worker. `assembles` is false when the solutions are
supplied from elsewhere and nothing writes into a copy of the system matrix,
which on a large problem is the biggest allocation here. `ordering` and
`pumpordering`, when given, are the fill reducing orderings of the system
matrix and of the pump Jacobian (see [`fillordering`](@ref)), which the
worker's first factorizations then take rather than choose.
"""
function LinearizedWorkspace(arrays::LinearizedArrays, sensitivity, lsys,
    Nports::Integer, Nmodes::Integer, Nnoisechannels::Integer,
    Nwpumpmodes::Integer, factorization; assembles::Bool = true,
    ordering = nothing, pumpordering = nothing)

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
    # one factorization of the pump Jacobian per worker; a sparse
    # factorization cannot be solved against from several threads at once
    sensitivitycache = if isnothing(sensitivity.reverse)
        nothing
    else
        c = FactorizationCache()
        # the canonical Jacobian when a direct current block is active: the
        # adjoint has to be taken through the system which was solved
        J = sensitivityjacobian(sensitivity.reverse.op)
        isnothing(pumpordering) || seedordering!(c, J, pumpordering)
        tryfactorize!(c, pumpjacobianfactorization(factorization), J)
        c
    end
    # a copy of the system matrix, because it is modified per frequency,
    # potentially by several workers at once
    A = assembles ? copy(lsys.Asparse) : lsys.Asparse
    cache = FactorizationCache()
    isnothing(ordering) || seedordering!(cache, A, ordering)
    return LinearizedWorkspace(
        cplx(n, np), zeros(Complex{Float64}, np), cplx(np, np),
        cplx(Nnoisechannels*Nmodes, np),
        phinforward,
        wantsop ? copy(lsys.Asparse) : lsys.Asparse,
        wantsop ? cplx(size(phinforward)...) : cplx(0, 0),
        cplx(np, np), sensitivitycache,
        isnothing(sensitivity.reverse) ? nothing :
            ReverseSensitivityBuffers(sensitivity.reverse, np),
        zeros(Complex{Float64}, np), zeros(Complex{Float64}, np),
        zeros(Float64, Nwpumpmodes),
        A, cache, cplx(np, np), cplx(Nnoisechannels*Nmodes, np),
        ScatteringNoiseWorkspace(), ScatteringWorkspace(),
        zeros(Float64, Nnoisechannels*Nmodes),
        NoiseReduction(zeros(Float64, np), zeros(Float64, np)),
        zeros(Float64, np), zeros(Float64, np),
        isempty(sensitivity.blockentries) ? nothing :
            WorkerBlockSensitivity(sensitivity.stamps,
                sensitivity.blockentries))
end

"""
    hblinsolve_inner!(ws::LinearizedWorkspace, arrays::LinearizedArrays,
        sensitivity, lsys, bnm, portindices, noiseportimpedanceindices,
        portimpedances, noiseportimpedances, nodeindices, componenttypes,
        w, wpumpmodes, Nmodes, Nnodes, wi, factorization;
        noiseplan = nothing, channeltemperatures = nothing,
        channelsigns = nothing, presolved = nothing,
        presolvedadjoint = nothing, presolvednoise = nothing, refine = true)

Solve the linearized problem at the frequencies `w[wi]`, using the
workspace `ws`, assembling each system matrix from the
[`HBLinearizedSystem`](@ref) `lsys` with [`assemblesystemmatrix!`](@ref)
and writing the results into the [`LinearizedArrays`](@ref) `arrays`
through per frequency views. An empty output array means that output was
not requested; small working matrices stand in for outputs which are
computed but not stored (`S` when only the quantum efficiency needs it).

`sensitivity` is a named tuple with the fixed operating point `stamps`,
the operating point `dAop` stamps of the forward contraction order, and
the [`ReverseSensitivity`](@ref) of the reverse order, or `nothing`.
`noiseplan`, `channeltemperatures` and `channelsigns` describe the noise
channels of the scattering blocks, the temperature of every channel and
the sign of each in the commutation relations (see
[`noisechannelsigns`](@ref)). `presolved`,
`presolvedadjoint` and `presolvednoise` are callbacks which replace the
assemble, factorize and solve of a frequency, the transposed solve, and
the noise scattering calculation with solutions computed elsewhere, which
is how the device sweep hands back its batches (see
[`devicesolutions`](@ref)).

Different frequency ranges may be computed in parallel: `lsys`,
`sensitivity` and `arrays` (through disjoint views) are shared, and each
task has its own `ws`.
"""
function hblinsolve_inner!(ws::LinearizedWorkspace, arrays::LinearizedArrays,
    sensitivity, lsys, bnm,
    portindices, noiseportimpedanceindices,
    portimpedances, noiseportimpedances, nodeindices,
    componenttypes, w, wpumpmodes, Nmodes, Nnodes, wi, factorization;
    noiseplan = nothing, channeltemperatures = nothing,
    channelsigns = nothing, presolved = nothing, presolvedadjoint = nothing,
    presolvednoise = nothing, refine::Bool = true)

    # everything downstream of the solve is the same whether the solution
    # came from the callbacks or from the factorization here
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
    scatteringwork = ws.scatteringwork
    occupation = ws.occupation
    rowsum = ws.rowsum
    rowcomp = ws.rowcomp

    # whether the transposed (adjoint) system must be solved: always for
    # the sensitivities, and otherwise when an output which reads the
    # adjoint solution (the noise scattering parameters, the quantum
    # efficiency, the commutation relations, the adjoint node outputs) was
    # requested
    needsadjoint = needsadjointsolve(arrays, noiseportimpedanceindices,
        noiseplan)

    # this worker's scattering block sensitivity state, bound before the
    # loop, where the signal frequency is `wsi`
    blocksens = ws.blocksens

    # The source current of each port mode in its own drive column, the only
    # source current its waves are credited with: a port which shares a node
    # with the driven one carries none.
    portsources = portsourcecurrents(bnm, portindices, nodeindices, Nmodes)
    portdrives = Diagonal(portsources)
    wantsS = !isempty(arrays.S) || !isempty(arrays.QE) ||
        !isempty(arrays.QEideal) || !isempty(arrays.CM) ||
        !isempty(arrays.Ssensitivity)
    # the scattering parameters and the noise scattering parameters both
    # divide by the incident waves of the port drives
    wantswaves = wantsS || !isempty(arrays.Snoise) || !isempty(arrays.Cnoise)

    for i in wi

        Sview = isempty(arrays.S) ? Sworking : view(arrays.S, :, :, i)
        Snoiseview = isempty(arrays.Snoise) ? Snoiseworking : view(arrays.Snoise, :, :, i)

        # the signal plus pump mode frequencies
        wsi = w[i]
        wmodes .= wsi .+ wpumpmodes

        # assemble the system matrix at this frequency,
        # Asparsecopy = AoLjnm + invLnm + im*Gnm*w - Cnm*w^2 with the per
        # column mode frequency, conjugating the negative frequency mode
        # entries and resolving the frequency dependent ones
        if isnothing(presolved)
            assemblesystemmatrix!(Asparsecopy, lsys, wmodes;
                scatteringwork = scatteringwork)

            # factorize the matrix, reusing the symbolic analysis; the
            # block factorization takes the mode count as its block size
            if factorization isa BlockFactorization
                tryfactorize!(cache, factorization, Asparsecopy;
                    blocksize = Nmodes, refine = refine ? 6 : 0)
            else
                tryfactorize!(cache, factorization, Asparsecopy)
            end

            trysolve!(phin, cache.factorization, bnm)
        else
            presolved(i, phin)
        end

        isempty(arrays.voltage) ||
            fluxtovoltage!(view(arrays.voltage, :, :, i), phin, wmodes, Nmodes)

        # copy the node fluxes for output; the auxiliary variables are
        # internal
        if !isempty(arrays.nodeflux)
            copy!(view(arrays.nodeflux,:,:,i), view(phin, 1:size(arrays.nodeflux,1), :))
        end

        if wantswaves
            calcinputwaves!(inputwave, portsources, portindices,
                portimpedances, componenttypes, wmodes)
        end

        # the scattering parameters
        if wantsS
            calcoutputwaves!(outputwave, phin, portdrives, portindices,
                portimpedances, nodeindices, componenttypes, wmodes)
            calcscatteringmatrix!(Sview, inputwave, outputwave)

            # the scalars which convert the adjoint contraction into the
            # scattering parameter derivatives
            if !isempty(arrays.Ssensitivity)
                calcsensitivityscaling!(sensitivitygamma, sensitivitybeta,
                    inputwave, portsources, portindices, portimpedances,
                    componenttypes, wmodes, Nmodes)
            end
        end

        if needsadjoint

            # keep the forward solution, which the adjoint solve overwrites
            if !isempty(arrays.Ssensitivity)
                copy!(phinforward, phin)
            end

            # Solve the transposed system, reusing the factorization. By the
            # adjoint identity the response at an output port to a source
            # anywhere in the circuit is that source contracted against the
            # transposed solution driven at the port, which is what the
            # noise and quantum efficiency calculations need. Without
            # scattering blocks this is also the solution with the conjugate
            # pump modulation matrix, related by the diagonal similarity of
            # `assemblesystemmatrix!`; the scattering rows break that
            # similarity, so with blocks it is the transposed system, whose
            # auxiliary port current rows the block noise channels read.
            if isnothing(presolvedadjoint)
                trysolvetranspose!(phin, cache.factorization, bnm)
            else
                presolvedadjoint(i, phin)
            end

            # copy the adjoint node fluxes for output
            if !isempty(arrays.nodefluxadjoint)
                copy!(view(arrays.nodefluxadjoint,:,:,i), view(phin, 1:size(arrays.nodefluxadjoint,1), :))
            end

            isempty(arrays.voltageadjoint) ||
                fluxtovoltage!(view(arrays.voltageadjoint, :, :, i), phin,
                    wmodes, Nmodes)

            # the noise scattering parameters. `presolvednoise` computes
            # them where the adjoint solution was computed and returns only
            # what the quantum efficiency and the commutation relations read,
            # two numbers per port mode rather than a row per noise channel;
            # the host forms the same reduction from its matrix below
            noise = ws.noise
            wantscnoise = !isempty(arrays.Cnoise)
            wantsreduction = !isempty(arrays.QE) || !isempty(arrays.CM)
            if !isempty(arrays.Snoise) || wantscnoise ||
                    !isempty(arrays.QE) || !isempty(arrays.CM)
                if isnothing(presolvednoise)
                    # the noise ports carry no source
                    calcoutputwaves!(noiseoutputwave, phin, nothing,
                        noiseportimpedanceindices, noiseportimpedances,
                        nodeindices, componenttypes, wmodes)
                    # the channels of the dissipative scattering blocks
                    # follow the lumped ones
                    if !isnothing(noiseplan)
                        scatteringnoisewaves!(noiseoutputwave, noiseplan,
                            lsys.scattering, phin, wmodes,
                            Nnoiseports*Nmodes, scatteringnoisework)
                    end
                    calcscatteringmatrix!(Snoiseview, inputwave, noiseoutputwave)
                    adjointnoisesigns!(Snoiseview, wmodes, Nmodes)
                else
                    noise = presolvednoise(i, inputwave, Snoiseview,
                        wantscnoise ? view(arrays.Cnoise,:,:,i) : nothing)
                end
            end

            # the scattering parameter sensitivities; `phin` now holds the
            # transposed solution the contraction needs
            if !isempty(arrays.Ssensitivity)
                # a scattering block stamp depends on the frequency through S
                # itself, so each worker rebuilds its own copy here
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
                    Nmodes)
                if !isnothing(sensitivity.reverse)
                    calcSsensitivityreverse!(view(arrays.Ssensitivity,:,:,:,i),
                        sensitivity.reverse, lsys, phinforward, phin,
                        sensitivitygamma, sensitivitybeta,
                        sensitivitycache, sensitivityrevbufs)
                end
            end

            # the occupation of each channel, which is where the temperature
            # enters; the noise scattering parameters carry none
            if wantsreduction || wantscnoise
                noiseoccupation!(occupation, channeltemperatures, wmodes,
                    Nmodes)
            end
            if wantscnoise && isnothing(presolvednoise)
                calcnoisecovariance!(view(arrays.Cnoise,:,:,i), Snoiseview,
                    occupation)
            end
            if wantsreduction && isnothing(presolvednoise)
                noisereduction!(noise, Snoiseview, wmodes, occupation,
                    channelsigns)
            end
            if !isempty(arrays.QE)
                calcqe!(view(arrays.QE,:,:,i), Sview, noise;
                    denom = rowsum, comp = rowcomp)
            end

            # the commutation relations
            if !isempty(arrays.CM)
                calccm!(view(arrays.CM,:,i), Sview, wmodes, noise;
                    comp = rowcomp)
            end
        else
            # the quantum efficiency and the commutation relations without
            # noise channels
            if !isempty(arrays.QE)
                calcqe!(view(arrays.QE,:,:,i), Sview;
                    denom = rowsum, comp = rowcomp)
            end
            if !isempty(arrays.CM)
                calccm!(view(arrays.CM,:,i), Sview, wmodes; comp = rowcomp)
            end
        end

        # the quantum efficiency of an ideal amplifier with the same gain
        if !isempty(arrays.QEideal)
            calcqeideal!(view(arrays.QEideal,:,:,i), Sview)
        end
    end
    return nothing
end
