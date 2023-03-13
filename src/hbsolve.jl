
"""
    NonlinearHB(nodeflux, Rbnm, Ljb, Lb, Ljbm, Nmodes, Nbranches, S)

A simple structure to hold the nonlinear harmonic balance solutions.

# Fields
- `w`: a tuple containing the the angular frequency of the pump in radians/s.
- `nodeflux`: the node fluxes resulting from inputs at each frequency and port.
- `Rbnm`: incidence matrix to convert between the node and branch basis.
- `Ljb`: sparse vector of Josephson junction inductances.
- `Lb`: sparse vector of linear inductances.
- `Ljbm`: sparse vector of linear inductances with each element duplicated Nmodes times.
- `Nmodes`: the number of signal and idler frequencies.
- `Nnodes`: the number of nodes in the circuit (including the ground node).
- `Nbranches`: the number of branches in the circuit.
- `nodes`: the vector of unique node strings.
- `ports`: vector of port numbers.
- `modes`: tuple of the pump mode indices where (1,) is the pump in the single
    pump case.
- `S`: the scattering matrix relating inputs and outputs for each combination
    of port and frequency (not currently functional).
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
end

"""
    LinearHB(S, Snoise, QE, QEideal, CM, nodeflux, voltage, Nmodes, Nnodes,
        Nbranches, signalindex, w)

A simple structure to hold the linearized harmonic balance solutions.

# Fields
- `S`: the scattering matrix relating inputs and outputs for each combination
    of port and frequency.
- `Snoise`: the scattering matrix relating inputs at the noise ports 
    (lossy devices) and outputs at the physical ports for each combination of
    port and frequency.
- `QE`: the quantum efficiency for each combination of port and frequency.
- `QEideal`: the quantum efficiency for an ideal amplifier with the same level
    of gain, for each combination of port and frequency. 
- `CM`: the commutation relations (equal to ±1), for each combination of port
    and frequency. 
- `nodeflux`: the node fluxes resulting from inputs at each frequency and port.
- `nodefluxadjoint`: the node fluxes resulting from inputs at each frequency and port.
- `voltage`: the node voltages resulting from inputs at each frequency and port.
- `Nmodes`: the number of signal and idler frequencies.
- `Nnodes`: the number of nodes in the circuit (including the ground node).
- `Nbranches`: the number of branches in the circuit.
- `nodes`: the vector of unique node strings.
- `ports`: vector of port numbers.
- `signalindex`: the index of the signal mode.
- `w`: the signal frequencies.
- `modes`: tuple of the signal mode indices where (0,) is the signal
"""
struct LinearHB
    S
    Snoise
    QE
    QEideal
    CM
    nodeflux
    nodefluxadjoint
    voltage
    Nmodes
    Nnodes
    Nbranches
    nodes
    ports
    signalindex
    w
    modes
end

"""
    HB(pump, Am, signal)

A simple structure to hold the nonlinear and linearized harmonic balance solutions.

# Fields
- `pump`: nonlinear harmonic balance solution for pump and pump harmonics
- `Am`: matrix encoding pump modulation
- `signal`: linearized harmonic balance solution
"""
struct HB
    pump
    signal
end


"""
    hbsolve(ws, wp, Ip, Nsignalmodes, Npumpmodes, circuit, circuitdefs;
        pumpports = [1], iterations = 1000, ftol = 1e-8,
        switchofflinesearchtol = 1e-5, alphamin = 1e-4,
        symfreqvar = nothing, nbatches = Base.Threads.nthreads(), sorting = :number,
        returnS = true, returnSnoise = false, returnQE = true, returnCM = true,
        returnnodeflux = false, returnvoltage = false, returnnodefluxadjoint = false,
        )

Calls the new harmonic balance solvers (which work for an arbitrary number of
modes and ports) `hbnlsolve` and `hblinsolve` using an identical syntax to
`hbsolveold`, which only supports four wave mixing processes involving single
strong tone and an arbitrary number of tone in the linearized solver. This
function is primarily for testing the new solvers and will eventually be
deprecated. 

This function attempts to mimic `hbsolveold`, but has a few important differences:
1. `Nmodulationharmonics` is no longer the number of signal modes but the number
    of modulation harmonics meaning `0` would give the signal only, `1` would
    three wave mixing processes, which this function does not enable. `2` would
    create idlers on either side of the signal mode at ws+2*wp and ws-2*wp.
    For other sets of idlers or three wave mixing processes, please use the new
    interface to `hbsolve`.
2. The outputs of the linearized harmonic balance solver `hblinsolve2` may not
    have the same ordering of signal modes as in `hblinsolve`. In `hblinsolve2`
    the signal mode is always at index 1 and the location of the other modes
    can be found by inspecting the contents of `modes`.
"""
function hbsolve(ws, wp, Ip, Nsignalmodes::Int, Npumpmodes::Int, circuit,
    circuitdefs; pumpports = [1], iterations = 1000, ftol = 1e-8,
    switchofflinesearchtol = 1e-5, alphamin = 1e-4,
    symfreqvar = nothing, nbatches = Base.Threads.nthreads(), sorting = :number,
    returnS = true, returnSnoise = false, returnQE = true, returnCM = true,
    returnnodeflux = false, returnvoltage = false, returnnodefluxadjoint = false,
    keyedarrays::Val{K} = Val(false)) where K

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

    # solve the nonlinear problem
    pump = hbnlsolve(w, sources, freq, indices, psc, cg, nm;
        iterations = iterations, x0 = nothing, ftol = ftol,
        switchofflinesearchtol = switchofflinesearchtol, alphamin = alphamin,
        symfreqvar = symfreqvar, keyedarrays = keyedarrays)

    # generate the signal modes
    signalfreq =truncfreqs(
        calcfreqsdft((Nsignalmodes,)),
        dc=true,odd=false,even=true,maxintermodorder=Inf,
    )

    # remove one of the signal modes if Nsignalmodes is even for compatibility
    # with old harmonic balance solver
    if mod(Nsignalmodes,2) == 0 && Nsignalmodes > 0
        signalfreq = JosephsonCircuits.removefreqs(signalfreq, [(Nsignalmodes,)])
    end

    # solve the linearized problem
    # i should make this a tuple
    signal = hblinsolve(ws, psc, cg, circuitdefs, signalfreq; pump = pump,
        symfreqvar = symfreqvar, nbatches = nbatches,
        returnS = returnS, returnSnoise = returnSnoise, returnQE = returnQE,
        returnCM = returnCM, returnnodeflux = returnnodeflux,
        returnnodefluxadjoint = returnnodefluxadjoint,
        returnvoltage = returnvoltage, keyedarrays = keyedarrays)

    return HB(pump, signal)
end


function hbsolve(ws, wp::NTuple{N,Any}, sources::Vector,
    Nmodulationharmonics::NTuple{M,Any}, Npumpharmonics::NTuple{N,Any},
    circuit, circuitdefs;dc = false, threewavemixing = false,
    fourwavemixing = true, maxintermodorder=Inf, iterations = 1000, ftol = 1e-8,
    switchofflinesearchtol = 1e-5, alphamin = 1e-4,
    symfreqvar = nothing, nbatches = Base.Threads.nthreads(), sorting = :number,
    returnS = true, returnSnoise = false, returnQE = true, returnCM = true,
    returnnodeflux = false, returnvoltage = false, returnnodefluxadjoint = false,
    keyedarrays::Val{K} = Val(true)) where {N,M,K}

    # calculate the frequency struct
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

    # solve the nonlinear problem
    pump = hbnlsolve(wp, sources, freq, indices, psc, cg, nm;
        iterations = iterations, x0 = nothing, ftol = ftol,
        switchofflinesearchtol = switchofflinesearchtol, alphamin = alphamin,
        symfreqvar = symfreqvar, keyedarrays = keyedarrays)

    # generate the signal modes
    signalfreq =truncfreqs(
        calcfreqsdft(Nmodulationharmonics),
        dc=true, odd=threewavemixing, even=fourwavemixing,
        maxintermodorder=maxintermodorder,
    )

    # solve the linearized problem
    # i should make this a tuple
    signal = hblinsolve(ws, psc, cg, circuitdefs, signalfreq; pump = pump,
        symfreqvar = symfreqvar, nbatches = nbatches,
        returnS = returnS, returnSnoise = returnSnoise, returnQE = returnQE,
        returnCM = returnCM, returnnodeflux = returnnodeflux,
        returnnodefluxadjoint = returnnodefluxadjoint,
        returnvoltage = returnvoltage, keyedarrays = keyedarrays)

    return HB(pump, signal)
end


"""
    hblinsolve(w, circuit,circuitdefs; Nmodulationharmonics = (0,),
        pump=nothing, symfreqvar=nothing, threewavemixing=false,
        fourwavemixing=true, maxintermodorder=Inf,
        nbatches::Integer = Base.Threads.nthreads(), returnS = true,
        returnSnoise = false, returnQE = true, returnCM = true,
        returnnodeflux = false, returnnodefluxadjoint = false, returnvoltage = false,
        )

Harmonic balance solver supporting an arbitrary number of small signals (weak
tones) linearized around `pump`, the solution of the nonlinear system consisting
of an arbitrary number of large signals (strong tones).

# Examples
```jldoctest
circuit = Array{Tuple{String,String,String,Union{Complex{Float64},Symbol,Int64}},1}(undef,0)
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

pump=hbnlsolve(
    (wp,),
    Npumpharmonics,
    [
        (mode=(0,),port=1,current=Idc),
        (mode=(1,),port=1,current=Ip),
    ],
    circuit,circuitdefs;dc=true,odd=fourwavemixing,even=threewavemixing)

signal = JosephsonCircuits.hblinsolve(ws,
    circuit, circuitdefs; Nmodulationharmonics = Nmodulationharmonics, pump = pump, symfreqvar=nothing,
    threewavemixing=false, fourwavemixing=true, returnnodeflux=true)
isapprox(signal.nodeflux,
    ComplexF64[9.901008591291e-12 - 6.40587007644028e-14im 2.164688307719963e-14 - 2.90852607344097e-16im 6.671563044645655e-14 - 8.585524364135119e-16im; 2.1633104519765224e-14 - 8.251861334047893e-16im 1.0099063486905209e-11 - 1.948847859339803e-13im -8.532003011745068e-15 + 3.234788465760295e-16im; 6.671648606599472e-14 + 7.892709980649199e-16im -8.53757633177974e-15 - 9.748395563374129e-17im 9.856580758892428e-12 + 5.859984004390703e-14im; 1.5888896262186103e-11 - 1.0303480614499543e-13im -2.557126237504446e-12 + 1.759201163407723e-14im -8.475819811683215e-12 + 5.3531443609574795e-14im; -2.5781681021577177e-13 + 4.757590640631487e-15im 2.36818731889176e-12 - 4.569646499606389e-14im 1.116372367616482e-13 - 2.039935997276492e-15im; -1.0210743447568219e-11 - 5.905490368441375e-14im 1.3377918536056493e-12 + 7.190105205618706e-15im 2.5392856657302323e-11 + 1.5143842454586225e-13im; 2.4781693042536835e-11 - 1.6057018472176702e-13im -2.5342360504077476e-12 + 1.7306764301173096e-14im -8.40554044664581e-12 + 5.269404591748149e-14im; -2.348528974341763e-13 + 3.949450668269274e-15im 1.1449271118157543e-11 - 2.2093702114766968e-13im 1.0261871618968225e-13 - 1.7240213938923877e-15im; -1.0140560031409567e-11 - 5.828587508192886e-14im 1.3288225860409326e-12 + 7.0954601524623594e-15im 3.423954321087654e-11 + 2.0403371894291513e-13im],
    atol = 1e-6)

# output
true
```
"""
function hblinsolve(w, circuit,circuitdefs; Nmodulationharmonics = (0,),
    pump=nothing, symfreqvar=nothing, threewavemixing=false,
    fourwavemixing=true, maxintermodorder=Inf,
    nbatches::Integer = Base.Threads.nthreads(), returnS = true,
    returnSnoise = false, returnQE = true, returnCM = true,
    returnnodeflux = false, returnnodefluxadjoint = false, returnvoltage = false,
    keyedarrays::Val{K} = Val(true)) where K

    # parse and sort the circuit
    psc = parsesortcircuit(circuit)

    # calculate the circuit graph
    cg = calccircuitgraph(psc)

    # generate the signal modes
    signalfreq =truncfreqs(
        calcfreqsdft(Nmodulationharmonics),
        dc=true,odd=threewavemixing,even=fourwavemixing,maxintermodorder=maxintermodorder,
    )

return hblinsolve(w, psc, cg, circuitdefs, signalfreq; pump = pump,
        symfreqvar = symfreqvar, nbatches = nbatches,
        returnS = returnS, returnSnoise = returnSnoise, returnQE = returnQE,
        returnCM = returnCM, returnnodeflux = returnnodeflux,
        returnnodefluxadjoint = returnnodefluxadjoint,
        returnvoltage = returnvoltage, keyedarrays = keyedarrays)
end

"""
    hblinsolve(w, psc::ParsedSortedCircuit,
        cg::CircuitGraph, circuitdefs, signalfreq::Frequencies{N}; pump=nothing,
        symfreqvar=nothing, nbatches::Integer = Base.Threads.nthreads(),
        returnS = true, returnSnoise = false, returnQE = true, returnCM = true,
        returnnodeflux = false, returnnodefluxadjoint = false, returnvoltage = false,
        )

Harmonic balance solver supporting an arbitrary number of small signals (weak
tones) linearized around `pump`, the solution of the nonlinear system consisting
of an arbitrary number of large signals (strong tones).

# Examples
```jldoctest
circuit = Array{Tuple{String,String,String,Union{Complex{Float64},Symbol,Int64}},1}(undef,0)
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
Npumpharmonics = (2,)
Nmodulationharmonics = (2,)
threewavemixing=false
fourwavemixing=true

frequencies = JosephsonCircuits.removeconjfreqs(
    JosephsonCircuits.truncfreqs(
        JosephsonCircuits.calcfreqsrdft(Npumpharmonics),
        dc=true, odd=true, even=false, maxintermodorder=Inf,
    )
)
fi = JosephsonCircuits.fourierindices(frequencies)
Nmodes = length(frequencies.modes)
psc = JosephsonCircuits.parsesortcircuit(circuit)
cg = JosephsonCircuits.calccircuitgraph(psc)
nm = JosephsonCircuits.numericmatrices(psc, cg, circuitdefs, Nmodes = Nmodes)
pump=hbnlsolve(
    (wp,),
    [
        (mode=(0,),port=1,current=Idc),
        (mode=(1,),port=1,current=Ip),
    ],
    frequencies, fi, psc, cg, nm)
signalfreq =JosephsonCircuits.truncfreqs(
    JosephsonCircuits.calcfreqsdft(Nmodulationharmonics),
    dc=true,odd=threewavemixing,even=fourwavemixing,maxintermodorder=Inf,
)
signal = JosephsonCircuits.hblinsolve(ws, psc, cg, circuitdefs, signalfreq;pump=pump,returnnodeflux=true)
isapprox(signal.nodeflux,
    ComplexF64[9.901008591291e-12 - 6.40587007644028e-14im 2.164688307719963e-14 - 2.90852607344097e-16im 6.671563044645655e-14 - 8.585524364135119e-16im; 2.1633104519765224e-14 - 8.251861334047893e-16im 1.0099063486905209e-11 - 1.948847859339803e-13im -8.532003011745068e-15 + 3.234788465760295e-16im; 6.671648606599472e-14 + 7.892709980649199e-16im -8.53757633177974e-15 - 9.748395563374129e-17im 9.856580758892428e-12 + 5.859984004390703e-14im; 1.5888896262186103e-11 - 1.0303480614499543e-13im -2.557126237504446e-12 + 1.759201163407723e-14im -8.475819811683215e-12 + 5.3531443609574795e-14im; -2.5781681021577177e-13 + 4.757590640631487e-15im 2.36818731889176e-12 - 4.569646499606389e-14im 1.116372367616482e-13 - 2.039935997276492e-15im; -1.0210743447568219e-11 - 5.905490368441375e-14im 1.3377918536056493e-12 + 7.190105205618706e-15im 2.5392856657302323e-11 + 1.5143842454586225e-13im; 2.4781693042536835e-11 - 1.6057018472176702e-13im -2.5342360504077476e-12 + 1.7306764301173096e-14im -8.40554044664581e-12 + 5.269404591748149e-14im; -2.348528974341763e-13 + 3.949450668269274e-15im 1.1449271118157543e-11 - 2.2093702114766968e-13im 1.0261871618968225e-13 - 1.7240213938923877e-15im; -1.0140560031409567e-11 - 5.828587508192886e-14im 1.3288225860409326e-12 + 7.0954601524623594e-15im 3.423954321087654e-11 + 2.0403371894291513e-13im],
    atol = 1e-6)

# output
true
```
"""
function hblinsolve(w, psc::ParsedSortedCircuit,
    cg::CircuitGraph, circuitdefs, signalfreq::Frequencies{N}; pump=nothing,
    symfreqvar=nothing, nbatches::Integer = Base.Threads.nthreads(),
    returnS = true, returnSnoise = false, returnQE = true, returnCM = true,
    returnnodeflux = false, returnnodefluxadjoint = false, returnvoltage = false,
    keyedarrays::Val{K} = Val(false)) where {N,K}

    Nsignalmodes = length(signalfreq.modes)
    # calculate the numeric matrices
    signalnm = numericmatrices(psc, cg, circuitdefs, Nmodes = Nsignalmodes)

    if pump == nothing

        allpumpfreq = calcfreqsrdft((0,))
        Amatrixmodes, Amatrixindices = hbmatind(allpumpfreq, signalfreq)
        Nwtuple = NTuple{length(allpumpfreq.Nw)+1,Int}((allpumpfreq.Nw...,length(signalnm.Ljb.nzval)))
        phimatrix = ones(Complex{Float64}, Nwtuple)
        wpumpmodes = calcmodefreqs((0.0,),signalfreq.modes)

    else

        pumpfreq = pump.frequencies

        # pumpfreq = JosephsonCircuits.calcfreqsrdft((2*Npumpmodes,))
        allpumpfreq = calcfreqsrdft(pumpfreq.Nharmonics)
        pumpindices = fourierindices(pumpfreq)
        Npumpmodes = length(pumpfreq.modes)

        Amatrixmodes, Amatrixindices = hbmatind(allpumpfreq, signalfreq)

        # calculate the dimensions of the array which holds the frequency
        # domain information for the fourier transform
        Nwtuple = NTuple{length(pumpfreq.Nw)+1,Int}((pumpfreq.Nw...,length(pump.Ljb.nzval)))

        # create an array to hold the frequency domain data for the
        # fourier transform
        phimatrix = zeros(Complex{Float64}, Nwtuple)

        # create an array to hold the time domain data for the RFFT. also generate
        # the plans.
        phimatrixtd, irfftplan, rfftplan = plan_applynl(phimatrix)

        # convert the branch flux vector to a matrix with the terms arranged
        # in the correct way for the inverse rfft including the appropriate
        # complex conjugates.
        branchflux = pump.Rbnm*pump.nodeflux
        phivectortomatrix!(
            branchflux[pump.Ljbm.nzind], phimatrix,
            pumpindices.vectomatmap,
            pumpindices.conjsourceindices,
            pumpindices.conjtargetindices,
            length(pump.Ljb.nzval)
        )

        # apply the sinusoidal nonlinearity when evaluaing the function
        applynl!(
            phimatrix,
            phimatrixtd,
            cos,
            irfftplan,
            rfftplan,
        )

        wpumpmodes = calcmodefreqs(pump.w,signalfreq.modes)
    end

    # this is the first signal frequency. we will use it for various setup tasks
    wmodes = w[1] .+ wpumpmodes
    wmodesm = Diagonal(repeat(wmodes,outer=psc.Nnodes-1));
    wmodes2m = Diagonal(repeat(wmodes.^2,outer=psc.Nnodes-1));

    Nfreq = prod(signalfreq.Nw)

    AoLjbmindices, conjindicessorted = calcAoLjbmindices(
        Amatrixindices, signalnm.Ljb, Nsignalmodes, cg.Nbranches, Nfreq
    )

    AoLjbm = calcAoLjbm2(phimatrix, 
        Amatrixindices, signalnm.Ljb, 1, Nsignalmodes, cg.Nbranches
    )

    AoLjnm = signalnm.Rbnm'*AoLjbm*signalnm.Rbnm

    # extract the elements we need
    Nnodes = psc.Nnodes
    typevector = psc.typevector
    nodeindexarraysorted = psc.nodeindexarraysorted
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

    # calculate the source currents
    Nports = length(portindices)

    # calculate the source terms in the branch basis
    bbm = zeros(Complex{Float64},Nbranches*Nsignalmodes,Nsignalmodes*Nports)

    # add a current source for each port and mode
    for (i,val) in enumerate(portindices)
        key = (nodeindexarraysorted[1,val],nodeindexarraysorted[2,val])
        for j = 1:Nsignalmodes
            bbm[(edge2indexdict[key]-1)*Nsignalmodes+j,(i-1)*Nsignalmodes+j] = 1
            # bbm2[(i-1)*Nmodes+j,(i-1)*Nmodes+j] = Lmean*1/phi0
        end
    end

    # calculate the source terms in the node basis
    bnm = transpose(Rbnm)*bbm

    # if there is a symbolic frequency variable, then we need to redo the noise
    # port calculation because calcnoiseportimpedanceindices() can't tell if a
    # symbolic expression is complex. 
    if symfreqvar == nothing
        noiseportimpedanceindices = signalnm.noiseportimpedanceindices
    else
        noiseportimpedanceindices = calcnoiseportimpedanceindices(
            psc.typevector, psc.nodeindexarraysorted,
            psc.mutualinductorvector,
            Symbolics.substitute.(signalnm.vvn, symfreqvar => wmodes[1]))
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


    # make arrays for the voltages, node fluxes, scattering parameters,
    # quantum efficiency, and commutatio relations. if we aren't returning a
    # matrix, set it to be an empty array.
    if returnS
        S = zeros(Complex{Float64},Nports*Nsignalmodes,Nports*Nsignalmodes,length(w))
    else
        S = zeros(Complex{Float64},0,0,0)
    end

    if returnSnoise
        Snoise = zeros(Complex{Float64},length(noiseportimpedanceindices)*Nsignalmodes,Nports*Nsignalmodes,length(w))
    else
        Snoise = zeros(Complex{Float64},0,0,0)
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
        nodeflux = zeros(Complex{Float64},Nsignalmodes*(Nnodes-1),Nsignalmodes*Nports,length(w))
    else
        nodeflux = Vector{Complex{Float64}}(undef,0)
    end

    if returnnodefluxadjoint
        nodefluxadjoint = zeros(Complex{Float64},Nsignalmodes*(Nnodes-1),Nsignalmodes*Nports,length(w))
    else
        nodefluxadjoint = Vector{Complex{Float64}}(undef,0)
    end

    if returnvoltage
        voltage = zeros(Complex{Float64},Nsignalmodes*(Nnodes-1),Nsignalmodes*Nports,length(w))
    else
        voltage = Vector{Complex{Float64}}(undef,0)
    end

    # generate the mode indices and find the signal index
    signalindex = 1
    Nmodes = Nsignalmodes

    # solve the linear system for the specified frequencies. the response for
    # each frequency is independent so it can be done in parallel; however
    # we want to reuse the factorization object and other input arrays. 
    # perform array allocations and factorization "nbatches" times.
    # parallelize using native threading
    batches = collect(Base.Iterators.partition(1:length(w),1+(length(w)-1)÷nbatches))
    Base.Threads.@threads for i in 1:length(batches)
        hblinsolve_inner!(S, Snoise, QE, CM, nodeflux, nodefluxadjoint, voltage, Asparse,
            AoLjnm, invLnm, Cnm, Gnm, bnm,
            AoLjnmindexmap, invLnmindexmap, Cnmindexmap, Gnmindexmap,
            Cnmfreqsubstindices, Gnmfreqsubstindices, invLnmfreqsubstindices,
            portindices, portimpedanceindices, noiseportimpedanceindices,
            portimpedances, noiseportimpedances, nodeindexarraysorted, typevector,
            w, wpumpmodes, Nmodes, Nnodes, symfreqvar, batches[i])
    end

    if returnQE
        # calculate the ideal quantum efficiency.
        QEideal = zeros(Float64,size(S))
        calcqeideal!(QEideal,S)
    else
        QEideal = zeros(Float64,0,0,0)
    end

    # if keyword argument keyedarrays = Val(true) then generate keyed arrays
    if returnS && K
        Sout = Stokeyed(S, signalfreq.modes, portnumbers, signalfreq.modes,
            portnumbers, w)
    else
        Sout = S
    end

    if returnQE && K
        QEout = Stokeyed(QE, signalfreq.modes, portnumbers, signalfreq.modes,
            portnumbers, w)
        QEidealout = Stokeyed(QEideal, signalfreq.modes, portnumbers, signalfreq.modes,
            portnumbers, w)
    else
        QEout = QE
        QEidealout = QEideal
    end

    return LinearHB(Sout, Snoise, QEout, QEidealout, CM, nodeflux, nodefluxadjoint,
        voltage, Nmodes, Nnodes, Nbranches, psc.uniquenodevectorsorted[2:end],
        portnumbers, signalindex, w, signalfreq.modes)
end


"""
    hblinsolve_inner!(S, Snoise, QE, CM, nodeflux, voltage, Asparse,
        AoLjnm, invLnm, Cnm, Gnm, bnm,
        AoLjnmindexmap, invLnmindexmap, Cnmindexmap, Gnmindexmap,
        Cnmfreqsubstindices, Gnmfreqsubstindices, invLnmfreqsubstindices,
        portindices, portimpedanceindices, noiseportimpedanceindices,
        portimpedances, noiseportimpedances, nodeindexarraysorted, typevector,
        w, indices, wp, Nmodes, Nnodes, symfreqvar, wi)

Solve the linearized harmonic balance problem for a subset of the frequencies
given by `wi`. This function is thread safe in that different frequencies can
be computed in parallel on separate threads.
"""
function hblinsolve_inner!(S, Snoise, QE, CM, nodeflux, nodefluxadjoint, voltage, Asparse,
    AoLjnm, invLnm, Cnm, Gnm, bnm,
    AoLjnmindexmap, invLnmindexmap, Cnmindexmap, Gnmindexmap,
    Cnmfreqsubstindices, Gnmfreqsubstindices, invLnmfreqsubstindices,
    portindices, portimpedanceindices, noiseportimpedanceindices,
    portimpedances, noiseportimpedances, nodeindexarraysorted, typevector,
    w, wpumpmodes, Nmodes, Nnodes, symfreqvar, wi)

    Nports = length(portindices)
    phin = zeros(Complex{Float64},Nmodes*(Nnodes-1),Nmodes*Nports)
    inputwave = Diagonal(zeros(Complex{Float64},Nports*Nmodes))
    outputwave = zeros(Complex{Float64},Nports*Nmodes,Nports*Nmodes)

    Nnoiseports = length(noiseportimpedanceindices)
    noiseoutputwave = zeros(Complex{Float64},Nnoiseports*Nmodes,Nports*Nmodes)

    # operate on a copy of Asparse because it may be modified by multiple threads
    # at the same time.
    Asparsecopy = copy(Asparse)

    # calculate the conjugate of AoLjnm
    AoLjnmconj = copy(AoLjnm)
    conj!(AoLjnmconj.nzval)

    # if using the KLU factorization and sparse solver then make a 
    # factorization for the sparsity pattern.
    cache = FactorizationCache(KLU.klu(Asparsecopy))

    # if the scattering matrix is empty define a new working matrix
    if isempty(S)
        Sview = zeros(Complex{Float64},Nports*Nmodes,Nports*Nmodes)
    end

    # if the noise scattering matrix is empty define a new working matrix
    if isempty(Snoise)
        Snoiseview = zeros(Complex{Float64},Nnoiseports*Nmodes,Nports*Nmodes)
    end

    # loop over the frequencies
    for i in wi

        # if the scattering matrix is not empty define a view
        if !isempty(S)
            Sview = view(S,:,:,i)
        end

        # if the noise scattering matrix is not empty define a view
        if !isempty(Snoise)
            Snoiseview = view(Snoise,:,:,i)
        end

        # calculate the frequency matrices
        ws = w[i]
        # wmodes = calcw(ws,indices,wp);
        wmodes = ws .+ wpumpmodes
        wmodesm = Diagonal(repeat(wmodes,outer=Nnodes-1));
        wmodes2m = Diagonal(repeat(wmodes.^2,outer=Nnodes-1));

        # perform the operation below in a way that doesn't allocate significant
        # memory, plus take the conjugates mentioned below.
        # Asparsecopy = (AoLjnm + invLnm - im.*Gnm*wmodesm - Cnm*wmodes2m)

        fill!(Asparsecopy.nzval,zero(eltype(Asparsecopy.nzval)))
        sparseadd!(Asparsecopy,1,AoLjnm,AoLjnmindexmap)

        # take the complex conjugate of the negative frequency terms in
        # the capacitance and conductance matrices. substitute in the symbolic
        # frequency variable if present.
        sparseaddconjsubst!(Asparsecopy,-1,Cnm,wmodes2m,Cnmindexmap,
            wmodesm .< 0,wmodesm,Cnmfreqsubstindices,symfreqvar)
        sparseaddconjsubst!(Asparsecopy,im,Gnm,wmodesm,Gnmindexmap,
            wmodesm .< 0,wmodesm,Gnmfreqsubstindices,symfreqvar)
        sparseaddconjsubst!(Asparsecopy,1,invLnm,Diagonal(ones(size(invLnm,1))),invLnmindexmap,
            wmodesm .< 0,wmodesm,invLnmfreqsubstindices,symfreqvar)

        # factor the sparse matrix
        factorklu!(cache, Asparsecopy)

        # solve the linear system
        ldiv!(phin,cache.factorization,bnm)

        # convert to node voltages. node flux is defined as the time integral of
        # node voltage so node voltage is derivative of node flux which can be
        # accomplished in the frequency domain by multiplying by j*w.
        if !isempty(voltage)
            voltage[:,:,i] .= im*wmodesm*phin
        end

        # copy the nodeflux for output
        if !isempty(nodeflux)
            copy!(view(nodeflux,:,:,i),phin)
        end

        # calculate the scattering parameters
        if !isempty(S) || !isempty(QE) || !isempty(CM)
            calcS!(Sview,inputwave,outputwave,phin,bnm,portimpedanceindices,portimpedanceindices,
                portimpedances,portimpedances,nodeindexarraysorted,typevector,wmodes,symfreqvar)
        end

        if (Nnoiseports > 0 || !isempty(nodefluxadjoint)) && (!isempty(Snoise) || !isempty(QE) || !isempty(CM) || !isempty(nodefluxadjoint))

            # solve the nonlinear system with the complex conjugate of the pump
            # modulation matrix
            fill!(Asparsecopy.nzval,zero(eltype(Asparsecopy.nzval)))
            sparseadd!(Asparsecopy,1,AoLjnmconj,AoLjnmindexmap)

            # take the complex conjugate of the negative frequency terms in
            # the capacitance and conductance matrices. substitute in the symbolic
            # frequency variable if present. 
            sparseaddconjsubst!(Asparsecopy,-1,Cnm,wmodes2m,Cnmindexmap,
                wmodesm .< 0,wmodesm,Cnmfreqsubstindices,symfreqvar)
            sparseaddconjsubst!(Asparsecopy,im,Gnm,wmodesm,Gnmindexmap,
                wmodesm .< 0,wmodesm,Gnmfreqsubstindices,symfreqvar)
            sparseaddconjsubst!(Asparsecopy,1,invLnm,Diagonal(ones(size(invLnm,1))),invLnmindexmap,
                wmodesm .< 0,wmodesm,invLnmfreqsubstindices,symfreqvar)

            # factor the sparse matrix
            factorklu!(cache, Asparsecopy)

            # solve the linear system
            ldiv!(phin,cache.factorization,bnm)

            # copy the nodeflux adjoint for output
            if !isempty(nodefluxadjoint)
                copy!(view(nodefluxadjoint,:,:,i),phin)
            end

            # calculate the noise scattering parameters
            if !isempty(Snoise)  || !isempty(QE) || !isempty(CM)
                calcSnoise!(Snoiseview,inputwave,noiseoutputwave,phin,bnm,portimpedanceindices,noiseportimpedanceindices,
                    portimpedances,noiseportimpedances,nodeindexarraysorted,typevector,wmodes,symfreqvar)
            end

            # calculate the quantum efficiency
            if !isempty(QE)
                calcqe!(view(QE,:,:,i),Sview,transpose(Snoiseview))
            end

            # calculate the commutation relations (Manley-Rowe relations)
            if !isempty(CM)
                calccm!(view(CM,:,i),Sview,transpose(Snoiseview), wmodes)
            end
        else
            # calculate the quantum efficiency
            if !isempty(QE)
                calcqe!(view(QE,:,:,i),Sview)
            end

            # calculate the commutation relations (Manley-Rowe relations)
            if !isempty(CM)
                calccm!(view(CM,:,i),Sview, wmodes)
            end
        end
    end
    return nothing
end


"""
    hbnlsolve(w::NTuple{N,Any}, Nharmonics::NTuple{N,Int}, sources,
        circuit, circuitdefs; iterations = 1000,
        maxintermodorder = Inf, dc = false, odd = true, even = false, x0 = nothing,
        ftol = 1e-8, switchofflinesearchtol = 1e-5, alphamin = 1e-4,
        symfreqvar = nothing, sorting= :number)

New version of the nonlinear harmonic balance solver suitable for arbitrary
numbers of ports, sources, and drives including direct current (zero frequency)
or flux pumping using a current source and a mutual inductor.

# Examples
```jldoctest
circuit = Array{Tuple{String,String,String,Union{Complex{Float64},Symbol,Int64}},1}(undef,0)
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
isapprox(out.nodeflux,
    ComplexF64[15.190314040027522 - 8.56492651167657e-24im, 2.991103820177504e-6 - 1.8501001011477133e-8im, -6.835392148510984 - 1.0356102442254259e-14im, 7.396422335315908e-6 - 4.5749403967992827e-8im, 6.835392148539885 - 1.0356102451770844e-14im, 1.008026285172782e-5 - 6.23498762664213e-8im],
    atol = 1e-6)

# output
true
```
"""
function hbnlsolve(w::NTuple{N,Any}, Nharmonics::NTuple{N,Int}, sources,
    circuit, circuitdefs; iterations = 1000,
    maxintermodorder = Inf, dc = false, odd = true, even = false, x0 = nothing,
    ftol = 1e-8, switchofflinesearchtol = 1e-5, alphamin = 1e-4,
    symfreqvar = nothing, sorting= :number, keyedarrays::Val{K} = Val(false)) where {N,K}

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

    return hbnlsolve(w, sources, freq, indices, psc, cg, nm;
        iterations = iterations, x0 = x0, ftol = ftol,
        switchofflinesearchtol = switchofflinesearchtol, alphamin = alphamin,
        symfreqvar = symfreqvar, keyedarrays = keyedarrays)
end

"""
    hbnlsolve(w::NTuple{N,Any}, sources, frequencies::Frequencies{N},
        indices::FourierIndices{N}, psc::ParsedSortedCircuit, cg::CircuitGraph,
        nm::CircuitMatrices; iterations = 1000, x0 = nothing,
        ftol = 1e-8, switchofflinesearchtol = 1e-5, alphamin = 1e-4,
        symfreqvar = nothing)

New version of the nonlinear harmonic balance solver suitable for arbitrary
numbers of ports, sources, and drives including direct current (zero frequency)
or flux pumping using a current source and a mutual inductor.

# Examples
```jldoctest
circuit = Array{Tuple{String,String,String,Union{Complex{Float64},Symbol,Int64}},1}(undef,0)
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
isapprox(out.nodeflux,
    ComplexF64[15.190314040027522 - 8.56492651167657e-24im, 2.991103820177504e-6 - 1.8501001011477133e-8im, -6.835392148510984 - 1.0356102442254259e-14im, 7.396422335315908e-6 - 4.5749403967992827e-8im, 6.835392148539885 - 1.0356102451770844e-14im, 1.008026285172782e-5 - 6.23498762664213e-8im],
    atol = 1e-6)

# output
true
```
"""
function hbnlsolve(w::NTuple{N,Any}, sources, frequencies::Frequencies{N},
    indices::FourierIndices{N}, psc::ParsedSortedCircuit, cg::CircuitGraph,
    nm::CircuitMatrices; iterations = 1000, x0 = nothing,
    ftol = 1e-8, switchofflinesearchtol = 1e-5, alphamin = 1e-4,
    symfreqvar = nothing, keyedarrays::Val{K} = Val(false)) where {N,K}

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
    typevector = psc.typevector
    nodeindexarraysorted = psc.nodeindexarraysorted
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

    # calculate the diagonal frequency matrices
    wmodesm = Diagonal(repeat(wmodes, outer = Nnodes-1))
    wmodes2m = Diagonal(repeat(wmodes.^2, outer = Nnodes-1))

    # calculate the source terms in the branch basis
    bbm = calcsources(modes, sources, portindices, portnumbers,
        nodeindexarraysorted, edge2indexdict, Lmean, Nnodes, Nbranches, Nmodes)

    # convert from the node basis to the branch basis
    bnm = transpose(Rbnm)*bbm

    # calculate the dimensions of the array which holds the frequency
    # domain information for the fourier transform
    Nwtuple = NTuple{length(Nw)+1,Int}((Nw...,length(Ljb.nzval)))

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
    AoLjbmindices, conjindicessorted = calcAoLjbmindices(Amatrixindices,
        Ljb, Nmodes, Nbranches, Nfreq)

    # right now i redo the calculation of AoLjbmindices, conjindicessorted in calcAoLjbm2
    AoLjbm = calcAoLjbm2(Amatrix, Amatrixindices, Ljb, Lmean, Nmodes, Nbranches)
    AoLjbmcopy = calcAoLjbm2(Amatrix, Amatrixindices, Ljb, Lmean, Nmodes, Nbranches)

    # convert to a sparse node matrix. Note: I was having problems with type 
    # instability when i used AoLjbm here instead of AoLjbmcopy. 
    AoLjnm = transpose(Rbnm)*AoLjbmcopy*Rbnm;

    if x0 == nothing
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
    # J .= AoLjnm + invLnm + im*Gnm*wmodesm - Cnm*wmodes2m
    # J = spaddkeepzeros(spaddkeepzeros(spaddkeepzeros(AoLjnm, invLnm), im*Gnm*wmodesm), - Cnm*wmodes2m)
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
    function fj!(F, J, x)
        calcfj2!(F, J, x, wmodesm, wmodes2m, Rbnm, Rbnmt, invLnm,
            Cnm, Gnm, bnm, Ljb, Ljbm, Nmodes,
            Nbranches, Lmean, AoLjbmvector, AoLjbm,
            AoLjnmindexmap, invLnmindexmap, Gnmindexmap, Cnmindexmap,
            AoLjbmindices, conjindicessorted,
            freqindexmap, conjsourceindices, conjtargetindices, phimatrix,
            AoLjnm, xbAoLjnm, AoLjbmRbnm, xbAoLjbmRbnm,
            phimatrixtd, irfftplan, rfftplan,
        )
        return nothing
    end

    # # use this for debugging purposes to return the function and the
    # # Jacobian
    # fj!(F,J,x)
    # return (F,J,x)

    # solve the nonlinear system
    nlsolve!(fj!, F, J, x; iterations = iterations, ftol = ftol,
        switchofflinesearchtol = switchofflinesearchtol, alphamin = alphamin)

    nodeflux = x

    # calculate the scattering parameters for the pump
    Nports = length(portindices)
    S = zeros(Complex{Float64}, Nports*Nmodes, Nports*Nmodes)
    inputwave = zeros(Complex{Float64}, Nports*Nmodes)
    outputwave = zeros(Complex{Float64},Nports*Nmodes)
    portimpedances = [vvn[i] for i in portimpedanceindices]
    if !isempty(S)
        calcS!(S,inputwave,outputwave,nodeflux,bnm/Lmean,portimpedanceindices,portimpedanceindices,
            portimpedances,portimpedances,nodeindexarraysorted,typevector,wmodes,symfreqvar)
    end

    #
    if !isempty(S) && K
        Sout = Stokeyed(S, modes, portnumbers, modes, portnumbers)
    else
        Sout = S
    end

    return NonlinearHB(w, frequencies, nodeflux, Rbnm, Ljb, Lb, Ljbm, Nmodes,
        Nbranches, psc.uniquenodevectorsorted[2:end], portnumbers, modes, Sout)

end


"""
    calcfj2(F,J,phin,wmodesm,wmodes2m,Rbnm,invLnm,Cnm,Gnm,bm,Ljb,Ljbindices,
        Ljbindicesm,Nmodes,Lmean,AoLjbm)
        
Calculate the residual and the Jacobian. These are calculated with one function
in order to reuse as much as possible.

Leave off the type signatures on F and J because the solver will pass a type of
Nothing if it only wants to calculate F or J. 

"""
function calcfj2!(F,
        J,
        nodeflux::AbstractVector,
        wmodesm::AbstractMatrix,
        wmodes2m::AbstractMatrix,
        Rbnm::AbstractArray{Int, 2},
        Rbnmt::AbstractArray{Int, 2},
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
        )

    # convert from a node flux to a branch flux
    phib = Rbnm*nodeflux

    if !(F == nothing)

        # convert the branch flux vector to a matrix with the terms arranged
        # in the correct way for the inverse rfft including the appropriate
        # complex conjugates.
        phivectortomatrix!(phib[Ljbm.nzind], phimatrix, freqindexmap,
            conjsourceindices, conjtargetindices, length(Ljb.nzval))

        # apply the sinusoidal nonlinearity when evaluaing the function
        applynl!(phimatrix, phimatrixtd, sin, irfftplan, rfftplan)

        # convert the sinphimatrix to a vector
        fill!(AoLjbmvector, 0)
        AoLjbmvectorview = view(AoLjbmvector, Ljbm.nzind)
        phimatrixtovector!(AoLjbmvectorview, phimatrix, freqindexmap,
            conjsourceindices, conjtargetindices, length(Ljb.nzval))

        for i in eachindex(AoLjbmvectorview)
            AoLjbmvectorview[i] = AoLjbmvectorview[i] * (Lmean/Ljbm.nzval[i])
        end
        # TODO: for multi-tone harmonic balance this needs conjugates (when any
        # of the frequencies are negative)
        F .= Rbnmt*AoLjbmvector .+ invLnm*nodeflux .+ im*Gnm*wmodesm*nodeflux .- Cnm*wmodes2m*nodeflux .- bnm

    end

    #calculate the Jacobian
    if !(J == nothing)

        # turn the phivector into a matrix again because applynl! overwrites
        # the frequency domain data
        phivectortomatrix!(phib[Ljbm.nzind], phimatrix, freqindexmap,
            conjsourceindices, conjtargetindices, length(Ljb.nzval))

        # apply a cosinusoidal nonlinearity when evaluating the Jacobian
        applynl!(phimatrix, phimatrixtd, cos, irfftplan, rfftplan)

        # calculate  AoLjbm
        updateAoLjbm2!(AoLjbm, phimatrix, AoLjbmindices, conjindicessorted,
            Ljb, Lmean)

        # convert to a sparse node matrix
        # AoLjnm = Rbnmt*AoLjbm*Rbnm
        # non allocating sparse matrix multiplication
        spmatmul!(AoLjbmRbnm, AoLjbm, Rbnm, xbAoLjbmRbnm)
        spmatmul!(AoLjnm, Rbnmt, AoLjbmRbnm, xbAoLjnm)

        # calculate the Jacobian. If J is sparse, keep it sparse. 
        # J .= AoLjnm + invLnm + im*Gnm*wmodesm - Cnm*wmodes2m
        # the code below adds the sparse matrices together with minimal
        # memory allocations and without changing the sparsity structure.
        fill!(J, 0)
        # TODO: for multi-tone harmonic balance this needs conjugates (when 
        # any of the frequencies are negative)
        sparseadd!(J, AoLjnm, AoLjnmindexmap)
        sparseadd!(J, invLnm, invLnmindexmap)
        sparseadd!(J, im, Gnm, wmodesm, Gnmindexmap)
        sparseadd!(J, -1, Cnm, wmodes2m, Cnmindexmap)
    end
    return nothing
end

"""
    calcAoLjbmindices(Amatrixindices,Ljb::SparseVector,Nmodes,Nbranches,Nfreq)

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
```jldoctest
@syms A11 A12 A21 A22 A31 A32 Lj1 Lj2
Amatrix = [A11 A12;A21 A22;A31 A32]
Amatrixindices = [1 0 0; 0 1 0; 0 0 1]
Ljb = JosephsonCircuits.SparseArrays.sparsevec([1,2],[Lj1,Lj2])
Lmean = 1
Nmodes = 3
Nbranches = 2
JosephsonCircuits.calcAoLjbm2(Amatrix, Amatrixindices, Ljb, Lmean, Nmodes, Nbranches).nzval

# output
6-element Vector{Any}:
 A11 / Lj1
 A11 / Lj1
 A11 / Lj1
 A12 / Lj2
 A12 / Lj2
 A12 / Lj2
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

Update the values in the sparse AoLjbm matrix in place.

"""
function updateAoLjbm2!(AoLjbm::SparseMatrixCSC,Am::Array, AoLjbmindices,
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
    calcsources(modes, sources, portindices, portnumbers, nodeindexarraysorted,
        edge2indexdict, Lmean, Nnodes, Nbranches, Nmodes)

Calculate the source terms in the branch basis. See also [`addsources!`](@ref).

# Examples
```jldoctest
modes = [(0,), (1,)]
sources = [(mode = (0,), port = 1, current = 0.0005), (mode = (1,), port = 1, current = 1.0e-10)]
portindices = [1]
portnumbers = [1]
nodeindexarraysorted = [2 2 2 2 0 2 3 4 3 3; 1 1 1 1 0 3 4 1 1 1]
edge2indexdict = Dict((1, 2) => 1, (3, 1) => 2, (1, 3) => 2, (4, 1) => 3, (2, 1) => 1, (1, 4) => 3, (3, 4) => 4, (4, 3) => 4)
Lmean = 1.005e-9 + 0.0im
Nnodes = 4
Nbranches = 4
Nmodes = 2
JosephsonCircuits.calcsources(modes, sources, portindices, portnumbers,
    nodeindexarraysorted, edge2indexdict, Lmean, Nnodes, Nbranches, Nmodes)

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
function calcsources(modes, sources, portindices, portnumbers, nodeindexarraysorted,
    edge2indexdict, Lmean, Nnodes, Nbranches, Nmodes)

    bbm = zeros(Complex{Float64}, Nbranches*Nmodes)

    addsources!(bbm, modes, sources, portindices, portnumbers,
        nodeindexarraysorted, edge2indexdict, Lmean, Nnodes, Nbranches,
        Nmodes)

    return bbm
end

"""
    addsources!(bbm, modes, sources, portindices, portnumbers,
        nodeindexarraysorted, edge2indexdict, Lmean, Nnodes, Nbranches, Nmodes)

Calculate the source terms in the branch basis. Overwrite bbm with the output.
See also [`calcsources`](@ref).
"""
function addsources!(bbm, modes, sources, portindices, portnumbers,
    nodeindexarraysorted, edge2indexdict, Lmean, Nnodes, Nbranches, Nmodes)

    # fill the vector with zeros
    fill!(bbm,0)

    # make a dictionary of ports
    portdict = Dict{eltype(portnumbers),eltype(portindices)}()
    for i in eachindex(portindices)
        portdict[portnumbers[i]] = portindices[i]
    end

    # make a dictionary of modes
    modedict = Dict{eltype(modes),Int}()
    for i in eachindex(modes)
        modedict[modes[i]] = i
    end

    for source in sources
        # pull out the necessary values from the named tuple
        port = source[:port]
        mode = source[:mode]
        current = source[:current]

        # check if the port is in the dictionary of ports
        if haskey(portdict,port)
            portindex = portdict[port]
            # check if the mode is in the dictionary of modes
            if haskey(modedict,mode)
                # if we find the mode and the port, set that branch in bbm
                # equal to the current scaled by the mean inductance and the
                # flux quantum.
                modeindex = modedict[mode]
                key = (nodeindexarraysorted[1, portindex], nodeindexarraysorted[2, portindex])
                bbm[(edge2indexdict[key]-1)*Nmodes+modeindex] += Lmean*current/phi0
            else
                throw(ArgumentError("Source mode $(mode) not found."))
            end
        else
            throw(ArgumentError("Source port $(port) not found."))
        end
    end

    return nothing
end


"""
    Stokeyed(S, outputmodes, outputportnumbers, inputmodes, inputportnumbers, w)

Convert a scattering parameter array to a keyed array. 

"""
function Stokeyed(S, outputmodes, outputportnumbers, inputmodes,
    inputportnumbers, w)
    Nfrequencies = length(w)
    
    return AxisKeys.KeyedArray(
        reshape(S, length(outputmodes), length(outputportnumbers),
            length(inputmodes), length(inputportnumbers), Nfrequencies),
        outputmode = outputmodes,
        outputport = outputportnumbers,
        inputmode = inputmodes,
        inputport = inputportnumbers,
        freqindex=1:Nfrequencies,
    )
end

function Stokeyed(S, outputmodes, outputportnumbers, inputmodes,
    inputportnumbers)
    
    return AxisKeys.KeyedArray(
        reshape(S, length(outputmodes), length(outputportnumbers),
            length(inputmodes), length(inputportnumbers)),
        outputmode = outputmodes,
        outputport = outputportnumbers,
        inputmode = inputmodes,
        inputport = inputportnumbers,
    )
end


# function scatteringparamstoarray(S,modes,portnumbers,w)
#     Nmodes = length(modes)
#     Nports = length(portnumbers)
#     Nfrequencies = length(w)
    
#     return KeyedArray(
#         reshape(S, Nmodes, Nports, Nmodes, Nports, Nfrequencies),
#         outputmode = modes,
#         outputport = portnumbers,
#         inputmode = modes,
#         inputport = portnumbers,
#         freqindex=1:Nfrequencies,
#     )
    
# end


# function fieldstokeyed(nodeflux, inputmodes, outputmodes,
#     inputportnumbers, outputportnumbers)

# nodeflux = KeyedArray(
#     reshape(
#         result.signal.nodeflux,
#         result.signal.Nmodes,
#         result.signal.Nnodes-1,
#         result.signal.Nmodes,
#         2,
#         length(result.signal.w),
#     ),
#     outputmode = result.signal.modes,
#     # i need to replace this with the node strings
#     # i also need to return the node strings, i think that's
#     # a good idea.
#     node=result.signal.nodes[2:end],
#     inputmode = result.signal.modes,
#     inputport = result.signal.ports,
#     freqindex=1:length(result.signal.w),
# )


# end

# function fieldstoarray(nodeflux, inputmodes, outputmodes,
#     inputportnumbers, outputportnumbers)


# end
