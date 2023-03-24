
"""
    hbsolveold(ws, wp, Ip, Nsignalmodes, Npumpmodes, circuit, circuitdefs;
        pumpports = [1], iterations = 1000, ftol = 1e-8,
        symfreqvar = nothing, nbatches = Base.Threads.nthreads(), sorting = :number,
        returnS = true, returnSnoise = false, returnQE = true, returnCM = true,
        returnnodeflux = false, returnvoltage = false)

Harmonic balance solver for single-pump four wave mixing processes in circuits
containing Josephson junctions, capacitors, inductors, and resistors. Dissipation
can be included through frequency dependent resistors or complex capacitance.

Returns user specified scattering parameters, quantum efficiency, and node
fluxes or voltages.

# Arguments
- `ws`: signal frequency or vector of signal frequencies in radians/second.
- `wp`: pump frequency in radians/second. This function only supports a single
    pump frequency.
- `Ip`: pump current or vector of pump currents in amps. Length of `Ip` must
    be equal to length of `pumpports`.
- `Nsignalmodes`: number of signal and idler modes.
- `Npumpmodes`: number of pump modes (pump harmonics).
- `circuit`: vector of tuples containing component names, nodes, and values.
- `circuitdefs`: dictionary defining the numerical values of circuit components.

# Keywords
- `pumpports = [1]`: vector of pump port numbers. Default is a single pump at port 1.
- `iterations = 1000`: number of iterations at which the nonlinear solver stops
    even if convergence criteria not reached. 
- `ftol = 1e-8`: relative or absolute tolerance at which nonlinear solver stops
    (whichever is reached first).
- `symfreqvar = nothing`: symbolic frequency variable which is set to `nothing`
    by default but should be set equal to the frequency variable like `w` if 
    there is frequency dependence.
- `nbatches = Base.Threads.nthreads()`: for the linearized harmonic balance solution,
    split the solutions for different frequencies into this many batches. Set
    equalt to the number of threads. Recommend configuring Julia to use 
    Sys.CPU_THREADS/2 threads. 
- `sorting = :number`: sort the ports by turning them into integers and sorting
    those integers. See [`sortnodes`](@ref) for other options if this fails.
- `returnS = true`: if `true`, return the scattering parameters for each set
    of ports and signal and idler frequencies.
- `returnSnoise = false`: if `true`, return the scattering parameters corresponding
    to inputs at the noise ports (lossy components) and outputs at the physical
    ports for the signal and idler frequencies. 
- `returnQE = true`: if `true`, return the quantum efficiency for each signal
    and idler at each combinaton of ports.
- `returnCM = true`: if `true`, return the commutation relations for each signal
    and idler at each combinaton of ports (should equal ±1).
- `returnnodeflux = false`: if `true`, return the node fluxes for each signal
    and idler at each node. Set to `false` by default to reduce memory usage.
- `returnnodefluxadjoint = false`: if `true`, return the node fluxes adjoint
    for each signal and idler at each node. Set to `false` by default to
    reduce memory usage.
- `returnvoltage = false`: if `true`, return the node voltages for each signal
    and idler at each node. Set to `false` by default to reduce memory usage.

# Examples
```
@variables Rleft Cc Lj Cj w L1
circuit = Tuple{String,String,String,Num}[]
push!(circuit,("P1","1","0",1))
push!(circuit,("R1","1","0",Rleft))
push!(circuit,("C1","1","2",Cc)) 
push!(circuit,("Lj1","2","0",Lj)) 
push!(circuit,("C2","2","0",Cj))
circuitdefs = Dict(
    Lj =>1000.0e-12,
    Cc => 100.0e-15,
    Cj => 1000.0e-15,
    Rleft => 50.0,
)
ws = 2*pi*(4.5:0.01:5.0)*1e9
wp = 2*pi*4.75001*1e9
Ip = 0.00565e-6
Nsignalmodes = 8
Npumpmodes = 8
result=JosephsonCircuits.hbsolveold(ws, wp, Ip, Nsignalmodes, Npumpmodes, circuit, circuitdefs,pumpports=[1])
using Plots;plot(ws/(2*pi*1e9),10*log10.(abs2.(result.signal.S[result.signal.signalindex,result.signal.signalindex,:])))
```
"""
function hbsolveold(ws, wp, Ip, Nsignalmodes, Npumpmodes, circuit, circuitdefs;
    pumpports = [1], iterations = 1000, ftol = 1e-8,
    symfreqvar = nothing, nbatches = Base.Threads.nthreads(), sorting = :number,
    returnS = true, returnSnoise = false, returnQE = true, returnCM = true,
    returnnodeflux = false,returnnodefluxadjoint = false, returnvoltage = false,
    returnvoltageadjoint = false, keyedarrays::Val{K} = Val(false)) where K

    # parse and sort the circuit
    psc = parsesortcircuit(circuit, sorting=sorting)

    # calculate the circuit graph
    cg = calccircuitgraph(psc)

    # solve the nonlinear system
    pump=hbnlsolveold(wp, Ip, Npumpmodes, psc, cg, circuitdefs, ports = pumpports,
        iterations = iterations, ftol = ftol,
        symfreqvar = symfreqvar)

    # the node flux
    # phin = pump.out.zero
    nodeflux = pump.nodeflux

    # convert from node flux to branch flux
    phib = pump.Rbnm*nodeflux

    # calculate the sine and cosine nonlinearities from the pump flux
    Am = sincosnloddtoboth(phib[pump.Ljbm.nzind], length(pump.Ljb.nzind), pump.Nmodes)

    # solve the linear system
    signal=hblinsolveold(ws, psc, cg, circuitdefs, wp = wp, Nmodes = Nsignalmodes,
        Am = Am, symfreqvar = symfreqvar, nbatches = nbatches,
        returnS = returnS, returnSnoise = returnSnoise, returnQE = returnQE,
        returnCM = returnCM, returnnodeflux = returnnodeflux,
        returnnodefluxadjoint = returnnodefluxadjoint,
        returnvoltage = returnvoltage, returnvoltageadjoint = returnvoltageadjoint,
        keyedarrays = keyedarrays)
    return HB(pump, signal)
end

"""
    hblinsolveold(w, circuit, circuitdefs; wp = 0.0, Nmodes = 1,
        Am = zeros(Complex{Float64},0,0), symfreqvar = nothing,
        nbatches = Base.Threads.nthreads(), sorting = :number, returnS = true,
        returnSnoise = false, returnQE = true, returnCM = true,
        returnnodeflux = false, returnvoltage = false)

Linearized harmonic balance solver for single-pump four wave mixing processes in circuits
containing Josephson junctions, capacitors, inductors, and resistors. Dissipation
can be included through frequency dependent resistors or complex capacitance.

Returns user specified scattering parameters, quantum efficiency, and node
fluxes or voltages.

# Arguments
- `w`: signal frequency or vector of signal frequencies in radians/second.
- `circuit`: vector of tuples containing component names, nodes, and values. 
- `circuitdefs`: dictionary defining the numerical values of circuit components.

# Keywords
- `wp = 0.0`: pump frequency in radians/second. This function only supports a single
    pump frequency.
- `Nmodes = 1`: number of signal and idler modes.
- `Am = zeros(Complex{Float64},0,0)`: 
- `symfreqvar = nothing`: symbolic frequency variable which is set to `nothing`
    by default but should be set equal to the frequency variable like `w` if 
    there is frequency dependence.
- `nbatches = Base.Threads.nthreads()`: for the linearized harmonic balance solution,
    split the solutions for different frequencies into this many batches. Set
    equalt to the number of threads. Recommend configuring Julia to use 
    Sys.CPU_THREADS/2 threads. 
- `sorting = :number`: sort the ports by turning them into integers and sorting
    those integers. See [`sortnodes`](@ref) for other options if this fails.
- `returnS = true`: if `true`, return the scattering parameters for each set
    of ports and signal and idler frequencies.
- `returnSnoise = false`: if `true`, return the scattering parameters corresponding
    to inputs at the noise ports (lossy components) and outputs at the physical
    ports for the signal and idler frequencies. 
- `returnQE = true`: if `true`, return the quantum efficiency for each signal
    and idler at each combinaton of ports.
- `returnCM = true`: if `true`, return the commutation relations for each signal
    and idler at each combinaton of ports (should equal ±1).
- `returnnodeflux = false`: if `true`, return the node fluxes for each signal
    and idler at each node. Set to `false` by default to reduce memory usage.
- `returnnodefluxadjoint = false`: if `true`, return the node fluxes adjoint
    for each signal and idler at each node. Set to `false` by default to
    reduce memory usage.
- `returnvoltage = false`: if `true`, return the node voltages for each signal
    and idler at each node. Set to `false` by default to reduce memory usage. 

# Examples
```
@variables Rleft Cc Lj Cj w L1
circuit = Tuple{String,String,String,Num}[]
push!(circuit,("P1","1","0",1))
push!(circuit,("R1","1","0",Rleft))
push!(circuit,("C1","1","2",Cc)) 
push!(circuit,("Lj1","2","0",Lj)) 
push!(circuit,("C2","2","0",Cj))
circuitdefs = Dict(
    Lj =>1000.0e-12,
    Cc => 100.0e-15,
    Cj => 1000.0e-15,
    Rleft => 50.0,
)
w = 2*pi*(4.5:0.01:5.0)*1e9
result=JosephsonCircuits.hblinsolveold(w, circuit, circuitdefs)
using Plots;plot(w/(2*pi*1e9),angle.(result.S[:]))
```
"""
function hblinsolveold(w, circuit, circuitdefs; wp = 0.0, Nmodes = 1,
    Am = zeros(Complex{Float64},0,0), symfreqvar = nothing,
    nbatches = Base.Threads.nthreads(), sorting = :number, returnS = true,
    returnSnoise = false, returnQE = true, returnCM = true,
    returnnodeflux = false, returnnodefluxadjoint = false,
    returnvoltage = false, returnvoltageadjoint = false,
    keyedarrays::Val{K} = Val(false)) where K

    # parse and sort the circuit
    psc = parsesortcircuit(circuit, sorting = sorting)

    # calculate the circuit graph
    cg = calccircuitgraph(psc)

    return hblinsolveold(w, psc, cg, circuitdefs; wp = wp, Nmodes = Nmodes,
        Am = Am, symfreqvar = symfreqvar, nbatches = nbatches, 
        returnS = returnS, returnSnoise = returnSnoise, returnQE = returnQE,
        returnCM = returnCM, returnnodeflux = returnnodeflux,
        returnnodefluxadjoint = returnnodefluxadjoint,
        returnvoltage = returnvoltage, returnvoltageadjoint = false,
        keyedarrays = keyedarrays)
end

function hblinsolveold(w, psc::ParsedSortedCircuit, cg::CircuitGraph,
    circuitdefs; wp = 0.0, Nmodes::Integer = 1, Am = zeros(Complex{Float64},0,0),
    symfreqvar = nothing, nbatches::Integer = Base.Threads.nthreads(),
    returnS = true, returnSnoise = false, returnQE = true, returnCM = true,
    returnnodeflux = false, returnnodefluxadjoint = false,
    returnvoltage = false, returnvoltageadjoint = false,
    keyedarrays::Val{K} = Val(false)) where K

    @assert nbatches >= 1

    # calculate the numeric matrices
    nm = numericmatrices(psc, cg, circuitdefs, Nmodes = Nmodes)

    # extract the elements we need
    Nnodes = psc.Nnodes
    typevector = psc.typevector
    nodeindexarraysorted = psc.nodeindexarraysorted
    Nbranches = cg.Nbranches
    edge2indexdict = cg.edge2indexdict
    Ljb = nm.Ljb
    Rbnm = nm.Rbnm
    Cnm = nm.Cnm
    Gnm = nm.Gnm
    invLnm = nm.invLnm
    portindices = nm.portindices
    portnumbers = nm.portnumbers
    portimpedanceindices = nm.portimpedanceindices
    vvn = nm.vvn

    # if there is a symbolic frequency variable, then we need to redo the noise
    # port calculation because calcnoiseportimpedanceindices() can't tell if a
    # symbolic expression is complex. 
    if isnothing(symfreqvar)
        noiseportimpedanceindices = nm.noiseportimpedanceindices
    else
        noiseportimpedanceindices = calcnoiseportimpedanceindices(
            psc.typevector, psc.nodeindexarraysorted,
            psc.mutualinductorvector,
            Symbolics.substitute.(nm.vvn, symfreqvar => w[1]))
    end

    # generate the mode indices and find the signal index
    indices = calcindices(Nmodes)
    signalindex = findall(indices .== 0)[1]

    # if Am is undefined, define it.
    if length(Am) == 0
        Am = zeros(Complex{Float64}, 2*Nmodes, length(Ljb.nzval))

        Am[1,:] .= 1.0
    end

    # calculate AoLjbm, the sparse branch AoLj matrix
    AoLjbm = calcAoLjbm(Am,Ljb,1,Nmodes,Nbranches)

    # convert to a sparse node matrix
    AoLjnm = transpose(Rbnm)*AoLjbm*Rbnm

    # calculate the source currents
    Nports = length(portindices)

    # calculate the source terms in the branch basis
    bbm = zeros(Complex{Float64},Nbranches*Nmodes,Nmodes*Nports)

    # add a current source for each port and mode
    for (i,val) in enumerate(portindices)
        key = (nodeindexarraysorted[1,val],nodeindexarraysorted[2,val])
        for j = 1:Nmodes
            bbm[(edge2indexdict[key]-1)*Nmodes+j,(i-1)*Nmodes+j] = 1
            # bbm2[(i-1)*Nmodes+j,(i-1)*Nmodes+j] = Lmean*1/phi0
        end
    end

    # calculate the source terms in the node basis
    bnm = transpose(Rbnm)*bbm

    # find the indices at which there are symbolic variables so we can
    # perform a substitution on only those. 
    Cnmfreqsubstindices  = symbolicindices(Cnm)
    Gnmfreqsubstindices  = symbolicindices(Gnm)
    invLnmfreqsubstindices  = symbolicindices(invLnm)

    # drop any zeros in AoLjnm
    dropzeros!(AoLjnm)
    # make something with the same shape as A so we can find the sparsity
    # pattern. substitute in the symbolic variables for this addition. 
    wmodes = calcw(w[1],indices,wp);
    Cnmcopy = freqsubst(Cnm,wmodes,symfreqvar)
    Gnmcopy = freqsubst(Gnm,wmodes,symfreqvar)
    invLnmcopy = freqsubst(invLnm,wmodes,symfreqvar)
    # Asparse = (AoLjnm + invLnmcopy + Gnmcopy + Cnmcopy)
    Asparse = spaddkeepzeros(spaddkeepzeros(spaddkeepzeros(AoLjnm,invLnmcopy),Gnmcopy),Cnmcopy)

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
        S = zeros(Complex{Float64},Nports*Nmodes,Nports*Nmodes,length(w))
    else
        S = zeros(Complex{Float64},0,0,0)
    end

    if returnSnoise
        Snoise = zeros(Complex{Float64},length(noiseportimpedanceindices)*Nmodes,Nports*Nmodes,length(w))
    else
        Snoise = zeros(Complex{Float64},0,0,0)
    end

    if returnQE
        QE = zeros(Float64,Nports*Nmodes,Nports*Nmodes,length(w))
    else
        QE = zeros(Float64,0,0,0)
    end

    if returnCM
        CM = zeros(Float64,Nports*Nmodes,length(w))
    else
        CM = zeros(Float64,0,0)
    end

    if returnnodeflux
        nodeflux = zeros(Complex{Float64},Nmodes*(Nnodes-1),Nmodes*Nports,length(w))
    else
        nodeflux = Vector{Complex{Float64}}(undef,0)
    end

    if returnnodefluxadjoint
        nodefluxadjoint = zeros(Complex{Float64},Nmodes*(Nnodes-1),Nmodes*Nports,length(w))
    else
        nodefluxadjoint = Vector{Complex{Float64}}(undef,0)
    end

    if returnvoltage
        voltage = zeros(Complex{Float64},Nmodes*(Nnodes-1),Nmodes*Nports,length(w))
    else
        voltage = Vector{Complex{Float64}}(undef,0)
    end

    if returnvoltageadjoint
        voltageadjoint = zeros(Complex{Float64},Nmodes*(Nnodes-1),Nmodes*Nports,length(w))
    else
        voltageadjoint = Vector{Complex{Float64}}(undef,0)
    end

    wpumpmodes = calcw(0.0,indices,wp);
    modes = [(2*i,) for i in indices]


    # solve the linear system for the specified frequencies. the response for
    # each frequency is independent so it can be done in parallel; however
    # we want to reuse the factorization object and other input arrays. 
    # perform array allocations and factorization "nbatches" times.
    # parallelize using native threading
    batches = collect(Base.Iterators.partition(1:length(w),1+(length(w)-1)÷nbatches))
    Base.Threads.@threads for i in 1:length(batches)
        hblinsolve_inner!(S, Snoise, QE, CM, nodeflux, nodefluxadjoint, voltage,
            voltageadjoint, Asparse, 
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
        Sout = Stokeyed(S, modes, portnumbers, modes,
            portnumbers, w)
    else
        Sout = S
    end

    if returnQE && K
        QEout = Stokeyed(QE, modes, portnumbers, modes,
            portnumbers, w)
        QEidealout = Stokeyed(QEideal, modes, portnumbers, modes,
            portnumbers, w)
    else
        QEout = QE
        QEidealout = QEideal
    end

    if returnCM && K
        CMout = CMtokeyed(CM, modes, portnumbers, w)
    else
        CMout = CM
    end

    return LinearHB(Sout, Snoise, QEout, QEidealout, CMout, nodeflux, nodefluxadjoint,
        voltage, voltageadjoint, Nmodes, Nnodes, Nbranches,
        psc.uniquenodevectorsorted[2:end], portnumbers, signalindex, w, modes)
end



"""
    hbnlsolveold(wp, Ip, Nmodes, circuit, circuitdefs; ports = [1],
        iterations = 1000, ftol = 1e-8, symfreqvar = nothing,
        sorting = :number)

Nonlinear harmonic balance solver for single-pump four wave mixing processes in circuits
containing Josephson junctions, capacitors, inductors, and resistors. Dissipation
can be included through frequency dependent resistors or complex capacitance.

# Arguments
- `wp`: pump frequency in radians/second. This function only supports a single
    pump frequency.
- `Ip`: pump current or vector of pump currents in amps. Length of `Ip` must
    be equal to length of `ports`.
- `Nmodes`: number of modes (harmonics).
- `circuit`: vector of tuples containing component names, nodes, and values.
- `circuitdefs`: dictionary defining the numerical values of circuit components.

# Keywords
- `ports = [1]`: vector of drive port numbers. Default is a single drive at port 1.
- `iterations = 1000`: number of iterations at which the nonlinear solver stops
    even if convergence criteria not reached.
- `ftol = 1e-8`: relative or absolute tolerance at which nonlinear solver stops
    (whichever is reached first).
- `symfreqvar = nothing`: symbolic frequency variable which is set to `nothing`
    by default but should be set equal to the frequency variable like `w` if
    there is frequency dependence.
- `sorting = :number`: sort the ports by turning them into integers and sorting
    those integers. See [`sortnodes`](@ref) for other options if this fails.

# Examples
```
@variables Rleft Cc Lj Cj w L1
circuit = Tuple{String,String,String,Num}[]
push!(circuit,("P1","1","0",1))
push!(circuit,("R1","1","0",Rleft))
push!(circuit,("C1","1","2",Cc)) 
push!(circuit,("Lj1","2","0",Lj)) 
push!(circuit,("C2","2","0",Cj))
circuitdefs = Dict(
    Lj =>1000.0e-12,
    Cc => 100.0e-15,
    Cj => 1000.0e-15,
    Rleft => 50.0,
)
wp = 2*pi*4.75001*1e9
Ip = 0.00565e-6
Nmodes = 8
hbnlsolve(wp, Ip, Nmodes, circuit, circuitdefs, ports=[1])
```
"""
function hbnlsolveold(wp, Ip, Nmodes, circuit, circuitdefs; ports = [1],
    iterations = 1000, ftol = 1e-8,
    switchofflinesearchtol = 1e-5, alphamin = 1e-4, symfreqvar = nothing,
    sorting = :number, keyedarrays::Val{K} = Val(false)) where K

    # parse and sort the circuit
    psc = parsesortcircuit(circuit, sorting = sorting)

    # calculate the circuit graph
    cg = calccircuitgraph(psc)

    return hbnlsolveold(wp, Ip, Nmodes, psc, cg, circuitdefs; ports = ports,
        iterations = iterations, ftol = ftol,
        switchofflinesearchtol = switchofflinesearchtol, alphamin = alphamin,
        symfreqvar = symfreqvar, keyedarrays = keyedarrays)

end

function hbnlsolveold(wp, Ip, Nmodes, psc::ParsedSortedCircuit, cg::CircuitGraph,
    circuitdefs; ports = [1], iterations = 1000, ftol = 1e-8,
    switchofflinesearchtol = 1e-5, alphamin = 1e-4,
    symfreqvar = nothing, keyedarrays::Val{K} = Val(false)) where K

    if length(ports) != length(Ip)
        throw(DimensionMismatch("Number of currents Ip must be equal to number of pump ports"))
    end

    # calculate the numeric matrices
    nm = numericmatrices(psc, cg, circuitdefs, Nmodes = Nmodes)

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
    noiseportimpedanceindices = nm.noiseportimpedanceindices
    vvn = nm.vvn
    Lmean = nm.Lmean
    # if there are no inductors, then Lmean will be zero so set it to be one
    if iszero(Lmean)
        Lmean = one(eltype(Lmean))
    end
    Lb = nm.Lb

    # make a sparse transpose (improves multiplication speed slightly)
    Rbnmt = sparse(transpose(Rbnm))

    # generate the pump frequencies
    # odd harmonics of the pump
    wmodes = zeros(Float64,Nmodes)
    for i = 1:Nmodes
        wmodes[i]=wp*(1 + 2*(i-1))
    end
    wmodesm = Diagonal(repeat(wmodes,outer=Nnodes-1))
    wmodes2m = Diagonal(repeat(wmodes.^2,outer=Nnodes-1))

    # calculate the source terms in the branch basis
    bbm = zeros(Complex{Float64},Nbranches*Nmodes)
    for (i,port) in enumerate(ports)
        for (j,portnumber) in enumerate(portnumbers)
            if portnumber == port
                portindex = portindices[j]
                key = (nodeindexarraysorted[1,portindex],nodeindexarraysorted[2,portindex])
                bbm[(edge2indexdict[key]-1)*Nmodes+1] = Lmean*Ip[i]/phi0
            end
        end
    end

    # convert from the node basis to the branch basis
    bnm = transpose(Rbnm)*bbm

    # define a random initial value for Am so we can get the sparsity structure of
    # the Jacobian. 
    Am = rand(Complex{Float64},2*Nmodes,length(Ljb.nzval))

    # calculate AoLjbm, the sparse branch AoLj matrix
    AoLjbm = calcAoLjbm(Am,Ljb,Lmean,Nmodes,Nbranches)
    AoLjbmcopy = calcAoLjbm(Am,Ljb,Lmean,Nmodes,Nbranches)

    # convert to a sparse node matrix. Note: I was having problems with type 
    # instability when i used AoLjbm here instead of AoLjbmcopy. 
    # AoLjnm = transpose(Rbnm)*AoLjbmcopy*Rbnm;
    AoLjbmRbnm = AoLjbmcopy*Rbnm
    xbAoLjbmRbnm = fill(false, size(AoLjbmcopy,1))
    AoLjnm = Rbnmt*AoLjbmRbnm
    xbAoLjnm = fill(false, size(Rbnmt,1))

    # define arrays of zeros for the function
    x = zeros(Complex{Float64},(Nnodes-1)*Nmodes)
    F = zeros(Complex{Float64},(Nnodes-1)*Nmodes)
    AoLjbmvector = zeros(Complex{Float64},Nbranches*Nmodes)

    # substitute in the symbolic variables
    Cnm = freqsubst(Cnm,wmodes,symfreqvar)
    Gnm = freqsubst(Gnm,wmodes,symfreqvar)
    invLnm = freqsubst(invLnm,wmodes,symfreqvar)

    # scale the matrices for numerical reasons
    Cnm *= Lmean
    Gnm *= Lmean
    invLnm *= Lmean

    # calculate the structure of the Jacobian
    # Jsparse = AoLjnm + invLnm - im*Gnm*wmodesm - Cnm*wmodes2m
    Jsparse = spaddkeepzeros(spaddkeepzeros(spaddkeepzeros(AoLjnm,invLnm),Gnm),Cnm)


    # make the index maps so we can efficiently add the sparse matrices 
    # together
    AoLjnmindexmap = sparseaddmap(Jsparse,AoLjnm)
    invLnmindexmap = sparseaddmap(Jsparse,invLnm)
    Gnmindexmap = sparseaddmap(Jsparse,Gnm)
    Cnmindexmap = sparseaddmap(Jsparse,Cnm)

    # build the function and Jacobian for solving the nonlinear system
    function fj!(F, Jsparse, x)
        calcfj!(F,Jsparse,x,wmodesm,wmodes2m,Rbnm,Rbnmt,invLnm,Cnm,Gnm,bnm,
            Ljb,Ljbm,Nmodes,
            Nbranches,Lmean,AoLjbmvector,AoLjbm,
            AoLjnmindexmap,invLnmindexmap,Gnmindexmap,Cnmindexmap,
            AoLjnm, xbAoLjnm, AoLjbmRbnm, xbAoLjbmRbnm)
        return nothing
    end

    # # use this for debugging purposes to return the function and the
    # # Jacobian
    # fj!(F,Jsparse,x)
    # return (F,Jsparse,x)

    # fj!(F,Jsparse,x)
    # @show bnm
    # @show Lmean
    # @show Jsparse

    # solve the nonlinear system. skip the nonlinear solve if the Jacobian
    # has no nonzero terms
    nlsolve!(fj!, F, Jsparse, x; iterations = iterations, ftol = ftol,
        switchofflinesearchtol = switchofflinesearchtol, alphamin = alphamin)

    nodeflux = x

    # calculate the scattering parameters for the pump
    Nports = length(portindices)
    S = zeros(Complex{Float64}, Nports*Nmodes, Nports*Nmodes)
    inputwave = zeros(Complex{Float64}, Nports*Nmodes)
    outputwave = zeros(Complex{Float64},Nports*Nmodes)
    portimpedances = [vvn[i] for i in portimpedanceindices]
    if !isempty(S)
        calcinputoutput!(inputwave,outputwave,nodeflux,bnm/Lmean,portimpedanceindices,portimpedanceindices,
            portimpedances,portimpedances,nodeindexarraysorted,typevector,wmodes,symfreqvar)
        calcscatteringmatrix!(S,inputwave,outputwave)
    end

    # calculate the frequency struct which we use in v2 of the solvers
    freq = removeconjfreqs(
        truncfreqs(
            calcfreqsrdft((2*Nmodes,)),
            dc=false, odd=true, even=false, maxintermodorder=Inf,
        )
    )
    modes = freq.modes

    return NonlinearHB((wp,), freq, nodeflux, Rbnm, Ljb, Lb, Ljbm, Nmodes,
        Nbranches, psc.uniquenodevectorsorted[2:end], portnumbers, modes, S)
end

"""
    calcfj(F,J,nodeflux,wmodesm,wmodes2m,Rbnm,invLnm,Cnm,Gnm,bm,Ljb,Ljbindices,
        Ljbindicesm,Nmodes,Lmean,AoLjbm)

Calculate the residual and the Jacobian. These are calculated with one function
in order to reuse the time domain nonlinearity calculation.

Leave off the type signatures on F and J because the solver will pass a type of
Nothing if it only wants to calculate F or J. 
"""
function calcfj!(F,
        J,
        nodeflux::AbstractVector,
        wmodesm::AbstractMatrix,
        wmodes2m::AbstractMatrix,
        Rbnm::AbstractArray{Int,2},
        Rbnmt::AbstractArray{Int,2},
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
        AoLjnm, xbAoLjnm, AoLjbmRbnm, xbAoLjbmRbnm)

    # convert from a node flux to a branch flux
    phib = Rbnm*nodeflux

    # calculate the sine and cosine nonlinearities 
    Am = sincosnloddtoboth(phib[Ljbm.nzind],nnz(Ljb),Nmodes)

    # calculate the residual
    if !isnothing(F)

        # calculate the function. use the sine terms. Am[2:2:2*Nmodes,:]
        # calculate  AoLjbm, this is just a diagonal matrix.
        # check if this is consistent with other calculations
        @inbounds for i in 1:nnz(Ljb)
            for j in 1:Nmodes
                AoLjbmvector[(Ljb.nzind[i]-1)*Nmodes + j] = Am[2*j,i]*(Lmean/Ljb.nzval[i])
            end
        end

        F .= Rbnmt*AoLjbmvector .+ invLnm*nodeflux .+ im*Gnm*wmodesm*nodeflux .- Cnm*wmodes2m*nodeflux .- bnm
    end

    # calculate the Jacobian
    if !isnothing(J)

        # calculate  AoLjbm
        updateAoLjbm!(AoLjbm,Am,Ljb,Lmean,Nmodes,Nbranches)

        # convert to a sparse node matrix
        # AoLjnm = Rbnmt*AoLjbm*Rbnm

        spmatmul!(AoLjbmRbnm, AoLjbm, Rbnm, xbAoLjbmRbnm)
        spmatmul!(AoLjnm, Rbnmt, AoLjbmRbnm, xbAoLjnm)

        # calculate the sparse Jacobian.
        # @time J .= AoLjnm .+ invLnm .+ im.*Gnm*wmodesm .- Cnm*wmodes2m
        # the code below adds the sparse matrices together with minimal
        # memory allocations and without changing the sparsity structure.
        fill!(J,0)
        sparseadd!(J,AoLjnm,AoLjnmindexmap)
        sparseadd!(J,invLnm,invLnmindexmap)
        sparseadd!(J,im,Gnm,wmodesm,Gnmindexmap)
        sparseadd!(J,-1,Cnm,wmodes2m,Cnmindexmap)
    end
    return nothing
end


"""
    calcAoLjbm(Am, Ljb::SparseVector, Lmean, Nmodes, Nbranches)

# Examples
```jldoctest
julia> @variables Lj1 Lj2 A11 A12 A21 A22 A31 A32;JosephsonCircuits.calcAoLjbm([A11;A21;A31],JosephsonCircuits.SparseArrays.sparsevec([1],[Lj1]),1,2,1)
2×2 SparseArrays.SparseMatrixCSC{Num, Int64} with 4 stored entries:
 A11 / Lj1  A31 / Lj1
 A31 / Lj1  A11 / Lj1

julia> @syms Lj1 Lj2 A11 A12 A21 A22 A31 A32;JosephsonCircuits.calcAoLjbm([A11;A21;A31],JosephsonCircuits.SparseArrays.sparsevec([1],[Lj1]),1,2,1).nzval
4-element Vector{Any}:
 A11 / Lj1
 A31 / Lj1
 conj(A31 / Lj1)
 A11 / Lj1

julia> @variables Lj1 Lj2 A11 A12 A21 A22 A31 A32;JosephsonCircuits.calcAoLjbm([A11 A12;A21 A22;A31 A32],JosephsonCircuits.SparseArrays.sparsevec([1,2],[Lj1,Lj2]),1,2,2)
4×4 SparseArrays.SparseMatrixCSC{Num, Int64} with 8 stored entries:
 A11 / Lj1  A31 / Lj1          ⋅          ⋅
 A31 / Lj1  A11 / Lj1          ⋅          ⋅
         ⋅          ⋅  A12 / Lj2  A32 / Lj2
         ⋅          ⋅  A32 / Lj2  A12 / Lj2

julia> @variables Lj1 Lj2 A11 A12 A21 A22 A31 A32;JosephsonCircuits.calcAoLjbm([A11;A21;A31],JosephsonCircuits.SparseArrays.sparsevec([1],[Lj1]),1,3,1)
3×3 SparseArrays.SparseMatrixCSC{Num, Int64} with 9 stored entries:
 A11 / Lj1  A31 / Lj1          0
 A31 / Lj1  A11 / Lj1  A31 / Lj1
         0  A31 / Lj1  A11 / Lj1
```
"""
function calcAoLjbm(Am, Ljb::SparseVector, Lmean, Nmodes, Nbranches)

    # define empty vectors for the rows, columns, and values
    I = Vector{eltype(Ljb.nzind)}(undef,nnz(Ljb)*Nmodes^2)
    J = Vector{eltype(Ljb.nzind)}(undef,nnz(Ljb)*Nmodes^2)

    type = promote_type(eltype(Am),eltype(1 ./Ljb.nzval))

    if type <: Symbolic
        type = Any
    end

    V = Vector{type}(undef,nnz(Ljb)*Nmodes^2)


    if size(Am,2) != nnz(Ljb)
        throw(DimensionMismatch("The second axis of Am must equal the number of nonzero elements in Ljb (the number of JJs)."))
    end

    if length(Ljb) > Nbranches
        throw(DimensionMismatch("The length of Ljb cannot be larger than the number of branches."))
    end

    # calculate  AoLjbm
    n = 1
    for i = 1:nnz(Ljb)
        for j = 1:Nmodes
            for k = 1:Nmodes

                # calculate the toeplitz matrices for each node 
                I[n]=j+(Ljb.nzind[i]-1)*Nmodes
                J[n]=k+(Ljb.nzind[i]-1)*Nmodes

                # assume terms we don't have pump data for are zero.
                index = 2*abs(j-k)+1
                if index > size(Am,1)
                    V[n] = 0
                else
                    V[n]=Am[index,i]*(Lmean/Ljb.nzval[i])
                end

                #take the complex conjugate of the upper half (not the diagonal)
                if j-k<0
                    V[n] = conj(V[n])
                end

                ## for debugging. calculate index, Ljb.nzind and i from I and J.
                # println("index: ",index," ", 2*abs(I[n]-J[n])+1)
                # println("Ljb.nzind: ",Ljb.nzind[i]," ", ((I[n]+J[n]-1) ÷ (2*Nmodes))+1)
                # println("i: ",i," ", searchsortedfirst(Ljb.nzind,((I[n]+J[n]-1) ÷ (2*Nmodes))+1))
                # println("j: ",j," ",I[n]-(Ljb.nzind[i]-1)*Nmodes," k: ",k," ",J[n]-(Ljb.nzind[i]-1)*Nmodes)

                n+=1
            end
        end
    end

    # assemble the sparse branch AoLj matrix
    return sparse(I,J,V,Nbranches*Nmodes,Nbranches*Nmodes)
end


"""
    updateAoLjbm!(AoLjbm::SparseMatrixCSC, Am, Ljb::SparseVector, Lmean, Nmodes,
        Nbranches)

# Examples
```jldoctest
@variables Lj1 Lj2 A11 A12 A21 A22 A31 A32;
AoLjbm = JosephsonCircuits.calcAoLjbm([A11 A12;A21 A22;A31 A32],JosephsonCircuits.SparseArrays.sparsevec([1,2],[Lj1,Lj2]),1,2,2);
AoLjbmcopy = copy(AoLjbm);
AoLjbmcopy.nzval .= 0;
JosephsonCircuits.updateAoLjbm!(AoLjbmcopy,[A11 A12;A21 A22;A31 A32],JosephsonCircuits.SparseArrays.sparsevec([1,2],[Lj1,Lj2]),1,2,2)
all(AoLjbmcopy.nzval .- AoLjbm.nzval .== 0)

# output
true
```
```jldoctest
@variables Lj1 Lj2 A11 A12 A21 A22 A31 A32;
AoLjbm = JosephsonCircuits.calcAoLjbm([A11 A12;A21 A22;A31 A32],JosephsonCircuits.SparseArrays.sparsevec([1,2],[Lj1,Lj2]),1,3,2);
AoLjbmcopy = copy(AoLjbm);
AoLjbmcopy.nzval .= 0;
JosephsonCircuits.updateAoLjbm!(AoLjbmcopy,[A11 A12;A21 A22;A31 A32],JosephsonCircuits.SparseArrays.sparsevec([1,2],[Lj1,Lj2]),1,3,2)
all(AoLjbmcopy.nzval .- AoLjbm.nzval .== 0)

# output
true
```
"""
function updateAoLjbm!(AoLjbm::SparseMatrixCSC, Am, Ljb::SparseVector, Lmean,
    Nmodes, Nbranches)

    # check that there are the right number of nonzero values. 
    # check that the dimensions are consistent with Nmode and Nbranches.

    if nnz(Ljb)*Nmodes^2 != nnz(AoLjbm)
        throw(DimensionMismatch("The number of nonzero elements in AoLjbm are not consistent with nnz(Ljb) and Nmodes."))
    end

    if size(Am,2) != nnz(Ljb)
        throw(DimensionMismatch("The second axis of Am must equal the number of nonzero elements in Ljb (the number of JJs)."))
    end

    if length(Ljb) > Nbranches
        throw(DimensionMismatch("The length of Ljb cannot be larger than the number of branches."))
    end

    # i want a vector length(Ljb) where the indices are the values Ljb.nzind
    # and the values are the indices of Ljb.nzind
    indexconvert = zeros(Int,length(Ljb))
    for (i,j) in enumerate(Ljb.nzind)
        indexconvert[j] = i
    end


    for l = 1:length(AoLjbm.colptr)-1
        for m in AoLjbm.colptr[l]:(AoLjbm.colptr[l+1]-1)

            i = indexconvert[((AoLjbm.rowval[m]+l-1) ÷ (2*Nmodes))+1]
            index = 2*abs(AoLjbm.rowval[m]-l)+1

            if index > size(Am,1)
                AoLjbm.nzval[m] = 0
            else
                AoLjbm.nzval[m]=Am[index,i]*(Lmean/Ljb.nzval[i])
            end

            #take the complex conjugate of the upper half (not the diagonal)
            if AoLjbm.rowval[m]-l<0
                AoLjbm.nzval[m] = conj(AoLjbm.nzval[m])
            end

        end
    end

    return nothing
end


"""
    sincosnloddtoboth(amodd::Array{Complex{Float64},1},Nbranches::Int,m::Int)

Applies the junction nonlinearity to a vector of branch fluxes of length Nbranches*m
where m is the number of odd pump harmonics (1w, 3w, 5w, etc). The ordering is
(mode 1, node 1), (mode 2, node 1) ... (mode 1, node 2) ... Returns even AND
odd terms in a 2d array with dimensions 2*m by Nbranches. 

# Examples
```jldoctest
julia> isapprox(JosephsonCircuits.sincosnloddtoboth([0.5+0.0im,0,0,0],1,4),ComplexF64[0.765197686557965 + 0.0im; 0.44005058574495276 + 0.0im; -0.11490348493140057 + 0.0im; -0.019563353994648498 + 0.0im; 0.0024766387010484486 + 0.0im; 0.0002497629794614272 + 0.0im; -2.084411456066653e-5 + 0.0im; -3.0046516347986037e-6 + 0.0im;;])
true

julia> isapprox(JosephsonCircuits.sincosnloddtoboth([0.02+0.0im,0,0.01+0.0im,0],2,2),ComplexF64[0.9996000399980445 + 0.0im 0.9999000024999694 + 0.0im; 0.019996000293322436 + 0.0im 0.009999500009166587 + 0.0im; -0.00019996666853328016 + 0.0im -4.999791669583568e-5 + 0.0im; -2.6664000106757513e-6 + 0.0im -3.3332500006440685e-7 + 0.0im])
true
```
"""
function sincosnloddtoboth(amodd::Array{Complex{Float64},1},
    Nbranches::Int, m::Int)

    if length(amodd) != Nbranches*m
        throw(DimensionMismatch("Length of node flux vector not consistent with number of modes number of frequencies"))
    end

    am = zeros(Complex{Float64},2*length(amodd))
    am[2:2:2*length(amodd)] .= amodd

    return sincosnl(reshape(am,2*m,Nbranches)) 
end


"""
    sincosnl(am::Array{Complex{Float64},2})

Applies the junction nonlinearity to a vector of Fourier coefficients of the
phases across the junction of size 2*m by (Nnodes-1) where m is the number of pump
harmonics (0w, 1w, 2w, 3w, etc). To save time, this calculates both the sine and
cosine nonlinearities at the same time. If the input is odd harmonics, the sine
terms will also be odd harmonics the cosine terms will be even harmonics.

# Examples
```jldoctest
julia> isapprox(JosephsonCircuits.sincosnl([0 0.001+0im;0 0]),ComplexF64[1.0 + 0.0im 1.000999499833375 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im])
true

julia> isapprox(JosephsonCircuits.sincosnl([0 0;0.001+0im 0;0 0;0 0; 0 0]), ComplexF64[0.99999900000025 + 0.0im 1.0 + 0.0im; 0.0009999995000000916 + 0.0im 0.0 + 0.0im; -4.999998333421463e-7 + 0.0im 0.0 + 0.0im; -1.6666664597909941e-10 + 0.0im 0.0 + 0.0im; 8.337774914934926e-14 + 0.0im 0.0 + 0.0im])
true

julia> isapprox(JosephsonCircuits.sincosnl([0 0;0.001+0im 0.25+0im;0 0;0 0; 0 0]),ComplexF64[0.99999900000025 + 0.0im 0.9384698079924576 + 0.0im; 0.0009999995000000916 + 0.0im 0.24226844566945333 + 0.0im; -4.999998333421463e-7 + 0.0im -0.03060435952740681 + 0.0im; -1.6666664597909941e-10 + 0.0im -0.0025556763673518224 + 0.0im; 8.337774914934926e-14 + 0.0im 0.0003214729527288296 + 0.0im])
true
```
"""
function sincosnl(am::Array{Complex{Float64},2})

    #choose the number of time points based on the number of fourier
    #coefficients
    #stepsperperiod = 2*size(am)[1]-2
    # changed to below because i wasn't using enough points when Nmodes=1.
    # the results contained only real values. 
    if size(am,1) == 2
        stepsperperiod = 2*size(am,1)-1
    else
        stepsperperiod = 2*size(am,1)-2
    end

    #transform back to time domain
    ift = FFTW.irfft(am,stepsperperiod,[1])*stepsperperiod

    #apply the nonlinear function
    nlift = cos.(ift) .+ sin.(ift)

    #fourier transform back to the frequency domain
    ftnlift = FFTW.rfft(nlift,[1])/stepsperperiod

    return ftnlift

end


"""
    calcindices(m::Integer)

The indices over which to calculate the idlers using the formula ws+2*i*wp 
where i is an index. This could be defined differently without causing any issues. 

# Examples
```jldoctest
julia> JosephsonCircuits.calcindices(7)
-3:3

julia> JosephsonCircuits.calcindices(8)
-4:3
```
"""
function calcindices(m::Integer)
    return -(m÷2):((m-1)÷2)
end

"""
    calcw(ws::Number, i::Integer, wp::Number)

Generate the signal and idler frequencies using the formula ws + 2*i*wp

Should I switch this to ws+i*wp so it can handle three wave mixing then 
always double i for four wave mixing?


# Examples
```jldoctest
julia> JosephsonCircuits.calcw(2*pi*4.0e9,-1,2*pi*5.0e9)/(2*pi*1.0e9)
-5.999999999999999

julia> JosephsonCircuits.calcw(2*pi*4.0e9,[-1,0,1],2*pi*5.0e9)/(2*pi*1.0e9)
3-element Vector{Float64}:
 -5.999999999999999
  4.0
 14.0

julia> JosephsonCircuits.calcw([2*pi*4.0e9,2*pi*4.1e9],[-1,0,1],2*pi*5.0e9)/(2*pi*1.0e9)
2×3 Matrix{Float64}:
 -6.0  4.0  14.0
 -5.9  4.1  14.1
```
"""
function calcw(ws, i::Integer, wp)
    return ws + 2*i*wp
end
function calcw(ws, i::AbstractVector, wp)

    w = zeros(typeof(ws),length(i))
    calcw!(ws,i,wp,w)
    return w
end
function calcw(ws::AbstractVector, i::AbstractVector, wp)

    w = zeros(eltype(ws),length(ws),length(i))
    calcw!(ws,i,wp,w)
    return w
end

"""
    calcw!(ws, i, wp, w)

Generate the signal and idler frequencies using the formula ws + 2*i*wp.
Overwrites w with output. 
"""
function calcw!(ws, i::AbstractVector, wp, w::AbstractVector)

    for j in 1:length(i)
        w[j] = calcw(ws,i[j],wp)
    end

    return nothing
end
function calcw!(ws::AbstractVector, i::AbstractVector, wp, w::AbstractMatrix)

    for j in 1:length(ws)
        for k in 1:length(i)
            w[j,k] = calcw(ws[j],i[k],wp)
        end
    end

    return nothing
end
