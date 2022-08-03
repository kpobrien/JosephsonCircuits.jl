
"""

    hbsolve()

"""

function hbsolve(ws,wp,Ip,Nsignalmodes,Npumpmodes,circuit,circuitdefs; pumpports=[1],
    solver =:klu, iterations=1000,ftol=1e-8,symfreqvar=nothing,
    nbatches = Base.Threads.nthreads(), sorting=:number, returnS = true,
    returnSnoise = false, returnQE = true,returnCM = true,
    returnnodeflux = false, returnvoltage = false,verbosity = 1)

    # parse and sort the circuit
    psc = parsesortcircuit(circuit,sorting=sorting)

    # calculate the circuit graph
    cg = calccircuitgraph(psc)

    # solve the nonlinear system    
    pump=hbnlsolve(wp,Ip,Npumpmodes,psc,cg,circuitdefs,ports=pumpports,
        solver=solver,iterations=iterations,ftol=ftol,symfreqvar=symfreqvar,
        verbosity = verbosity)

    # the node flux
    # phin = pump.out.zero
    phin = pump.phin    

    # convert from node flux to branch flux
    phib = pump.Rbnm*phin

    # calculate the sine and cosine nonlinearities from the pump flux
    Am = sincosnloddtoboth(phib[pump.Ljbm.nzind],length(pump.Ljb.nzind),pump.Nmodes)

    # solve the linear system
    signal=hblinsolve(ws,psc,cg,circuitdefs,wp=wp,Nmodes=Nsignalmodes,
        Am=Am,solver=solver,symfreqvar=symfreqvar,nbatches = nbatches,
        returnS = returnS, returnSnoise = returnSnoise, returnQE = returnQE,
        returnCM = returnCM,returnnodeflux = returnnodeflux,
        returnvoltage = returnvoltage,verbosity = verbosity)
    return HB(pump,Am,signal)
end

struct HB
    pump
    Am
    signal
end


"""

    hblinsolve(w,wp,Nmodes,circuit,circuitdefs,Am)

"""
function hblinsolve(w,circuit,circuitdefs;wp=0,Nmodes=1,Am=zeros(Complex{Float64},0,0),
    solver=:klu,symfreqvar=nothing, nbatches = Base.Threads.nthreads(), 
    sorting=:number,returnS = true, returnSnoise = false, returnQE = true,
    returnCM = true,returnnodeflux = false, returnvoltage = false,verbosity = 1)

    # parse and sort the circuit
    psc = parsesortcircuit(circuit,sorting=sorting)

    # calculate the circuit graph
    cg = calccircuitgraph(psc)

    return hblinsolve(w,psc,cg,circuitdefs;wp=wp,Nmodes=Nmodes,Am=Am,
        solver=solver,symfreqvar=symfreqvar, nbatches = nbatches, 
        returnS = returnS, returnSnoise = returnSnoise, returnQE = returnQE,
        returnCM = returnCM,returnnodeflux = returnnodeflux,
        returnvoltage = returnvoltage,verbosity = verbosity)
end

function hblinsolve(w,psc::ParsedSortedCircuit,cg::CircuitGraph,
    circuitdefs;wp=0.0,Nmodes::Integer=1,Am=zeros(Complex{Float64},0,0),
    solver=:klu,symfreqvar=nothing, nbatches::Integer = Base.Threads.nthreads(),
    returnS = true, returnSnoise = false, returnQE = true, returnCM = true,
    returnnodeflux = false, returnvoltage = false,verbosity = 1)

    @assert nbatches >= 1

    # calculate the numeric matrices
    nm=numericmatrices(psc,cg,circuitdefs,Nmodes=Nmodes)

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
    portimpedanceindices = nm.portimpedanceindices
    noiseportimpedanceindices = nm.noiseportimpedanceindices
    vvn = nm.vvn
    # portdict = nm.portdict
    # resistordict = nm.resistordict
    # capacitornoiseports  = nm.capacitornoiseports
    # resistornoiseports = nm.resistornoiseports

    # generate the mode indices and find the signal index
    indices = calcindices(Nmodes)
    signalindex = findall(indices .== 0)[1]

    # if Am is undefined, define it. 
    if length(Am) == 0
        Am = zeros(Complex{Float64},2*Nmodes,length(Ljb.nzval))

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
    # capacitornoiseportsfreqsubstindices  = symbolicindices(capacitornoiseports)
    # resistornoiseportsfreqsubstindices = symbolicindices(resistornoiseports)

    # drop any zeros in AoLjnm
    dropzeros!(AoLjnm)
    # make something with the same shape as A so we can find the sparsity
    # pattern. substitute in the symbolic variables for this addition. 
    wmodes = calcw(w[1],indices,wp);
    Cnmcopy = freqsubst(Cnm,wmodes,symfreqvar)
    Gnmcopy = freqsubst(Gnm,wmodes,symfreqvar)
    invLnmcopy = freqsubst(invLnm,wmodes,symfreqvar)
    Asparse = (AoLjnm + invLnmcopy + Gnmcopy + Cnmcopy)

    # make the index maps so we can efficiently add the sparse matrices
    # together without allocations or changing the sparsity structure. 
    Cnmindexmap = sparseaddmap(Asparse,Cnmcopy)
    Gnmindexmap = sparseaddmap(Asparse,Gnmcopy)
    invLnmindexmap = sparseaddmap(Asparse,invLnmcopy)
    AoLjnmindexmap = sparseaddmap(Asparse,AoLjnm)

    # # solve for node fluxes and scattering parameters
    # portbranches = collect(keys(portdict))
    # portimpedance = collect(values(resistordict))
    # noiseportbranches = collect(keys(resistornoiseports))
    # noiseportimpedance = collect(values(resistornoiseports))

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
        nodeflux = Vector{Complex{Float64}}(undef,0)
    else
        nodeflux = Vector{Complex{Float64}}(undef,0)
    end        

    if returnvoltage
        voltage = Vector{Complex{Float64}}(undef,0)
    else
        voltage = Vector{Complex{Float64}}(undef,0)    
    end    

    # solve the linear system for the specified frequencies. the response for
    # each frequency is independent so it can be done in parallel; however
    # we want to reuse the factorization object and other input arrays. 
    # perform array allocations and factorization "nbatches" times.
    # parallelize using native threading
    batches = collect(Base.Iterators.partition(1:length(w),1+(length(w)-1)÷nbatches))
    Base.Threads.@threads for i in 1:length(batches)
        hblinsolve_inner!(S,Snoise,QE,CM,nodeflux,voltage,Asparse,AoLjnm,invLnm,Cnm,Gnm,bnm,
            AoLjnmindexmap,invLnmindexmap,Cnmindexmap,Gnmindexmap,
            Cnmfreqsubstindices,Gnmfreqsubstindices,invLnmfreqsubstindices,
            portindices,portimpedanceindices,noiseportimpedanceindices,
            portimpedances,noiseportimpedances,nodeindexarraysorted,typevector,
            w,indices,wp,Nmodes,Nnodes,solver,symfreqvar,batches[i])
    end

    if returnQE
        # calculate the ideal quantum efficiency.
        QEideal = zeros(Float64,size(S))
        calcqeideal!(QEideal,S)
    else
        QEideal = zeros(Float64,0,0,0)
    end

    return LinearHB(S,Snoise,QE,QEideal,CM,nodeflux,voltage,Nmodes,Nnodes,Nbranches,signalindex,w)

end

struct LinearHB
    S
    Snoise
    QE
    QEideal
    CM
    nodeflux
    voltage
    Nmodes
    Nnodes
    Nbranches
    signalindex
    w
end


function hblinsolve_inner!(S,Snoise,QE,CM,nodeflux,voltage,Asparse,AoLjnm,invLnm,Cnm,Gnm,bnm,
    AoLjnmindexmap,invLnmindexmap,Cnmindexmap,Gnmindexmap,
    Cnmfreqsubstindices,Gnmfreqsubstindices,invLnmfreqsubstindices,
    portindices,portimpedanceindices,noiseportimpedanceindices,
    portimpedances,noiseportimpedances,nodeindexarraysorted,typevector,
    w,indices,wp,Nmodes,Nnodes,solver,symfreqvar,wi)

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

    if solver ==:klu
        # if using the KLU factorization and sparse solver then make a 
        # factorization for the sparsity pattern. 
        F = KLU.klu(Asparsecopy)
    else
        error("Error: Unknown solver")
    end

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
        wmodes = calcw(ws,indices,wp);
        wmodesm = Diagonal(repeat(wmodes,outer=Nnodes-1));
        wmodes2m = Diagonal(repeat(wmodes.^2,outer=Nnodes-1));

        # perform the operation below in a way that doesn't allocate significant
        # memory, plus take the conjugates mentioned below. 
        # A = (AoLjnm + invLnm - im.*Gnm*wmodesm - Cnm*wmodes2m)

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

        # solve the linear system
        if solver == :klu
            try 
                # update the factorization. the sparsity structure does 
                # not change so we can reuse the factorization object.
                KLU.klu!(F,Asparsecopy)

                # solve the linear system
                ldiv!(phin,F,bnm)
            catch e
                if isa(e, SingularException)
                    # reusing the symbolic factorization can sometimes
                    # lead to numerical problems. if the first linear
                    # solve fails try factoring and solving again
                    F = KLU.klu(Asparsecopy)
                    ldiv!(phin,F,bnm)
                else
                    throw(e)
                end
            end
        else
            error("Error: Unknown solver")
        end

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

        if Nnoiseports > 0 && (!isempty(Snoise) || !isempty(QE) || !isempty(CM))

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

            # solve the linear system
            if solver == :klu
                try 
                    # update the factorization. the sparsity structure does 
                    # not change so we can reuse the factorization object.
                    KLU.klu!(F,Asparsecopy)

                    # solve the linear system
                    ldiv!(phin,F,bnm)
                catch e
                    if isa(e, SingularException)
                        # reusing the symbolic factorization can sometimes
                        # lead to numerical problems. if the first linear
                        # solve fails try factoring and solving again
                        F = KLU.klu(Asparsecopy)
                        ldiv!(phin,F,bnm)
                    else
                        throw(e)
                    end
                end
            else
                error("Error: Unknown solver")
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
    hbnlsolve(wp,Ip,Nmodes,circuit,circuitdefs)

I should make it so i can easily sweep pump frequency and pump current
without parsing this again. at the very least i should define this
so that i can pass in the already parsed circuit and already generated
numeric matrices.

"""
function hbnlsolve(wp,Ip,Nmodes,circuit,circuitdefs;ports=[1],solver=:klu,
    iterations=1000,symfreqvar=nothing,ftol=1e-8,sorting=:number,
    verbosity=1)

    # parse and sort the circuit
    psc = parsesortcircuit(circuit,sorting=sorting)

    # calculate the circuit graph
    cg = calccircuitgraph(psc)

    return hbnlsolve(wp,Ip,Nmodes,psc,cg,circuitdefs;ports=ports,
        solver=solver,iterations=iterations,ftol=ftol,symfreqvar=symfreqvar,
        verbosity=verbosity)

end

function hbnlsolve(wp,Ip,Nmodes,psc::ParsedSortedCircuit,cg::CircuitGraph,
    circuitdefs;ports=[1],solver=:klu,iterations=1000,ftol = 1e-8,symfreqvar=nothing,
    verbosity=1)


    if length(ports) != length(Ip)
        error("Number of currents Ip must be equal to number of pump ports")
    end

    # calculate the numeric matrices
    nm=numericmatrices(psc,cg,circuitdefs,Nmodes=Nmodes)

    # extract the elements we need
    Nnodes = psc.Nnodes
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
    Lb = nm.Lb

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
    for (i,portindex) in enumerate(portindices)
        portnumber = portnumbers[i]
        key = (nodeindexarraysorted[1,portindex],nodeindexarraysorted[2,portindex])
        if portnumber in ports
            # bbm[(edge2indexdict[key]-1)*Nmodes+1] = Lmean*Ip/phi0
            bbm[(edge2indexdict[key]-1)*Nmodes+1] = Lmean*Ip[portindex]/phi0
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
    AoLjnm = transpose(Rbnm)*AoLjbmcopy*Rbnm;

    # define arrays of zeros for the function
    x = zeros(Complex{Float64},(Nnodes-1)*Nmodes)
    F = zeros(Complex{Float64},(Nnodes-1)*Nmodes)
    AoLjbmvector = zeros(Complex{Float64},Nbranches*Nmodes)

    # make a sparse transpose (improves multiplication speed slightly)
    Rbnmt = sparse(transpose(Rbnm))

    # substitute in the symbolic variables
    Cnm = freqsubst(Cnm,wmodes,symfreqvar)
    Gnm = freqsubst(Gnm,wmodes,symfreqvar)
    invLnm = freqsubst(invLnm,wmodes,symfreqvar)

    # scale the matrices for numerical reasons
    Cnm *= Lmean
    Gnm *= Lmean
    invLnm *= Lmean

    # calculate the structure of the Jacobian
    Jsparse = (AoLjnm + invLnm - im*Gnm*wmodesm - Cnm*wmodes2m)

    # make the index maps so we can efficiently add the sparse matrices 
    # together
    AoLjnmindexmap = sparseaddmap(Jsparse,AoLjnm)
    invLnmindexmap = sparseaddmap(Jsparse,invLnm)
    Gnmindexmap = sparseaddmap(Jsparse,Gnm)
    Cnmindexmap = sparseaddmap(Jsparse,Cnm)

    # perform a factorization. this will be updated later for each 
    # interation
    factorization = KLU.klu(Jsparse)

    deltax = copy(x)

    Nsamples = 100
    samples = Float64[]
    fmin = Float64[]
    fvals = Float64[]
    fpvals = Float64[]
    dfdalphavals = Float64[]
    alphas = Float64[]
    normF = Float64[]

    alpha1 = 0.0

    # perform Newton's method with linesearch based on Nocedal and Wright
    # chapter 3 section 5. 
    for n = 1:iterations

        if alpha1 == 1.0
            # if alpha was 1, we don't need to update the function 
            # because we have already calculated that in the last
            # loop. just update the jacobian. since we set alpha1=0
            # before the loop, this will never be called on the first
            # iteration. 
            calcfj!(nothing,Jsparse,x,wmodesm,wmodes2m,Rbnm,Rbnmt,invLnm,Cnm,Gnm,bnm,
            Ljb,Ljbm,Nmodes,
            Nbranches,Lmean,AoLjbmvector,AoLjbm,
            AoLjnmindexmap,invLnmindexmap,Gnmindexmap,Cnmindexmap)
        else
            # update the residual function and the Jacobian
            calcfj!(F,Jsparse,x,wmodesm,wmodes2m,Rbnm,Rbnmt,invLnm,Cnm,Gnm,bnm,
            Ljb,Ljbm,Nmodes,
            Nbranches,Lmean,AoLjbmvector,AoLjbm,
            AoLjnmindexmap,invLnmindexmap,Gnmindexmap,Cnmindexmap)
        end

        push!(normF,norm(F))

        # Jsparse.nzval .*= invnormJ
        # F .*= invnormJ

        # solve the linear system
        try 
            # update the factorization. the sparsity structure does not change
            # so we can reuse the factorization object.
            KLU.klu!(factorization,Jsparse)

            # solve the linear system            
            ldiv!(deltax,factorization,F)
        catch e
            if isa(e, SingularException)
                # reusing the symbolic factorization can sometimes lead to
                # numerical problems. if the first linear solve fails
                # try factoring and solving again
                factorization = KLU.klu(Jsparse)
                ldiv!(deltax,factorization,F)
            else
                throw(e)
            end
        end

        # multiply deltax by -1
        rmul!(deltax,-1)

        # calculate the objective function and the derivative of the objective
        # with respect to the scalar variable alpha which parameterizes the
        # path between the old x and the new x. 
        # Note: the dot product takes the complex conjugate of the first vector
        f = real(0.5*dot(F,F))
        dfdalpha = real(dot(F,Jsparse*deltax))

        # # evaluate the objective function at Nsample points
        # for alpha in range(0,1,Nsamples)
        #     calcfj!(F,nothing,x - alpha*x1,wmodesm,wmodes2m,Rbnm,Rbnmt,invLnm,Cnm,Gnm,bnm,
        #     Ljb,Ljbm,Nmodes,
        #     Nbranches,Lmean,AoLjbmvector,AoLjbm,
        #     AoLjnmindexmap,invLnmindexmap,Gnmindexmap,Cnmindexmap)    
        #     push!(samples,real(0.5*dot(F,F)))
        # end

        # evaluate the function at the trial point
        calcfj!(F,nothing,x+deltax,wmodesm,wmodes2m,Rbnm,Rbnmt,invLnm,Cnm,Gnm,bnm,
        Ljb,Ljbm,Nmodes,
        Nbranches,Lmean,AoLjbmvector,AoLjbm,
        AoLjnmindexmap,invLnmindexmap,Gnmindexmap,Cnmindexmap)

        # F .*= invnormJ

        fp = real(0.5*dot(F,F))

        # coefficients of the quadratic equation a*alpha^2+b*alpha+c to interpolate
        # f vs alpha
        a = -dfdalpha + fp - f
        b = dfdalpha
        c = f
        alpha1 = -b/(2*a)
        f1fit = -b*b/(4*a) + c


        if f1fit > fp
            f1fit = fp
            alpha1 = 1.0
        end
        # if the fitted alpha overshoots the size of the interval (from 0 to 1),
        # then set alpha to 1 and make a full length step. 
        if alpha1 > 1 || alpha1 <= 0
            alpha1 = 1.0
            f1fit = fp
        end

        # if a is zero, alpha1 will be NaN
        if abs2(a) == 0 
            alpha1 = 1.0
        end

        # switch to newton once the norm is small enough
        switchofflinesearchtol = 1e-5
        # switchofflinesearchtol = 0
        normx = norm(x)
        if sqrt(fp)/normx <= switchofflinesearchtol && sqrt(f)/normx <= switchofflinesearchtol
            alpha1 = 1.0
            # println("norm(F)/norm(phi): ",sqrt(fp)/norm(x))
        end
        # alpha1 = 1

        # update x
        x .+= deltax*alpha1

        if norm(F,Inf) <= ftol || ( norm(x) > 0 && norm(F)/norm(x) < ftol)
            # terminate iterations if infinity norm or relative norm are less
            # than ftol. check that norm(x) is greater than zero to avoid
            # divide by zero errors. 
            # println("converged to: infinity norm of : ",norm(F,Inf)," after ",n," iterations")
            # println("norm(F)/norm(phi): ",norm(F)/norm(x))
            break
        end

        if n == iterations
            println("Warning: Solver did not converge with infinity norm of : ",norm(F,Inf)," after maximum iterations of ", n)
            println("norm(F)/norm(x): ", norm(F)/norm(x))
        end
    end


    phin = x
    out = nothing


    # calculate the scattering parameters for the pump
    # Nports = length(cdict[:P])
    Nports = length(portindices)
    # input = zeros(Complex{Float64},Nports*Nmodes)
    # input = Diagonal(zeros(Complex{Float64},Nports*Nmodes))

    input = zeros(Complex{Float64},Nports*Nmodes)    
    output = zeros(Complex{Float64},Nports*Nmodes)
    phibports = zeros(Complex{Float64},Nports*Nmodes)
    S = zeros(Complex{Float64},Nports*Nmodes,Nports*Nmodes)

    # calculate the scattering parameters
    if any(Ip .> 0)
        # calcS!(S,input/Lmean,output)
    end

    return NonlinearHB(out,phin,Rbnm,Ljb,Lb,Ljbm,Nmodes,Nbranches,S)
    # return (samples=samples,fmin=fmin,alphas=alphas,
    #     fvals=fvals,fpvals=fpvals,dfdalphavals=dfdalphavals)
    # return normF

end

struct NonlinearHB
    out
    phin
    Rbnm
    Ljb
    Lb
    Ljbm
    Nmodes
    Nbranches
    S
end


"""
    calcfj(F,J,phin,wmodesm,wmodes2m,Rbnm,invLnm,Cnm,Gnm,bm,Ljb,Ljbindices,
        Ljbindicesm,Nmodes,Lmean,AoLjbm)
        
Calculate the residual and the Jacobian. These are calculated with one function
in order to reuse the time domain nonlinearity calculation.

Leave off the type signatures on F and J because the solver will pass a type of
Nothing if it only wants to calculate F or J. 

"""
function calcfj!(F,
        J,
        phin::AbstractVector,
        wmodesm::AbstractMatrix, 
        wmodes2m::AbstractMatrix, 
        Rbnm::AbstractArray{Int64,2}, 
        Rbnmt::AbstractArray{Int64,2}, 
        invLnm::AbstractMatrix, 
        Cnm::AbstractMatrix, 
        Gnm::AbstractMatrix, 
        bnm::AbstractVector, 
        Ljb::SparseVector, 
        Ljbm::SparseVector,
        Nmodes::Int64, 
        Nbranches::Int64,
        Lmean,
        AoLjbmvector::AbstractVector,
        AoLjbm,AoLjnmindexmap,invLnmindexmap,Gnmindexmap,Cnmindexmap)

    # convert from a node flux to a branch flux
    phib = Rbnm*phin

    # calculate the sine and cosine nonlinearities 
    Am = sincosnloddtoboth(phib[Ljbm.nzind],nnz(Ljb),Nmodes)

    # calculate the residual
    if !(F == nothing)

        # calculate the function. use the sine terms. Am[2:2:2*Nmodes,:]
        # calculate  AoLjbm, this is just a diagonal matrix. 
        # check if this is consistent with other calculations
        @inbounds for i = 1:nnz(Ljb)
            for j = 1:Nmodes
                AoLjbmvector[(Ljb.nzind[i]-1)*Nmodes + j] = Am[2*j,i]*(Lmean/Ljb.nzval[i])

            end
        end

        F .= Rbnmt*AoLjbmvector + (invLnm + im*Gnm*wmodesm - Cnm*wmodes2m)*phin - bnm
    end

    # calculate the Jacobian
    if !(J == nothing)

        # calculate  AoLjbm
        updateAoLjbm!(AoLjbm,Am,Ljb,Lmean,Nmodes,Nbranches)

        # convert to a sparse node matrix
        AoLjnm = Rbnmt*AoLjbm*Rbnm

        # calculate the sparse Jacobian.
        # @time J .= AoLjnm .+ invLnm .- im.*Gnm*wmodesm .- Cnm*wmodes2m
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
    calcAoLjbm(Am,Ljb::SparseVector,Lmean,Nmodes,Nbranches)


# Examples
```jldoctest
julia> @variables Lj1 Lj2 A11 A12 A21 A22 A31 A32;JosephsonCircuits.calcAoLjbm([A11;A21;A31],JosephsonCircuits.SparseArrays.sparsevec([1],[Lj1]),1,2,1)
2×2 SparseArrays.SparseMatrixCSC{Num, Int64} with 4 stored entries:
 A11 / Lj1  A31 / Lj1
 A31 / Lj1  A11 / Lj1

julia> @variables Lj1 Lj2 A11 A12 A21 A22 A31 A32;JosephsonCircuits.calcAoLjbm([A11 A12;A21 A22;A31 A32],JosephsonCircuits.SparseArrays.sparsevec([1,2],[Lj1,Lj2]),1,2,2)
4×4 SparseArrays.SparseMatrixCSC{Num, Int64} with 8 stored entries:
 A11 / Lj1  A31 / Lj1          ⋅          ⋅
 A31 / Lj1  A11 / Lj1          ⋅          ⋅
         ⋅          ⋅  A12 / Lj2  A32 / Lj2
         ⋅          ⋅  A32 / Lj2  A12 / Lj2
```
"""
function calcAoLjbm(Am,Ljb::SparseVector,Lmean,Nmodes,Nbranches)

    # define empty vectors for the rows, columns, and values
    I = Vector{eltype(Ljb.nzind)}(undef,nnz(Ljb)*Nmodes^2)
    J = Vector{eltype(Ljb.nzind)}(undef,nnz(Ljb)*Nmodes^2)

    type = promote_type(eltype(Am),eltype(1 ./Ljb.nzval))

    if type <: Symbolic
        type = Any
    end

    V = Vector{type}(undef,nnz(Ljb)*Nmodes^2)


    if size(Am,2) != nnz(Ljb)
        throw(DimensionError("The second axis of Am must equal the number of nonzero elements in Ljb (the number of JJs)."))
    end

    if length(Ljb) > Nbranches
        throw(DimensionError("The length of Ljb cannot be larger than the number of branches."))
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
    updateAoLjbm!(AoLjbm::SparseMatrixCSC,Am,Ljb::SparseVector,Lmean,Nmodes,Nbranches)

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
"""
function updateAoLjbm!(AoLjbm::SparseMatrixCSC,Am,Ljb::SparseVector,Lmean,Nmodes,Nbranches)

    # check that there are the right number of nonzero values. 
    # check that the dimensions are consistent with Nmode and Nbranches.

    if nnz(Ljb)*Nmodes^2 != nnz(AoLjbm)
        throw(DimensionError("The number of nonzero elements in AoLjbm are not consistent with nnz(Ljb) and Nmodes."))
    end

    if size(Am,2) != nnz(Ljb)
        throw(DimensionError("The second axis of Am must equal the number of nonzero elements in Ljb (the number of JJs)."))
    end

    if length(Ljb) > Nbranches
        throw(DimensionError("The length of Ljb cannot be larger than the number of branches."))
    end

    # i want a vector length(Ljb) where the indices are the values Ljb.nzind
    # and the values are the indices of Ljb.nzind
    indexconvert = zeros(Int64,length(Ljb))
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
    sincosnloddtoboth(amodd::Array{Complex{Float64},1},Nbranches::Int64,m::Int64)

Applies the junction nonlinearity to a vector of branch fluxes of length Nbranches*m
where m is the number of odd pump harmonics (1w, 3w, 5w, etc). The ordering is 
(mode 1, node 1), (mode 2, node 1) ... (mode 1, node 2) ... Returns even AND
odd terms in a 2d array with dimensions 2*m by Nbranches. 

# Examples
```
julia> JosephsonCircuits.sincosnloddtoboth([0.5+0.0im,0,0,0],1,4)
8×1 Matrix{ComplexF64}:
      0.765197686557965 + 0.0im
    0.44005058574495276 + 0.0im
   -0.11490348493140044 + 0.0im
    -0.0195633539946485 + 0.0im
  0.0024766387010484933 - 0.0im
  0.0002497629794613717 + 0.0im
  -2.084411456060309e-5 + 0.0im
 -3.0046516347986037e-6 + 0.0im

julia> JosephsonCircuits.sincosnloddtoboth([0.02+0.0im,0,0.01+0.0im,0],2,2)
4×2 Matrix{ComplexF64}:
       0.9996+0.0im       0.9999+0.0im
     0.019996-0.0im    0.0099995-0.0im
 -0.000199967+0.0im  -4.99979e-5+0.0im
   -2.6664e-6+0.0im  -3.33325e-7+0.0im
```
"""
function sincosnloddtoboth(amodd::Array{Complex{Float64},1},
    Nbranches::Int64,m::Int64)

    if length(amodd) != Nbranches*m
        error("Length of node flux vector not consistent with number of modes
        number of frequencies")
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
```
julia> JosephsonCircuits.sincosnl([0 0;0.001+0im 0;0 0;0 0; 0 0])
5×2 Matrix{ComplexF64}:
     0.999999+0.0im  1.0+0.0im
        0.001-0.0im  0.0-0.0im
      -5.0e-7+0.0im  0.0+0.0im
 -1.66667e-10+0.0im  0.0+0.0im
  8.33777e-14+0.0im  0.0+0.0im

julia> JosephsonCircuits.sincosnl([0 0;0.001+0im 0.25+0im;0 0;0 0; 0 0])
5×2 Matrix{ComplexF64}:
     0.999999+0.0im      0.93847+0.0im
        0.001-0.0im     0.242268-0.0im
      -5.0e-7+0.0im   -0.0306044+0.0im
 -1.66667e-10+0.0im  -0.00255568+0.0im
  8.33777e-14+0.0im  0.000321473+0.0im
```
"""
function sincosnl(am::Array{Complex{Float64},2})

    #choose the number of time points based on the number of fourier
    #coefficients
    #stepsperperiod = 2*size(am)[1]-2
    # changed to below because i wasn't using enough points when Nmodes=1.
    # the results contained only real values. 
    if size(am)[1] == 2
        stepsperperiod = 2*size(am)[1]-1
    else
        stepsperperiod = 2*size(am)[1]-2
    end        

    #transform back to time domain
    ift = FFTW.irfft(am,stepsperperiod,[1])*stepsperperiod

    #apply the nonlinear function
    nlift = cos.(ift) .+ sin.(ift)
    
    #fourier transform back to the frequency domain
    ftnlift = FFTW.rfft(nlift,[1])/stepsperperiod
    
    return ftnlift

end


function applynl(am::Array{Complex{Float64}},f::Function)

    #choose the number of time points based on the number of fourier
    #coefficients
    # changed to below because i wasn't using enough points when Nmodes=1.
    # the results contained only real values. 
    if size(am,1) == 2
        stepsperperiod = 2*size(am,1)-1
    else
        stepsperperiod = 2*size(am,1)-2
    end        

    #transform back to time domain
    ift = FFTW.irfft(am,stepsperperiod,1:length(size(am))-1)

    # normalize the fft
    normalization = prod(size(ift)[1:end-1])
    ift .*= normalization

    #apply the nonlinear function
    nlift = f.(ift)
    
    #fourier transform back to the frequency domain
    ftnlift = FFTW.rfft(nlift,1:length(size(am))-1)/normalization
    
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
    calcw(ws::Number,i::Integer,wp::Number)

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
function calcw(ws,i::Integer,wp)
    return ws + 2*i*wp
end
function calcw(ws,i::AbstractVector,wp)

    w = zeros(typeof(ws),length(i))
    calcw!(ws,i,wp,w)
    return w
end
function calcw(ws::AbstractVector,i::AbstractVector,wp)

    w = zeros(eltype(ws),length(ws),length(i))
    calcw!(ws,i,wp,w)
    return w
end

"""
    calcw!(ws,i,wp,w)

Generate the signal and idler frequencies using the formula ws + 2*i*wp.
Overwrites w with output. 
"""
function calcw!(ws,i::AbstractVector,wp,w::AbstractVector)

    for j in 1:length(i)
        w[j] = calcw(ws,i[j],wp)
    end

    return nothing
end
function calcw!(ws::AbstractVector,i::AbstractVector,wp,w::AbstractMatrix)

    for j in 1:length(ws)
        for k in 1:length(i)
            w[j,k] = calcw(ws[j],i[k],wp)
        end
    end

    return nothing
end
