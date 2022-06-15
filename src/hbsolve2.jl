

"""

    hbsolve2()

New version of the harmonic balance solver suitable for arbitrary numbers of 
ports, sources, and pumps. Still under development. 
"""

function hbsolve2(ws,wp,Ip,Nsignalmodes,Npumpmodes,circuit,circuitdefs; pumpports=[1],
    solver =:klu, iterations=100,symfreqvar = nothing, sorting=:number,
    nbatches = Base.Threads.nthreads())
    # solve the nonlinear system    
    # use the old syntax externally

    w = (wp,)
    Nharmonics = (2*Npumpmodes,)
    sources = ((w=w[1],port=pumpports[1],current=Ip),)
    pump=hbnlsolve2(w,Nharmonics,sources,circuit,circuitdefs,
        solver=solver,odd=true,dc=false,even=false,symfreqvar = symfreqvar,
        sorting=sorting)
    # pump=hbnlsolve2(w,Nharmonics,sources,circuit,circuitdefs,
    #     solver=solver,odd=true,dc=false,even=false)
    # the node flux
    phin = pump.out.zero

    # convert from node flux to branch flux
    phib = pump.Rbnm*phin

    # calculate the sine and cosine nonlinearities from the pump flux
    # Am = sincosnloddtoboth(phib[pump.Ljbm.nzind],length(pump.Ljb.nzind),pump.Nmodes)
    Am = sincosnloddtoboth(phib[pump.Ljbm.nzind],length(pump.Ljb.nzind),Npumpmodes)

    # solve the linear system
    signal=hblinsolve(ws,wp,Nsignalmodes,circuit,circuitdefs,Am=Am,solver=solver,
        symfreqvar = symfreqvar, nbatches = nbatches, sorting = sorting)
    # signal=hblinsolve(ws,wp,Nsignalmodes,circuit,circuitdefs,Am=Am,solver=solver)

    # return (pump=pump, Am=Am, signal=signal)
    return HB(pump,Am,signal)
end

"""
    hbnlsolve(wp,Ip,Nmodes,circuit,circuitdefs)
"""

function hbnlsolve2(w::Tuple,Nharmonics::Tuple,sources::NamedTuple,circuit,circuitdefs;ports=[1],
    solver=:klu,iterations=100,maxintermodorder=Inf,dc=false,odd=true,even=false,
    x0=nothing,symfreqvar = nothing,  sorting=:number)

    return hbnlsolve2(w,Nharmonics,(sources,),circuit,circuitdefs;
    solver=solver,iterations=iterations,maxintermodorder=maxintermodorder,
    dc=dc,odd=odd,even=even,x0=x0,freqvar = freqvar, sorting=sorting)
end

# function hbnlsolve2(w::Tuple,Nharmonics::Tuple,sources::NamedTuple,circuit,circuitdefs;ports=[1],
#     solver=:klu,iterations=100,maxintermodorder=Inf,dc=false,odd=true,even=false,
#     x0=nothing)

#     return hbnlsolve2(w,Nharmonics,(sources,),circuit,circuitdefs;
#     solver=solver,iterations=iterations,maxintermodorder=maxintermodorder,
#     dc=dc,odd=odd,even=even,x0=x0)
# end

function hbnlsolve2(w::Tuple,Nharmonics::Tuple,sources::Tuple,circuit,circuitdefs;
    solver=:klu,iterations=100,maxintermodorder=Inf,dc=false,odd=true,even=false,
    x0=nothing,symfreqvar = nothing, sorting= sorting)

# function hbnlsolve2(w::Tuple,Nharmonics::Tuple,sources::Tuple,circuit,circuitdefs;
#     solver=:klu,iterations=100,maxintermodorder=Inf,dc=false,odd=true,even=false,
#     x0=nothing)

    # calculate the 
    Nw,coords,values,dropcoords,dropvalues = calcfrequencies(w,Nharmonics,
        maxintermodorder=maxintermodorder,dc=dc,even=even,odd=odd)
    Nt=NTuple{length(Nw),Int64}(ifelse(i == 1, 2*val-1, val) for (i,val) in enumerate(Nw))

    dropdict = Dict(dropcoords .=> dropvalues)

    values2 = calcfrequencies2(Nt,coords,values);

    freqindexmap,conjsourceindices,conjtargetindices = calcphiindices(Nt,dropdict)

    indices = calcrdftsymmetries(Nt)

    # assign the frequencies
    wmodes = values2[:]
    Nmodes = length(wmodes)

    # parse the circuit components
    c0 = parsecircuit(circuit)

    numbervector = valuevectortonumber(c0.valuevector,circuitdefs)

    # sort the nodes
    uniquenodevectorsorted,nodeindexarraysorted = sortnodes(c0.uniquenodevector,
        c0.nodeindexvector,sorting=:number)

    branchvector = extractbranches(c0.typevector,nodeindexarraysorted)

    # calculate the graph of inductive components glelist, the
    # superconducting spanning tree selist, and the list of loop
    # indices celist. 
    g = calcgraphs(branchvector,length(uniquenodevectorsorted))

    portdict=componentdictionaryP(c0.typevector,nodeindexarraysorted,
        c0.mutualinductorvector,numbervector)

    resistordict=componentdictionaryR(c0.typevector,nodeindexarraysorted,
        c0.mutualinductorvector,numbervector)


    # # calculate the capacitance, inductance, and inverse inductances matrices
    # m = calcmatrices(c0.typevector,nodeindexarraysorted,numbervector,c0.namedict,
    #     c0.mutualinductorvector,g.edge2indexdict,g.Rbn,Nmodes,
    #     length(uniquenodevectorsorted),g.Nbranches)

    # Nnodes = length(uniquenodevectorsorted)
    # Nbranches = g.Nbranches
    # Lmean = m.Lmean
    # invLnm = m.invLnm*Lmean
    # Gnm = m.Gnm*Lmean
    # Cnm = m.Cnm*Lmean
    # Ljb = m.Ljb
    # Ljbm = m.Ljbm
    # Lb = m.Lb
    # Rbnm = m.Rbnm


    typevector = c0.typevector
    namedict = c0.namedict
    mutualinductorvector = c0.mutualinductorvector
    edge2indexdict = g.edge2indexdict
    Rbn = g.Rbn
    Nnodes = length(c0.uniquenodevector)
    Nbranches = g.Nbranches

    # calculate Lmean
    Lmean = calcLmean(typevector,numbervector)
    # Lmean = 1.0

    # capacitance matrix
    Cnm = calcCn(typevector,nodeindexarraysorted,numbervector,Nmodes,Nnodes)

    # conductance matrix
    Gnm = calcGn(typevector,nodeindexarraysorted,numbervector,Nmodes,Nnodes)

    # branch inductance vectors
    Lb = calcLb(typevector,nodeindexarraysorted,numbervector,edge2indexdict,1,Nbranches)
    Lbm = calcLb(typevector,nodeindexarraysorted,numbervector,edge2indexdict,Nmodes,Nbranches)

    # branch Josephson inductance vectors
    Ljb = calcLjb(typevector,nodeindexarraysorted,numbervector,edge2indexdict,1,Nbranches)
    Ljbm = calcLjb(typevector,nodeindexarraysorted,numbervector,edge2indexdict,Nmodes,Nbranches)

    # mutual branch inductance matrix
    Mb = calcMb(typevector,nodeindexarraysorted,numbervector,namedict,
        mutualinductorvector,edge2indexdict,1,Nbranches)

    # inverse nodal inductance matrix from branch inductance vector and branch
    # inductance matrix
    invLnm = calcinvLn(Lb,Mb,Rbn,Nmodes)

    # expand the size of the indicidence matrix
    Rbnm = diagrepeat(Rbn,Nmodes)


    # calculate the diagonal frequency matrices
    wmodesm = Diagonal(repeat(wmodes,outer=Nnodes-1))
    wmodes2m = Diagonal(repeat(wmodes.^2,outer=Nnodes-1))

    # calculate the source terms in the branch basis
    bbm = zeros(Complex{Float64},Nbranches*Nmodes)    
    # for (key,val) in sort(collect(cdict[:P]), by=x->x[2])

    for source in sources
        for (key,val) in portdict
            if val == source[:port]
                # now i have to find the index.
                # this will depend on the frequency index
                # i should calculate that in the right way now. 
                for i = 1:length(values)
                    if values[i] == source[:w]
                        bbm[(g.edge2indexdict[key]-1)*Nmodes+i] = Lmean*source[:current]/phi0
                        break
                    end
                end
                break
            end
        end

    end


    # convert from the node basis to the branch basis
    bnm = transpose(Rbnm)*bbm

    Nwtuple=NTuple{length(Nt)+1,Int64}(
            if i==1 
                (Nt[i] ÷ 2) +1 
            elseif i <= length(Nt)
                Nt[i] 
            else 
                length(Ljb.nzval)
            end 
            for i = 1:length(Nt) + 1
        )
    phimatrix = zeros(Complex{Float64},Nwtuple)


    Amatrix = rand(Complex{Float64},Nwtuple)

    # calculate AoLjbm, the sparse branch AoLj matrix
    AoLjbm = calcAoLjbm2(Amatrix,Ljb,Lmean,Nmodes,Nbranches,freqindexmap)
    AoLjbmcopy = calcAoLjbm2(Amatrix,Ljb,Lmean,Nmodes,Nbranches,freqindexmap)

    # convert to a sparse node matrix. Note: I was having problems with type 
    # instability when i used AoLjbm here instead of AoLjbmcopy. 
    AoLjnm = transpose(Rbnm)*AoLjbmcopy*Rbnm;

    if x0 == nothing
        x = zeros(Complex{Float64},(Nnodes-1)*Nmodes)
    else
        x = x0
    end
    F = zeros(Complex{Float64},(Nnodes-1)*Nmodes)
    AoLjbm2 = zeros(Complex{Float64},Nbranches*Nmodes)

    # make a sparse transpose (improves multiplication speed slightly)
    Rbnmt = sparse(transpose(Rbnm))

    # substitute in the mode frequencies for components which have frequency
    # defined symbolically.
    Cnm = freqsubst(Cnm,wmodes,symfreqvar)
    Gnm = freqsubst(Gnm,wmodes,symfreqvar)
    invLnm = freqsubst(invLnm,wmodes,symfreqvar)

    # scale the matrices for numerical reasons
    Cnm *= Lmean
    Gnm *= Lmean
    invLnm *= Lmean

    # Calculate something with the same sparsity structure as the Jacobian.
    # Don't bother multiplying by the diagonal frequency matrices since they
    # won't change the sparsity structure. 
    # Jsparse = (AoLjnm + invLnm - im.*Gnm*wmodesm - Cnm*wmodes2m)
    Jsparse = (AoLjnm + invLnm - im*Gnm - Cnm)

    # make the index maps so we can add the sparse matrices together without
    # memory allocations. 
    AoLjnmindexmap = sparseaddmap(Jsparse,AoLjnm)
    invLnmindexmap = sparseaddmap(Jsparse,invLnm)
    Gnmindexmap = sparseaddmap(Jsparse,Gnm)
    Cnmindexmap = sparseaddmap(Jsparse,Cnm)


    function FJsparse2!(F,J,x)
        calcfj2!(F,J,x,wmodesm,wmodes2m,Rbnm,Rbnmt,invLnm,Cnm,Gnm,bnm,
        Ljb,Ljbm,Nmodes,
        Nbranches,Lmean,AoLjbm2,AoLjbm,
        AoLjnmindexmap,invLnmindexmap,Gnmindexmap,Cnmindexmap,
        freqindexmap,conjsourceindices,conjtargetindices,phimatrix)
        return F,J
    end

    if solver == :klu
        # perform a factorization. this will be updated later for each 
        # interation
        FK = KLU.klu(Jsparse)

        odsparse = NLsolve.OnceDifferentiable(NLsolve.only_fj!(FJsparse2!),x,F,Jsparse)

        # if the sparsity structure doesn't change, we can cache the 
        # factorization. this is a significant speed improvement.
        # out=NLsolve.nlsolve(odsparse,method = :trust_region,autoscale=false,x,iterations=iterations,linsolve=(x, A, b) ->(KLU.klu!(FK,A);ldiv!(x,FK,b)) )
        out=NLsolve.nlsolve(odsparse,method = :trust_region,autoscale=false,x,iterations=iterations,linsolve=(x, A, b) ->(KLU.klu!(FK,A);ldiv!(x,FK,b)) )

    else
        error("Error: Unknown solver")
    end

    if out.f_converged == false
        println("Nonlinear solver not converged. You may need to supply a better
        guess at the solution vector, increase the number of pump harmonics, or
        increase the number of iterations.")
    end

    # calculate the scattering parameters for the pump
    # Nports = length(cdict[:P])
    Nports = length(portdict)
    # input = zeros(Complex{Float64},Nports*Nmodes)
    input = Diagonal(zeros(Complex{Float64},Nports*Nmodes))

    output = zeros(Complex{Float64},Nports*Nmodes)
    phibports = zeros(Complex{Float64},Nports*Nmodes)
    inputval = zero(Complex{Float64})
    S = zeros(Complex{Float64},Nports*Nmodes,Nports*Nmodes)

    phin = out.zero

    # calculate the branch fluxes
    # calcphibports!(phibports,phin,portdict,Nmodes)

    # i don't think i need to do the scattering parameter analysis for the two
    # tone harmonic balance.

    # # calculate the input and output voltage waves at each port
    # calcinput!(input,Ip,phibports,portdict,resistordict,wmodes)
    # calcoutput!(output,phibports,portdict,resistordict,wmodes)

    # calculate the scattering parameters
    # calcS!(S,input,output,phibports,portdict,resistordict,wmodes)
    # calcS!(S,input,output)

    # return (
    #     out = out,
    #     wmodes = wmodes, 
    #     wmodesm = wmodesm,
    #     wmodes2m=wmodes2m,
    #     bbm = bbm,
    #     bnm = bnm,
    #     portdict = portdict,
    #     values = values,
    #     # Avector = Avector,
    #     # Amatrix = Amatrix,
    #     freqindexmap = freqindexmap,
    #     conjsourceindices = conjsourceindices,
    #     conjtargetindices = conjtargetindices,
    #     Nwtuple=Nwtuple,
    #     Nt = Nt,
    #     Nw = Nw,
    #     dropcoords = dropcoords,
    #     dropvalues = dropvalues,
    #     AoLjnm = AoLjnm,
    #     Rbnm = Rbnm,
    #     Ljbm = Ljbm,
    #     Ljb = Ljb,
    #     )
    return NonlinearHB(out,phin,Rbnm,Ljb,Lb,Ljbm,Nmodes,Nbranches,S)

    # return (out=out,
    #     wmodesm = wmodesm,
    #     wmodes2m = wmodes2m,
    #     Rbnm = Rbnm,
    #     invLnm = invLnm,
    #     Cnm = Cnm,
    #     Gnm = Gnm,
    #     bnm = bnm,
    #     Ljb = Ljb,
    #     Lb = Lb,
    #     Ljbm = Ljbm,
    #     Nmodes = Nmodes,
    #     Nbranches = Nbranches,
    #     Lmean = Lmean,
    #     AoLjbm = AoLjbm,
    #     AoLjnm = AoLjnm,
    #     II = II,
    #     JJ = JJ,
    #     KK = KK,
    #     # J = J,
    #     S=S,
    #     )

end


"""
    calcfj(F,J,phin,wmodesm,wmodes2m,Rbnm,invLnm,Cnm,Gnm,bm,Ljb,Ljbindices,
        Ljbindicesm,Nmodes,Lmean,AoLjbm)
        
Calculate the residual and the Jacobian. These are calculated with one function
in order to reuse the time domain nonlinearity calculation.

Leave off the type signatures on F and J because the solver will pass a type of
Nothing if it only wants to calculate F or J. 

"""
function calcfj2!(F,
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
        AoLjbm,AoLjnmindexmap,invLnmindexmap,Gnmindexmap,Cnmindexmap,
        freqindexmap,conjsourceindices,conjtargetindices,phimatrix)

    # convert from a node flux to a branch flux
    phib = Rbnm*phin

    # phib[Ljbm.nzind] are the branch fluxes for each of the JJ's
    phivectortomatrix!(phib[Ljbm.nzind],phimatrix,freqindexmap,conjsourceindices,
        conjtargetindices,length(Ljb.nzval))

    if !(F == nothing)

        Amsin = applynl(phimatrix,(x) -> sin(x))
        Nfreq = prod(size(Amsin)[1:end-1])

        # calculate the function. use the sine terms. Am[2:2:2*Nmodes,:]
        # calculate  AoLjbm, this is just a diagonal matrix. 
        for i = 1:nnz(Ljb)
            for j = 1:Nmodes
                AoLjbmvector[(Ljb.nzind[i]-1)*Nmodes + j] = Amsin[freqindexmap[j]+Nfreq*(i-1)]*(Lmean/Ljb.nzval[i])
            end
        end

        F .= Rbnmt*AoLjbmvector + (invLnm + im*Gnm*wmodesm - Cnm*wmodes2m)*phin - bnm
    end

    #calculate the Jacobian
    if !(J == nothing)

        # calculate  AoLjbm
        Amcos = applynl(phimatrix,(x) -> cos(x))
        # AoLjbm = calcAoLjbm2(Amcos,Ljb,Ljbindices,Lmean,Nmodes,Nbranches,II,JJ,KK,freqindexmap)
        updateAoLjbm2!(AoLjbm,Amcos,Ljb,Lmean,Nmodes,Nbranches,freqindexmap)

        # convert to a sparse node matrix
        AoLjnm = Rbnmt*AoLjbm*Rbnm

        # calculate the Jacobian. If J is sparse, keep it sparse. 
        # J .= AoLjnm + invLnm - im*Gnm*wmodesm - Cnm*wmodes2m
        # the code below adds the sparse matrices together with minimal
        # memory allocations and without changing the sparsity structure.
        fill!(J,zero(eltype(J)))
        sparseadd!(J,AoLjnm,AoLjnmindexmap)
        sparseadd!(J,invLnm,invLnmindexmap)
        sparseadd!(J,im,Gnm,wmodesm,Gnmindexmap)
        sparseadd!(J,-1,Cnm,wmodes2m,Cnmindexmap)
    end
    return nothing
end



function calcAoLjbm2(Am,Ljb::SparseVector,Lmean,Nmodes,Nbranches,freqindexmap)

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

    Nfreq = prod(size(Am)[1:end-1])

    # calculate  AoLjbm
    n = 1
    for i = 1:nnz(Ljb)
        for j = 1:Nmodes
            for k = 1:Nmodes

                # calculate the toeplitz matrices for each node 
                I[n]=j+(Ljb.nzind[i]-1)*Nmodes
                J[n]=k+(Ljb.nzind[i]-1)*Nmodes

                # assume terms we don't have pump data for are zero.
                # index = 2*abs(j-k)+1
                index = abs(freqindexmap[j]-freqindexmap[k])+1
                if index > size(Am,1)
                    V[n] = 0
                else
                    # V[n]=Am[index,i]*(Lmean/Ljb.nzval[i])
                    V[n]=Am[index+Nfreq*(i-1)]*(Lmean/Ljb.nzval[i])

                end

                #take the complex conjugate of the upper half (not the diagonal)
                if j-k<0
                    V[n] = conj(V[n])
                end

                ## for debugging. calculaet index, Ljb.nzind and i from I and J.
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

Update the values in the sparse AoLjbm matrix in place.

"""
function updateAoLjbm2!(AoLjbm::SparseMatrixCSC,Am,Ljb::SparseVector,Lmean,Nmodes,Nbranches,freqindexmap)

    Nfreq = prod(size(Am)[1:end-1])


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
            j = AoLjbm.rowval[m]-(Ljb.nzind[i]-1)*Nmodes
            k = l-(Ljb.nzind[i]-1)*Nmodes

            index = abs(freqindexmap[j]-freqindexmap[k])+1

            if index > size(Am,1)
                AoLjbm.nzval[m] = 0
            else
                # AoLjbm.nzval[m]=Am[index,i]*(Lmean/Ljb.nzval[i])
                AoLjbm.nzval[m]=Am[index+Nfreq*(i-1)]*(Lmean/Ljb.nzval[i])

            end

            #take the complex conjugate of the upper half (not the diagonal)
            if AoLjbm.rowval[m]-l<0
                AoLjbm.nzval[m] = conj(AoLjbm.nzval[m])
            end

        end
    end

    return nothing
end
