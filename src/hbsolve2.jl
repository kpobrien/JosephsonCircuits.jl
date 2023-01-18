

"""

    hbsolve2()

New version of the harmonic balance solver suitable for arbitrary numbers of 
ports, sources, and pumps. Still under development. 
"""

# function hbsolve2(ws,wp,Ip,Nsignalmodes,Npumpmodes,circuit,circuitdefs; pumpports=[1],
#     solver =:klu, iterations=1000,ftol=1e-8,symfreqvar = nothing, sorting=:number,
#     nbatches = Base.Threads.nthreads())
function hbsolve2(ws, wp, Ip, Nsignalmodes, Npumpmodes, circuit, circuitdefs;
    pumpports = [1], solver = :klu, iterations = 1000, ftol = 1e-8,
    symfreqvar = nothing, nbatches = Base.Threads.nthreads(), sorting = :number,
    returnS = true, returnSnoise = false, returnQE = true, returnCM = true,
    returnnodeflux = false, returnvoltage = false, verbosity = 1)
    # solve the nonlinear system
    # use the old syntax externally

    # parse and sort the circuit
    psc = parsesortcircuit(circuit, sorting=sorting)

    # calculate the circuit graph
    cg = calccircuitgraph(psc)

    w = (wp,)
    Nharmonics = (2*Npumpmodes,)
    sources = ((w=w[1],port=pumpports[1],current=Ip),)
    pump=hbnlsolve2(w,Nharmonics,sources,circuit,circuitdefs,
        solver=solver,odd=true,dc=false,even=false,symfreqvar = symfreqvar,
        sorting=sorting,ftol=ftol)

    # the node flux
    # phin = pump.out.zero
    nodeflux = pump.nodeflux    

    # convert from node flux to branch flux
    phib = pump.Rbnm*nodeflux

    # calculate the sine and cosine nonlinearities from the pump flux
    Am = sincosnloddtoboth(phib[pump.Ljbm.nzind], length(pump.Ljb.nzind), pump.Nmodes)

    # solve the linear system
    signal=hblinsolve(ws, psc, cg, circuitdefs, wp = wp, Nmodes = Nsignalmodes,
        Am = Am, solver = solver, symfreqvar = symfreqvar, nbatches = nbatches,
        returnS = returnS, returnSnoise = returnSnoise, returnQE = returnQE,
        returnCM = returnCM, returnnodeflux = returnnodeflux,
        returnvoltage = returnvoltage, verbosity = verbosity)
    return HB(pump, Am, signal)
end

"""
    hbnlsolve(wp,Ip,Nmodes,circuit,circuitdefs)
"""

function hbnlsolve2(w::Tuple,Nharmonics::Tuple,sources::NamedTuple,circuit,circuitdefs;ports=[1],
    solver=:klu,iterations=1000,maxintermodorder=Inf,dc=false,odd=true,even=false,
    x0=nothing,ftol=1e-8,symfreqvar = nothing,  sorting=:number)

    return hbnlsolve2(w,Nharmonics,(sources,),circuit,circuitdefs;
    solver=solver,iterations=iterations,maxintermodorder=maxintermodorder,
    dc=dc,odd=odd,even=even,x0=x0,ftol=1e-8,symfreqvar = symfreqvar, sorting=sorting)
end

# function hbnlsolve2(w::Tuple,Nharmonics::Tuple,sources::NamedTuple,circuit,circuitdefs;ports=[1],
#     solver=:klu,iterations=100,maxintermodorder=Inf,dc=false,odd=true,even=false,
#     x0=nothing)

#     return hbnlsolve2(w,Nharmonics,(sources,),circuit,circuitdefs;
#     solver=solver,iterations=iterations,maxintermodorder=maxintermodorder,
#     dc=dc,odd=odd,even=even,x0=x0)
# end

function hbnlsolve2(w::Tuple,Nharmonics::Tuple,sources::Tuple,circuit,circuitdefs;
    solver=:klu,iterations=1000,maxintermodorder=Inf,dc=false,odd=true,even=false,
    x0=nothing,ftol=1e-8,symfreqvar = nothing, sorting= :number)

# function hbnlsolve2(w::Tuple,Nharmonics::Tuple,sources::Tuple,circuit,circuitdefs;
#     solver=:klu,iterations=100,maxintermodorder=Inf,dc=false,odd=true,even=false,
#     x0=nothing)

    # calculate the 
    Nw,coords,values,dropcoords,dropvalues = calcfrequencies(w,Nharmonics,
        maxintermodorder=maxintermodorder,dc=dc,even=even,odd=odd)
    Nt=NTuple{length(Nw),Int}(ifelse(i == 1, 2*val-1, val) for (i,val) in enumerate(Nw))

    dropdict = Dict(dropcoords .=> dropvalues)

    values2 = calcfrequencies2(Nt,coords,values);

    freqindexmap,conjsourceindices,conjtargetindices = calcphiindices(Nt,dropdict)

    indices = calcrdftsymmetries(Nt)

    # assign the frequencies
    wmodes = values2[:]
    Nmodes = length(wmodes)



    # parse and sort the circuit
    psc = parsesortcircuit(circuit,sorting=sorting)

    # calculate the circuit graph
    cg = calccircuitgraph(psc)


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
    # portdict = nm.portdict

    portindices = nm.portindices
    portnumbers = nm.portnumbers
    portimpedanceindices = nm.portimpedanceindices

    # noiseportimpedanceindices = nm.noiseportimpedanceindices

    # resistordict = nm.resistordict
    Lmean = nm.Lmean
    # Lmean = 1.0
    Lb = nm.Lb

    # calculate the diagonal frequency matrices
    wmodesm = Diagonal(repeat(wmodes,outer=Nnodes-1))
    wmodes2m = Diagonal(repeat(wmodes.^2,outer=Nnodes-1))

    # calculate the source terms in the branch basis
    bbm = zeros(Complex{Float64},Nbranches*Nmodes)  
    for source in sources
        # for (key,val) in portdict
        for (i,portindex) in enumerate(portindices)
            portnumber = portnumbers[i]
            key = (nodeindexarraysorted[1,portindex],nodeindexarraysorted[2,portindex])

            if portnumber == source[:port]
                # now i have to find the index.
                # this will depend on the frequency index
                # i should calculate that in the right way now. 
                for i = 1:length(values)
                    if values[i] == source[:w]
                        bbm[(edge2indexdict[key]-1)*Nmodes+i] = Lmean*source[:current]/phi0
                        break
                    end
                end
                break
            end
        end
    end


    # convert from the node basis to the branch basis
    bnm = transpose(Rbnm)*bbm

    Nwtuple=NTuple{length(Nt)+1,Int}(
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
    AoLjnm.nzval .= rand(Complex{Float64},length(AoLjnm.nzval))
    Jsparse = (AoLjnm + invLnm - im.*Gnm*wmodesm - Cnm*wmodes2m)
    # Jsparse = (AoLjnm + invLnm - im*Gnm - Cnm)


    # return (Jsparse,AoLjnm)

    # make the index maps so we can add the sparse matrices together without
    # memory allocations. 
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

    # perform Newton's method with linesearch based on Nocedal and Wright
    # chapter 3 section 5. 
    for n = 1:iterations

        # update the residual function and the Jacobian
        # calcfj!(F,Jsparse,x,wmodesm,wmodes2m,Rbnm,Rbnmt,invLnm,Cnm,Gnm,bnm,
        # Ljb,Ljbm,Nmodes,
        # Nbranches,Lmean,AoLjbmvector,AoLjbm,
        # AoLjnmindexmap,invLnmindexmap,Gnmindexmap,Cnmindexmap)

        calcfj2!(F,Jsparse,x,wmodesm,wmodes2m,Rbnm,Rbnmt,invLnm,Cnm,Gnm,bnm,
        Ljb,Ljbm,Nmodes,
        Nbranches,Lmean,AoLjbm2,AoLjbm,
        AoLjnmindexmap,invLnmindexmap,Gnmindexmap,Cnmindexmap,
        freqindexmap,conjsourceindices,conjtargetindices,phimatrix)

        push!(normF,norm(F))

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
        # dfdalpha = real(dot(F,Jsparse*deltax))
        dfdalpha = real(dot(F,Jsparse,deltax))

        # # evaluate the objective function at Nsample points
        # for alpha in range(0,1,Nsamples)
        #     calcfj!(F,nothing,x - alpha*x1,wmodesm,wmodes2m,Rbnm,Rbnmt,invLnm,Cnm,Gnm,bnm,
        #     Ljb,Ljbm,Nmodes,
        #     Nbranches,Lmean,AoLjbmvector,AoLjbm,
        #     AoLjnmindexmap,invLnmindexmap,Gnmindexmap,Cnmindexmap)    
        #     push!(samples,real(0.5*dot(F,F)))
        # end

        # # evaluate the function at the trial point
        # calcfj!(F,nothing,x+deltax,wmodesm,wmodes2m,Rbnm,Rbnmt,invLnm,Cnm,Gnm,bnm,
        # Ljb,Ljbm,Nmodes,
        # Nbranches,Lmean,AoLjbmvector,AoLjbm,
        # AoLjnmindexmap,invLnmindexmap,Gnmindexmap,Cnmindexmap)


        calcfj2!(F,nothing,x,wmodesm,wmodes2m,Rbnm,Rbnmt,invLnm,Cnm,Gnm,bnm,
        Ljb,Ljbm,Nmodes,
        Nbranches,Lmean,AoLjbm2,AoLjbm,
        AoLjnmindexmap,invLnmindexmap,Gnmindexmap,Cnmindexmap,
        freqindexmap,conjsourceindices,conjtargetindices,phimatrix)

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
            alpha1 = 1
        end
        # if the fitted alpha overshoots the size of the interval (from 0 to 1),
        # then set alpha to 1 and make a full length step. 
        if alpha1 > 1 || alpha1 <= 0
            alpha1 = 1
            f1fit = fp
        end

        # switch to newton once the norm is small enough
        switchofflinesearchtol = 1e-3
        if fp <= switchofflinesearchtol && f <= switchofflinesearchtol && f1fit <= switchofflinesearchtol
            alpha1 = 1
        end
        # alpha1 = 1

        # update x
        x .+= deltax*alpha1

        # x .-= minusdeltax*alpha1

        # push!(alphas,alphafit)
        # push!(fmin,fminfit)
        # push!(fvals,f)
        # push!(fpvals,fp)
        # push!(dfdalphavals,dfdalpha)

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
    # Nports = length(portdict)
    Nports = length(portindices)

    # input = zeros(Complex{Float64},Nports*Nmodes)
    input = Diagonal(zeros(Complex{Float64},Nports*Nmodes))

    output = zeros(Complex{Float64},Nports*Nmodes)
    phibports = zeros(Complex{Float64},Nports*Nmodes)
    inputval = zero(Complex{Float64})
    S = zeros(Complex{Float64},Nports*Nmodes,Nports*Nmodes)

    # phin = out.zero

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
    return NonlinearHB(phin,Rbnm,Ljb,Lb,Ljbm,Nmodes,Nbranches,S)

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

    # println(size(Am))
    # println(nnz(Ljb))
    # println(Nbranches)
    # println(length(Ljb))

    # if size(Am,3) != nnz(Ljb)
    #     throw(DimensionMismatch("The second axis of Am must equal the number of nonzero elements in Ljb (the number of JJs)."))
    # end

    # if length(Ljb) != Nbranches
    #     throw(DimensionMismatch("The length of Ljb should be the same as the number of branches."))
    # end



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
s, conjindicessorted = JosephsonCircuits.calcAoLjbmindices(
    Amatrixindices,
    Ljb,
    Nmodes,
    Nbranches,
    Nfreq);
s

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
Amatrixindices = [1 -2 -3 -4; 2 1 -2 -3; 3 2 1 -2; 4 3 2 1]
Ljb = JosephsonCircuits.SparseArrays.sparsevec([1,3],[1.0,1.0])
Nmodes = 4
Nbranches = length(Ljb)
Nfreq = 4
s, conjindicessorted = JosephsonCircuits.calcAoLjbmindices(
    Amatrixindices,
    Ljb,
    Nmodes,
    Nbranches,
    Nfreq);
for c in conjindicessorted;s.nzval[c] = -s.nzval[c];end;s

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
function calcAoLjbmindices(Amatrixindices::Matrix,Ljb::SparseVector,Nmodes,Nbranches,Nfreq)


    # evaluate Amatrixindices to find the number of entries of each type
    nentries = 0
    nconjentries = 0
    nzeros = 0
    for j in 1:Nmodes
        for k in 1:Nmodes
            if Amatrixindices[j,k] == 0
                nzeros+=1
            else
                nentries+=1
                if Amatrixindices[j,k] < 0
                    nconjentries+=1
                end
            end
        end
    end

    # make these into a sparse matrix. skip any zeros
    conjindices = Vector{Int}(undef,nconjentries*nnz(Ljb))
    I = Vector{Int}(undef,nentries*nnz(Ljb))
    J = Vector{Int}(undef,nentries*nnz(Ljb))
    V = Vector{Int}(undef,nentries*nnz(Ljb))
    Vsort = Vector{Int}(undef,nentries*nnz(Ljb))

    # generate the contents of the sparse matrix 
    n=1
    nconj=1
    for i in 1:nnz(Ljb)
        for j in 1:Nmodes
            for k in 1:Nmodes
                if Amatrixindices[j,k] != 0
                    I[n] = j+(Ljb.nzind[i]-1)*Nmodes
                    J[n] = k+(Ljb.nzind[i]-1)*Nmodes
                    Vsort[n] = n
                    index = abs(Amatrixindices[j,k])+Nfreq*(i-1)
                    V[n] = index
                    if Amatrixindices[j,k] < 0
                        conjindices[nconj] = n
                        nconj+=1
                    end
                    n+=1
                end
            end
        end
    end

    # create the sparse matrix
    s = JosephsonCircuits.SparseArrays.sparse(I,J,Vsort,Nbranches*Nmodes,Nbranches*Nmodes)

    # find the sorting of nzvals in the sparse matrix and apply that same
    # sorting to 
    Vsort2 = copy(s.nzval)
    conjindicessorted = Vsort2[conjindices]

    s.nzval .= V[Vsort2]
    return s, conjindicessorted

end

"""

Update the values in the sparse AoLjbm matrix in place.

"""
function updateAoLjbm2!(AoLjbm::SparseMatrixCSC,Am,Ljb::SparseVector,Lmean,Nmodes,Nbranches,freqindexmap)

    Nfreq = prod(size(Am)[1:end-1])


    # check that there are the right number of nonzero values. 
    # check that the dimensions are consistent with Nmode and Nbranches.

    # if nnz(Ljb)*Nmodes^2 != nnz(AoLjbm)
    #     throw(DimensionError("The number of nonzero elements in AoLjbm are not consistent with nnz(Ljb) and Nmodes."))
    # end

    # if size(Am,2) != nnz(Ljb)
    #     throw(DimensionError("The second axis of Am must equal the number of nonzero elements in Ljb (the number of JJs)."))
    # end

    # if length(Ljb) > Nbranches
    #     throw(DimensionError("The length of Ljb cannot be larger than the number of branches."))
    # end

    # i want a vector length(Ljb) where the indices are the values Ljb.nzind
    # and the values are the indices of Ljb.nzind
    indexconvert = zeros(Int,length(Ljb))
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
