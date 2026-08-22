
# Deprecated `connectS` and `connectS!` intra and interconnection functions.
# These are deprecated in order to allow the `intraconnectS` and
# `interconnectS` functions to accept noise correlation matrices without the
# ambiguity of whether the second matrix is a second scattering
# parameter matrix or a noise correlation matrix.
function connectS(Sa::AbstractArray{T,N}, k::Int, l::Int;
    nbatches::Int = Base.Threads.nthreads()) where {T,N}
    Base.depwarn(lazy"connectS(Sa::AbstractArray, k::Int, l::Int)` is deprecated, use `intraconnectS(Sa, k, l)` instead.", :connectS; force=true)
    return intraconnectS(Sa, k, l; nbatches = nbatches)
end

function connectS(Sa::AbstractArray{T,N}, Sb::AbstractArray{T,N}, k::Int, l::Int;
    nbatches::Int = Base.Threads.nthreads()) where {T,N}
    Base.depwarn(lazy"connectS(Sa::AbstractArray, Sb::AbstractArray, k::Int, l::Int)` is deprecated, use `interconnectS(Sa, Sb, k, l)` instead.", :connectS; force=true)
    return interconnectS(Sa, Sb, k, l; nbatches = nbatches)
end

function connectS!(Sout, Sa, k::Int, l::Int;
    nbatches::Int = Base.Threads.nthreads())
    Base.depwarn(lazy"connectS!(Sout, Sa, k::Int, l::Int)` is deprecated, use `intraconnectS!(Sout, Sa, k, l)` instead.", :connectS!; force=true)
    return intraconnectS!(Sout, Sa, k, l; nbatches = nbatches)
end

function connectS!(Sout, Sa, Sb, k::Int, l::Int;
    nbatches::Int = Base.Threads.nthreads())
    Base.depwarn(lazy"connectS!(Sout, Sa, Sb, k::Int, l::Int)` is deprecated, use `interconnectS!(Sout, Sa, Sb, k, l)` instead.", :connectS!; force=true)
    return interconnectS!(Sout, Sa, Sb, k, l; nbatches = nbatches)
end


#     hbsolve(ws, wp, Ip, Nsignalmodes::Int, Npumpmodes::Int, circuit,
#         circuitdefs; pumpports = [1], iterations = 1000, ftol = 1e-8,
#         switchofflinesearchtol = 1e-5, alphamin = 1e-4,
#         symfreqvar = nothing, nbatches = Base.Threads.nthreads(), sorting = :number,
#         returnS = true, returnSnoise = false, returnQE = true, returnCM = true,
#         returnnodeflux = false, returnvoltage = false, returnnodefluxadjoint = false,
#         returnvoltageadjoint = false, keyedarrays::Bool = false,
#         sensitivitynames::Vector{String} = String[], returnSsensitivity = false,
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
    symfreqvar = nothing, nbatches = Base.Threads.nthreads(), sorting = :number,
    returnS::Bool = true, returnSnoise::Bool = false, returnQE::Bool = true,
    returnCM::Bool = true, returnnodeflux::Bool = false,
    returnvoltage::Bool = false, returnnodefluxadjoint::Bool = false,
    returnvoltageadjoint::Bool = false, keyedarrays::Bool = false,
    sensitivitynames::Vector{String} = String[],
    returnSsensitivity::Bool = false, returnZ = nothing,
    returnZadjoint = nothing, returnZsensitivity = nothing,
    returnZsensitivityadjoint = nothing,
    factorization = KLUfactorization())

    Base.depwarn(lazy"""
    Calls the new harmonic balance solvers, [`hbnlsolve`](@ref) and
    [`hblinsolve`](@ref), which work for an arbitrary number of modes and ports),
    using an identical syntax to the legacy harmonic balance solver, which only
    supported four wave mixing processes involving single strong tone and an
    arbitrary number of tone in the linearized solver. This function is
    primarily for testing the new solvers and is now deprecated. Please switch
    to the new syntax.
        """, :hbsolve; force=true)

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
    nonlinear = hbnlsolve(w, sources, freq, indices, psc, cg, nm;
        iterations = iterations, x0 = nothing, ftol = ftol,
        symfreqvar = symfreqvar, keyedarrays = keyedarrays,
        sensitivitynames = sensitivitynames, factorization = factorization)

    # generate the signal modes
    signalfreq =truncfreqs(
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
        returnSsensitivity = returnSsensitivity, returnZ = returnZ,
        returnZadjoint = returnZadjoint, returnZsensitivity = returnZsensitivity,
        returnZsensitivityadjoint = returnZsensitivityadjoint,
        factorization = factorization)

    return HB(nonlinear, linearized)
end
