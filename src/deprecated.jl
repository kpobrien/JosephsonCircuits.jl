# Deprecated entry points. Each one warns once and forwards to its
# replacement.

# `connectS(Sa, k, l)` and `connectS(Sa, Sb, k, l)` were split into
# `intraconnectS` and `interconnectS` so that each can also take noise
# covariance matrices: with one name, a second matrix argument would be
# ambiguous between a second scattering matrix and a noise covariance.
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
#         circuitdefs; pumpports = [1], keyword arguments...)
#
# The original `hbsolve` signature: a single pump at the scalar frequency
# `wp`, applied to `pumpports` with the currents `Ip`, four wave mixing
# only, and mode counts given as integers rather than tuples. It is
# translated into a call of the current solvers and warns that it is
# deprecated. Note that the signal modes of the result are ordered as
# `hblinsolve` orders them (the signal at index 1, the rest as listed in
# `modes`), which need not match the order the original solver used.
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

    # the single pump as a one element frequency tuple
    w = (wp,)
    Nharmonics = (2*Npumpmodes,)

    # one source per pump port, all at the pump frequency
    sources = [(mode = (1,), port = pumpports[1], current = Ip[1])]
    @assert length(pumpports) == length(Ip)
    if length(pumpports) > 1
        for i in 2:length(pumpports)
            push!(sources,(mode = (1,), port = pumpports[i], current = Ip[i]))
        end
    end

    # the pump harmonics: odd harmonics only, which is four wave mixing
    freq = removeconjfreqs(
        truncfreqs(
            calcfreqsrdft(Nharmonics),
            dc=false, odd=true, even=false, maxintermodorder=Inf,
        )
    )

    indices = fourierindices(freq)

    Nmodes = length(freq.modes)

    psc = compile(circuit; sorting = sorting)
    # nothing here reads the loops of the circuit graph
    cg = calccircuitgraph(psc; loops = false)
    nm=numericmatrices(psc, cg, circuitdefs, Nmodes = Nmodes)

    nonlinear = hbnlsolve(w, sources, freq, indices, psc, cg, nm;
        iterations = iterations, x0 = nothing, ftol = ftol,
        symfreqvar = symfreqvar, keyedarrays = keyedarrays,
        sensitivitynames = sensitivitynames,
        method = NewtonKrylov(preconditioner = BlockDiagonal(
            factorization = factorization)))

    # the signal modes: the signal and the even pump harmonics on either
    # side of it
    signalfreq =truncfreqs(
        calcfreqsdft((Nsignalmodes,)),
        dc=true,odd=false,even=true,maxintermodorder=Inf,
    )

    # the original solver kept one mode fewer when Nsignalmodes is even;
    # drop the highest one to match
    if mod(Nsignalmodes,2) == 0 && Nsignalmodes > 0
        signalfreq = JosephsonCircuits.removefreqs(
            signalfreq,
            [(Nsignalmodes,)],
        )
    end

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
