# =========================================================================
# The nonlinear harmonic balance solve.
# =========================================================================

"""
    hbnlsolve(w::NTuple{N,Number}, Nharmonics::NTuple{N,Int}, sources,
        circuit, circuitdefs; iterations = 1000,
        maxintermodorder = Inf, dc = false, odd = true, even = false,
        x0 = nothing, ftol = 1e-8, symfreqvar = nothing, sorting= :number)

Harmonic balance solver supporting an arbitrary number of large signals
(strong tones or pumps) and arbitrary numbers of ports, sources, and drives
including direct current (zero frequency) or flux pumping using a current
source and a mutual inductor. Use `hblinsolve` to linearize the system of
equations about the operating point found with `hbnlsolve`.

The system is solved in a modified nodal analysis (MNA) formulation in the
node flux basis: resistors with constant real values (including complex
storage with zero imaginary part) and mutually coupled inductor branches
are assigned auxiliary branch current variables with their constitutive
relations kept as explicit equations, which is algebraically equivalent to
the nodal formulation wherever the latter is well posed. Promoting the
coupled inductors keeps the system matrix entries bounded as the coupling
coefficient approaches one, where the nodal inverse inductance entries
diverge as `1/(1-k^2)`. The system is nondimensionalized by the solver
inductance scale `Z0/w0` (see [`calcsolverscale`](@ref)), the geometric
mean port impedance over the geometric mean nonzero drive frequency, so
the residual tolerance `ftol` is independent of the unit system and the
auxiliary variables have magnitudes comparable to the node fluxes. One gauge fixing equation per floating
inductive/Josephson subnetwork and zero-frequency mode makes circuits
which are structurally singular at DC in a purely nodal formulation
(nodes or subnetworks with no inductive path to ground) exactly solvable
without workaround inductors; if the net direct current injected into
such a subnetwork is nonzero, no periodic solution exists and an
`ArgumentError` is thrown. The reported residual norms are those of the
augmented system, a supplied `x0` may have either the node flux length or
the augmented length, and the returned structure contains only the node
fluxes and the original incidence matrix. Commensurate drive frequencies
whose retained intermodulation products reach (numerically) zero
frequency are rejected with an `ArgumentError`. See `src/mna.jl`.
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
$(_DOC_FTOL)
$(_DOC_METHOD)
$(_DOC_STAGEDKWARGS)
$(_DOC_ANDERSON)
- `symfreqvar = nothing`: the symbolic frequency variable, eg `w`.
$(_DOC_SORTING)

# Returns
- `NonlinearHB`: A simple structure to hold the harmonic balance solutions.
    See [`NonlinearHB`](@ref).

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
    circuit, circuitdefs; rtol = 0.0,
    iterations = 1000,
    maxharmonics::NTuple{N,Int} = Nharmonics,
    maxintermodorder = Inf, dc::Bool = false, odd::Bool = true,
    even::Bool = false, x0 = nothing, ftol = 1e-8,
    switchofflinesearchtol = nothing, alphamin = nothing,
    method = :newtonkrylov, andersondepth::Integer = method == :quasinewton ? 5 : 0,
    symfreqvar = nothing, sorting = :number, keyedarrays::Bool = true,
    sensitivitynames::Vector{String} = String[],
    returnoperatingpoint::Bool = false,
    krylovcouplingmodes = :none,
    krylovrecycle::Integer = 0, krylovharvest::Integer = 8,
    krylovkwargs::NamedTuple = (;),
    stagedkwargs::NamedTuple = (;),
    factorization = nothing, backend = CPU(),
    precision::Type{<:AbstractFloat} = Float64, debugJacobian = false,
    linearsolver::AbstractHBLinearSolver = InternalGMRES(),
    returnsystem::Bool = false, assemblejacobian::Bool = true,
    ) where {N}

    if method === :staged
        return stagedhbnlsolve(w, Nharmonics, sources, circuit, circuitdefs;
            iterations = iterations, maxharmonics = maxharmonics,
            maxintermodorder = maxintermodorder, dc = dc, odd = odd,
            even = even, ftol = ftol, symfreqvar = symfreqvar,
            sorting = sorting, keyedarrays = keyedarrays,
            sensitivitynames = sensitivitynames,
            returnoperatingpoint = returnoperatingpoint,
            krylovcouplingmodes = krylovcouplingmodes,
            krylovrecycle = krylovrecycle, krylovharvest = krylovharvest,
            krylovkwargs = krylovkwargs, factorization = factorization,
            backend = backend, precision = precision, stagedkwargs...)
    end

    # calculate the frequency struct
    freq = removeconjfreqs(
        truncfreqs(
            calcfreqsrdft(Nharmonics); dc = dc, odd = odd, even = even,
            maxintermodorder = maxintermodorder, maxharmonics = maxharmonics,
        )
    )

    indices = fourierindices(freq)

    Nmodes = length(freq.modes)

    # parse, graph, and assemble; a typed circuit takes the compiled path
    psc, cg, nm = preparecircuit(circuit, circuitdefs;
        sorting = sorting, Nmodes = Nmodes)


    return hbnlsolve(w, sources, freq, indices, psc, cg, nm;
        rtol = rtol,
        iterations = iterations, x0 = x0, ftol = ftol,
        switchofflinesearchtol = switchofflinesearchtol, alphamin = alphamin,
        method = method, andersondepth = andersondepth,
        symfreqvar = symfreqvar, keyedarrays = keyedarrays,
        sensitivitynames = sensitivitynames,
        returnoperatingpoint = returnoperatingpoint,
        krylovcouplingmodes = krylovcouplingmodes,
        krylovrecycle = krylovrecycle, krylovharvest = krylovharvest,
        krylovkwargs = krylovkwargs,
        factorization = factorization, backend = backend,
        precision = precision, debugJacobian = debugJacobian,
        linearsolver = linearsolver, returnsystem = returnsystem,
        assemblejacobian = assemblejacobian,
        )
end

# A fully numeric circuit (such as the output of a circuit builder) needs no
# component definitions, so `circuitdefs` is optional.
function hbnlsolve(w::NTuple{N,Number}, Nharmonics::NTuple{N,Int}, sources,
    circuit; kwargs...) where {N}
    return hbnlsolve(w, Nharmonics, sources, circuit, Dict{Any,Any}();
        kwargs...)
end

"""
    hbnlsolve(w::NTuple{N,Number}, sources, frequencies::Frequencies{N},
        indices::FourierIndices{N}, psc::CompiledCircuit, cg::CircuitGraph,
        nm::CircuitMatrices; iterations = 1000, x0 = nothing,
        ftol = 1e-8, symfreqvar = nothing)

New version of the nonlinear harmonic balance solver suitable for arbitrary
numbers of ports, sources, and drives including direct current (zero frequency)
or flux pumping using a current source and a mutual inductor.

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
psc = JosephsonCircuits.compile(circuit)
cg = JosephsonCircuits.calccircuitgraph(psc)
nm = JosephsonCircuits.numericmatrices(psc, cg, circuitdefs, Nmodes = Nmodes)

out=hbnlsolve(
    (wp,),
    [
        (mode=(0,),port=1,current=Idc),
        (mode=(1,),port=1,current=Ip),
    ],
    frequencies, fi, psc, cg, nm)
isapprox(out.nodeflux[:],
    ComplexF64[15.190314040027522 - 8.56492651167657e-24im, 2.991103820177504e-6 - 1.8501001011477133e-8im, -6.835392148510984 - 1.0356102442254259e-14im, 7.396422335315908e-6 - 4.5749403967992827e-8im, 6.835392148539885 - 1.0356102451770844e-14im, 1.008026285172782e-5 - 6.23498762664213e-8im],
    atol = 1e-6)

# output
true
```

The system is solved in a modified nodal analysis formulation in the node
flux basis; see the primary [`hbnlsolve`](@ref) docstring and `src/mna.jl`
for details.
"""
function hbnlsolve(w, sources, frequencies::Frequencies,
    indices::FourierIndices, psc::CompiledCircuit, cg::CircuitGraph,
    nm::CircuitMatrices;
    iterations = 1000, x0 = nothing,
    ftol = 1e-8, rtol = 0.0, switchofflinesearchtol = nothing, alphamin = nothing,
    method = :newtonkrylov, andersondepth::Integer = method == :quasinewton ? 5 : 0,
    symfreqvar = nothing, keyedarrays::Bool = true,
    sensitivitynames::Vector{String} = String[],
    returnoperatingpoint::Bool = false,
    krylovcouplingmodes = :none,
    krylovrecycle::Integer = 0, krylovharvest::Integer = 8,
    krylovkwargs::NamedTuple = (;),
    factorization = nothing, backend = CPU(),
    precision::Type{<:AbstractFloat} = Float64, debugJacobian = false,
    linearsolver::AbstractHBLinearSolver = InternalGMRES(),
    returnsystem::Bool = false, assemblejacobian::Bool = true,
    )

    # A solver object carries its own options; a symbol is the legacy
    # spelling. This must happen before setup, which reads `method` to
    # decide whether the real representation entry points are built.
    solverobject = method
    method = solvermethod(method)
    andersondepth = solverobject isa QuasiNewton ?
        solverobject.andersondepth : andersondepth
    krylovkwargs = merge(krylovkwargs, solverkwargs(solverobject))

    # the factorization defaults to the one matching the backend: KLU on the
    # host, cuDSS on a device, where a host factorization could not be applied
    # to the device vectors of the Krylov iteration anyway. An explicit
    # factorization is always honored.
    factorization = if !isnothing(factorization)
        factorization
    elseif backend isa CPU
        KLUfactorization()
    else
        CUDSSFactorization()
    end

    # deprecation warnings for switchofflinesearchtol and alphamin.
    if !isnothing(switchofflinesearchtol)
        Base.depwarn(lazy"The `switchofflinesearchtol` kwarg is deprecated and no longer used (and no longer necessary). Please remove it to avoid errors in future versions.", :hbnlsolve; force=true)
    end

    if !isnothing(alphamin)
        Base.depwarn(lazy"The `alphamin` kwarg is deprecated and no longer used (and no longer necessary). Please remove it to avoid errors in future versions.", :hbnlsolve; force=true)
    end

    # reject non-finite physical inputs before mode construction and
    # source assembly: the frequency canonicalization and source
    # compatibility diagnostics compare against bounds built from these
    # values, and non-finite inputs would otherwise reach them (an
    # infinite value compares as within an infinite bound, and a NaN
    # silently fails every comparison) or fail later without explanation.
    all(isfinite, w) || throw(ArgumentError("All drive frequencies must be finite."))
    for source in sources
        if !isfinite(abs(source[:current]))
            throw(ArgumentError("All source currents must be finite."))
        end
    end

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
    Amatrixconjindices = indices.hbconjmatindices
    # The harmonic balance residual is computed with cyclic Fourier
    # transforms, so its exact Jacobian couples modes whose differences
    # alias back onto the sampled grid; hbmatind with `alias = true`
    # includes those couplings (the sum couplings of hbconjmatind always
    # alias). The exact real Jacobian of method = :newton is assembled from
    # the aliased difference indices, so it is the exact derivative of the
    # residual for multi-tone problems as well; this was established with
    # the matrix-free Jacobian-vector products of HBSystem as the ground
    # truth, which the assembled Jacobian matches to machine precision. The
    # complex holomorphic Jacobian of method = :quasinewton deliberately
    # keeps the truncated (non-aliased) indices: it is an approximation to
    # the exact Jacobian either way, the truncation does not change its
    # iteration counts in practice, and the aliased couplings would densify
    # it and slow its factorization.
    Amatrixindicesaliased = hbmatind(frequencies; alias = true)[2]

    # generate the frequencies of the modes
    Nmodes = length(modes)
    wmodes = calcmodefreqs(w,modes)

    # only the all-zero mode tuple may represent zero frequency: its
    # frequency is exactly 0.0 by construction, which makes the exact
    # iszero classification in the gauge fixing and source compatibility
    # machinery sound. any other retained tuple whose physical frequency
    # cancels to zero, or to within roundoff of the magnitudes being
    # combined (commensurate drive frequencies), would duplicate the DC
    # coordinate with different conjugacy assumptions and produce
    # vanishing capacitor and resistor stamps, so it is rejected with an
    # explanation rather than allowed to form a (nearly) singular system.
    for m in 1:Nmodes
        if any(!iszero, modes[m])
            terms = [float(real(modes[m][j]*w[j])) for j in eachindex(w)]
            if isnumericallyzero(wmodes[m], terms)
                throw(ArgumentError("The mode tuple $(modes[m]) has a "*
                    "(numerically) zero physical frequency for the drive "*
                    "frequencies $(w./(2*pi)) Hz, so it would duplicate "*
                    "the DC coordinate: the drive frequencies are "*
                    "commensurate at the retained intermodulation order. "*
                    "Choose incommensurate drive frequencies (for example "*
                    "by a small offset) or reduce the number of "*
                    "harmonics."))
            end
        end
    end

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
    Cnmcopy = nm.Cnm
    Gnmcopy = nm.Gnm
    invLnmcopy = nm.invLnm
    portindices = nm.portindices
    portnumbers = nm.portnumbers
    portimpedances = nm.portimpedances
    vvn = nm.vvn
    # if there are no inductors, then Lmean will be zero so set it to be one
    Lmean = if iszero(nm.Lmean)
        one(eltype(nm.Lmean))
    else
        nm.Lmean
    end
    # fail immediately, with the actual cause, if any component value
    # contains symbolic variables which were not assigned numerical values
    # in circuitdefs (values depending only on the symbolic frequency
    # variable are frequency dependent components and are accepted).
    checkcomponentvaluesdefined(componentnames, vvn, symfreqvar)

    # nondimensionalize the system with the solver inductance scale Z0/w0
    # (see calcsolverscale) instead of the mean inductance: the scaled
    # entries are then of order one for circuits driven near their
    # characteristic impedance and frequency, the auxiliary branch currents
    # of the modified nodal analysis formulation have magnitudes comparable
    # to the node fluxes even in circuits without inductors, and the
    # residual tolerance ftol is independent of the unit system. The scale
    # multiplies rows only, so the node fluxes and all physical outputs are
    # unchanged. The local name Lmean is kept for continuity with the
    # matrix and plan interfaces.
    Lmean = calcsolverscale(w, componenttypes, vvn, portimpedances, Lmean)
    Lb = nm.Lb

    # find the indices associated with the components for which we will
    # calculate sensitivities
    sensitivityindices = zeros(Int,length(sensitivitynames))
    for i in eachindex(sensitivitynames)
        sensitivityindices[i] = componentnamedict[sensitivitynames[i]]
        if componenttypes[sensitivityindices[i]] == :S
            throw(ArgumentError(lazy"Sensitivities with respect to scattering block components are not supported; got $(sensitivitynames[i])."))
        end
    end

    # calculate the diagonal frequency matrices
    wmodesm = Diagonal(repeat(wmodes, outer = Nnodes-1))
    wmodes2m = Diagonal(repeat(wmodes.^2, outer = Nnodes-1))

    # calculate the source terms in the branch basis
    bbm = calcsources(modes, sources, portindices, portnumbers,
        nodeindices, edge2indexdict, Lmean, Nnodes, Nbranches, Nmodes)

    # convert from the node basis to the branch basis
    bnm = transpose(Rbnm)*bbm

    # calculate the dimensions of the array which holds the frequency
    # domain information for the fourier transform
    Nwtuple = NTuple{length(Nw)+1,Int}((Nw..., length(Ljb.nzval)))

    # create an array to hold the frequency domain data for the
    # fourier transform
    # on the backend: the time domain array and every workspace of the system
    # derive from it with similar, so they follow it there.
    phimatrix = tobackend(backend, zeros(Complex{precision}, Nwtuple))

    # create an array to hold the time domain data for the RFFT. also generate
    # the plans, which come from the package extension on a device backend.
    phimatrixtd, irfftplan, rfftplan = plan_applynl(phimatrix, backend)

    # the number of frequency entries per Josephson junction in phimatrix
    Nfreq = prod(Nwtuple[1:end-1])

    x = if isnothing(x0)
        zeros(Complex{Float64}, (Nnodes-1)*Nmodes)
    else
        copy(x0)
    end
    F = zeros(Complex{Float64}, (Nnodes-1)*Nmodes)
    # substitute in the mode frequencies for components which have frequency
    # defined symbolically.
    Cnm = freqsubst(Cnmcopy, wmodes, symfreqvar)
    Gnm = freqsubst(Gnmcopy, wmodes, symfreqvar)
    # Inductor branches which participate in mutual coupling are promoted
    # to auxiliary branch current variables by the modified nodal analysis
    # formulation below, instead of being eliminated through the inverse of
    # the branch inductance matrix: the inverse inductance entries of the
    # nodal formulation diverge as 1/(1-k^2) as the coupling coefficient k
    # approaches one (and do not exist at |k| = 1), while the branch
    # inductance entries of the promoted constitutive equations remain
    # bounded. numericmatrices already excludes the coupled branches from
    # the inverse inductance matrix.
    Mb = nm.Mb
    coupledbranches = mnacoupledbranches(Mb)
    invLnm = freqsubst(invLnmcopy, wmodes, symfreqvar)


    # take the complex conjugate of the terms associated with modes with
    # negative frequencies. this is the same operation hblinsolve performs
    # with sparseaddconjsubst!. these matrices are used in both the residual
    # and the Jacobian in calcfj2!.
    conjnegfreq!(Cnm, wmodes)
    conjnegfreq!(Gnm, wmodes)
    conjnegfreq!(invLnm, wmodes)

    # scale the matrices for numerical reasons
    rmul!(Cnm,Lmean)
    rmul!(Gnm,Lmean)
    rmul!(invLnm,Lmean)


    # Set up the modified nodal analysis (MNA) formulation, which assigns
    # auxiliary branch current variables to the port resistors, keeping their
    # constitutive relations as explicit equations, and adds one gauge
    # fixing equation per floating component of the static flux-stiffness
    # graph and zero-frequency mode. Eliminating the auxiliary variables
    # recovers the nodal equations exactly, so the results are identical
    # whenever the nodal formulation is well posed. Because the MNA
    # equations are linear, the whole augmentation reduces to padding the
    # incidence matrix with empty columns and folding the constant matrix
    # Amna into the static linear term, after which the plan machinery,
    # HBSystem, and the solvers operate on the augmented system unchanged.
    # See mna.jl for details.
    Rbnmout = Rbnm
    Nnodal = (Nnodes-1)*Nmodes
    # reject inductor and junction values for which the static
    # classification and the matrix stamps are ill defined
    checkstaticstiffnessvalues(componenttypes, vvn)
    # No components are promoted to auxiliary branch currents here: port
    # resistors stay as node conductances (the two forms are algebraically
    # identical, the wave extraction is nodal, and the augmented system is
    # smaller without them). The list remains as the hook for elements a
    # conductance cannot express, such as voltage sources, should they be
    # added; the augmentation machinery (mnaaugmentation, calcAmna) is
    # shared with the gauge equations and stays.
    mnaindices = Int[]
    Nauxr = length(mnaindices)*Nmodes
    Nauxscattering = countscatteringports(psc)*Nmodes
    Naux = Nauxr + length(coupledbranches)*Nmodes + Nauxscattering
    # gauge fixing is based on the connected components of the static
    # flux-stiffness graph (edges are inductors and Josephson junctions),
    # which handles floating inductive and Josephson subnetworks as well
    # as individually isolated nodes. before adding the gauge equations,
    # check that the net direct current injected into each floating
    # component is zero, because otherwise no periodic solution exists and
    # the gauge equation would silently absorb the incompatible source
    # into the flux reference.
    floatingcomponents = calcstaticfluxcomponents(componenttypes,
        nodeindices, vvn, Nnodes)

    # Direct current through resistors. The zero frequency mode of a
    # periodic flux carries no voltage, so a resistor is an open circuit at
    # DC. The missing coordinate is the average voltage, constant on each
    # static flux component because an inductor or zero-voltage junction is
    # a short there. For linear conductances and prescribed current sources
    # that block does not depend on the nonlinear state, so it is solved
    # here and its resistor current subtracted from the zero frequency
    # source: an exact elimination. `bnmsource` is kept because the port and
    # scattering outputs need the current which was applied, not the
    # corrected one.
    # The average voltages are unknowns with their own transport rows, and
    # the resistor and block currents they drive reach the nodes through the
    # coupling, so the system is built on the applied source: nothing is
    # subtracted from it here and so nothing can be subtracted twice.
    # a copy, not an alias: `bnm` is handed to the `HBSystem`, which owns it
    # from then on and writes into it when the drive is scaled. The direct
    # current rows are built from the drive as applied, so they must not
    # move when that happens.
    bnmsource = copy(bnm)
    dcplan = dcconductanceplan(floatingcomponents, Gnm, wmodes, Nmodes,
        Nnodes)

    # Nothing injects direct current in most circuits, and then every
    # average voltage and every block port current is zero: the explicit
    # block would carry a subsystem whose answer is the zero it starts at
    # and charge for it on every residual, product and preconditioner
    # application. Skipping it leaves the `i = 0` rows the scattering stamp
    # writes, which are then the right answer rather than a simplification.
    dcexplicit = !isnothing(dcplan) && dcinjected(dcplan, bnmsource, Nmodes)

    # The average voltages, their transport rows and the blocks' zero
    # frequency rows live in a wrapper around the system rather than in it,
    # so the raw `HBSystem` is not the system being solved: its scattering
    # rows still say `i = 0` and its resistors still carry no direct
    # current. Handing that object out, or differentiating it, would define
    # a second and different operating point. Refuse instead, until the
    # composite system exists to hand out.
    if dcexplicit
    end
    # the canonical Jacobian is assembled as a host `SparseMatrixCSC`; on a
    # device backend `Jr` is not one, and there is no device assembly yet
    if dcexplicit && method == :newton && !(backend isa CPU)
        throw(ArgumentError("a circuit which injects direct current is solved with an assembled canonical Jacobian under method = :newton, and that assembly is host only. Use method = :newtonkrylov on this backend, which is matrix free, or run :newton on the CPU."))
    end
    if dcexplicit && !(method in (:newtonkrylov, :newton, :external))
        throw(ArgumentError(lazy"a circuit which injects direct current is solved with the average voltages as unknowns, which needs the real system; $(method) solves the complex holomorphic one. Use method = :newtonkrylov, :newton, or an ExternalSolver."))
    end
    dcsol = nothing
    # the converged canonical state, which is the point the whole system is
    # differentiated at when the direct current block is active
    dccanonical = nothing

    gaugeindices = calcdcgaugeindices(floatingcomponents, wmodes, Nmodes)
    Amna = calcAmna(mnaindices, nodeindices, vvn, gaugeindices, wmodes,
        Nmodes, Nnodes, Lmean)
    # pad with the coupled inductor auxiliary variables and add their
    # constitutive equations and Kirchhoff current law couplings, which are
    # real and frequency independent (see calcAmnaind).
    Amna = mnapad(Amna, length(coupledbranches)*Nmodes + Nauxscattering)
    AmnaL = calcAmnaind(coupledbranches, Lb, Mb, cg.Rbn, Nmodes,
        Nnodal + Nauxr, Nnodal + Naux, Lmean)
    Amna = spaddkeepzeros(Amna, AmnaL)
    # the hybrid (wave to modified nodal analysis) contribution of the
    # scattering block components: the pump mode frequencies are fixed, so
    # the constitutive equations im*w_m*Lmean*B(w_m)*phi - C(w_m)*i = 0 and
    # the Kirchhoff current law couplings of the auxiliary port currents
    # form a constant matrix, folded into the augmentation exactly like the
    # promoted resistor equations (same signed frequency conjugation
    # convention and Lmean scaling; see scatteringlinearterm). The
    # residual, Jacobian, and solver machinery operate on the augmented
    # system unchanged.
    # the stamped blocks, kept because the explicit direct current path
    # needs each block's auxiliary base to find its zero frequency current
    stampedblocks = StampedScatteringBlock[]
    if Nauxscattering > 0
        # The pump mode frequencies are fixed, so this is a constant matrix
        # folded into the augmentation before the solve begins, exactly as
        # the promoted resistor equations are. Nothing about it is per
        # iteration, so it needs nothing of the backend: the augmented
        # system it produces is what every backend already solves.
        ssys = scatteringstampsystem(psc.scatteringblocks, Nmodes;
            auxoffset = Nnodal + Naux - Nauxscattering,
            Ntotal = Nnodal + Naux, scale = Lmean)
        append!(stampedblocks, ssys.blocks)
        Amna = spaddkeepzeros(Amna, scatteringlinearterm(psc, wmodes,
            Nmodes; auxoffset = Nnodal + Naux - Nauxscattering,
            Ntotal = Nnodal + Naux, scale = Lmean,
            blocks = psc.scatteringblocks))
    end
    # pad the system matrices and vectors with the auxiliary variables.
    # the incidence matrix gains empty columns so the branch fluxes and
    # the Josephson junction terms are unaffected, and the constant MNA
    # equations are folded into the static linear term, which is added
    # with unit coefficient and is its own contribution to the Jacobian.
    Rbnm = hcat(Rbnm, spzeros(eltype(Rbnm), size(Rbnm,1), Naux))
    invLnm = spaddkeepzeros(mnapad(invLnm, Naux), Amna)
    Cnm = mnapad(Cnm, Naux)
    Gnm = mnapad(Gnm, Naux)
    bnm = vcat(bnm, zeros(eltype(bnm), Naux))
    wmodesm = Diagonal(repeat(wmodes,
        outer = Nnodes-1+length(mnaindices)+length(coupledbranches)+
            countscatteringports(psc)))
    wmodes2m = Diagonal(repeat(wmodes.^2,
        outer = Nnodes-1+length(mnaindices)+length(coupledbranches)+
            countscatteringports(psc)))
    if length(x) == Nnodal
        # accept a nodal initial guess, materializing keyed arrays or
        # other array types into a plain vector. transform it into the
        # selected gauge (a physically irrelevant common DC shift of a
        # floating component would otherwise enter the gauge fixing rows)
        # and initialize the auxiliary currents consistently with the node
        # fluxes, so the initial augmented residual equals the initial
        # nodal Kirchhoff current law residual.
        x = vcat(Vector{Complex{Float64}}(vec(x)),
            zeros(Complex{Float64}, Naux))
        mnagaugenormalize!(x, floatingcomponents, wmodes, Nmodes)
        mnainitialauxall!(x, Amna, Nnodal, Nauxr, coupledbranches, Lb, Mb,
            cg.Rbn, Nmodes, Lmean)
    elseif length(x) == Nnodal + Naux
        # accept a full augmented state, with the layout documented in
        # calcAmna, applying the same gauge normalization and reconciling
        # the auxiliary currents with the constitutive relations.
        x = Vector{Complex{Float64}}(vec(x))
        mnagaugenormalize!(x, floatingcomponents, wmodes, Nmodes)
        mnainitialauxall!(x, Amna, Nnodal, Nauxr, coupledbranches, Lb, Mb,
            cg.Rbn, Nmodes, Lmean)
    else
        throw(DimensionMismatch(lazy"The initial value x0 has length $(length(x)) but the solver expects $(Nnodal) node flux unknowns, optionally followed by $(Naux) auxiliary current unknowns."))
    end
    F = zeros(Complex{Float64}, Nnodal + Naux)

    modelayout = ModeLayout(selfconjmodes(frequencies), Nnodal + Naux)
    Fr = complex_to_real(F, modelayout.isreal)
    xr = complex_to_real(x, modelayout.isreal)

    # Build the Jacobian matrix and precomputed assembly plan for the chosen
    # method. The plans map the Fourier coefficients of cos(phi(t)) and the
    # frequency dependent linear term matrices directly into the nonzeros of
    # the Jacobian, so no complex branch matrices, incidence matrix
    # multiplications, or real conversions are performed while iterating.
    # Only the machinery for the chosen method is built, unless the debug
    # output is requested, in which case both are built so they can be
    # compared.
    # the complex Jacobian, as :newton's real one below: on a backend the
    # structure is stored transposed, because a device factorization is
    # compressed by rows, and the values live on the backend too. The host
    # `Jx` is still built either way, because `debugJacobian` compares
    # against it and the sensitivity calculation reads its structure.
    Jx, complexjacobianplan = if method == :quasinewton || debugJacobian
        plancomplexjacobian(Amatrixindices, Ljb, Lmean, Rbnm, Nmodes,
            Nbranches, Nfreq, invLnm, Gnm, Cnm)
    else
        nothing, nothing
    end
    devicex = method == :quasinewton && !(backend isa CPU)
    Jxb, complexjacobianplan = if devicex
        Jxt = sparse(transpose(Jx))
        nodesandsignsx = branchnodesandsigns(Rbnm, Nmodes, Nbranches)
        cjp = planstructurecomplexjacobian(Jxt, Float64, Amatrixindices, Ljb,
            Lmean, nodesandsignsx, invLnm, Gnm, Cnm, wmodesm, wmodes2m,
            Nmodes, Nfreq, backend; transposed = true)
        # the pattern goes to the backend as well, so the sparse products of
        # the line search read it there rather than from the host
        patt = DeviceSparsePattern(
            tobackend(backend, SparseArrays.getcolptr(Jxt)),
            tobackend(backend, rowvals(Jxt)), size(Jxt, 1), size(Jxt, 2))
        DeviceValuedSparseMatrix(patt,
            tobackend(backend, zeros(eltype(Jx), nnz(Jxt)))), cjp
    else
        Jx, complexjacobianplan
    end

    # method = :newtonkrylov deliberately does not appear here. Its steps come
    # from the matrix-free Jacobian-vector product and its preconditioner
    # builds and assembles its own, much sparser, restricted plan, so neither
    # the full Jacobian nor its plan is ever used. On a multi-tone problem
    # that plan is the single largest object in the solve (millions of stored
    # entries), so not building it is a substantial saving in both time and
    # memory.
    Jr, realjacobianplan = if method == :newton || debugJacobian ||
            ((returnsystem || method == :external) && assemblejacobian) ||
            returnoperatingpoint
        # on a backend the structure is built there, transposed, because a
        # device factorization is compressed by rows, and the assembly writes
        # its values there too; on a host it is the matrix KLU factorizes
        devicej = !(backend isa CPU)
        Jrs, nodesandsigns = realjacobianstructure(Amatrixindicesaliased,
            Amatrixconjindices, Ljb, Rbnm, Nmodes, Nbranches, invLnm, Gnm,
            Cnm, modelayout, modelayout, Float64;
            transposed = devicej, backend = backend)
        rjp = planstructurerealjacobian(Jrs, Float64, Amatrixindicesaliased,
            Amatrixconjindices, Ljb, Lmean, nodesandsigns, invLnm, Gnm, Cnm,
            wmodesm, wmodes2m, modelayout, modelayout, Nmodes, Nfreq,
            backend; transposed = devicej)
        (devicej ? DeviceValuedSparseMatrix(Jrs,
            tobackend(backend, zeros(Float64, nnz(Jrs)))) : Jrs), rjp
    else
        nothing, nothing
    end

    # the unified evaluation object for the nonlinear system: the residual,
    # the matrix-free Jacobian-vector and Hessian-vector products, and the
    # assembled Jacobians are all evaluated through it, in the complex or
    # the equivalent real representation.
    # the real form of the linear term is needed exactly when the solve
    # applies the real representation entry points. That is stated directly
    # rather than derived from realjacobianplan: :newtonkrylov solves in the
    # real representation but takes its preconditioner from a
    # ModeCouplingPreconditioner, which builds its own restricted plan, so it
    # needs no full realjacobianplan and the two conditions are not the same.
    # method = :quasinewton stays in the complex representation and skips it.
    # an external solver and the `returnsystem` path both work in the real
    # representation, because the residual is not complex differentiable
    realrepresentation = method == :newton || method == :newtonkrylov ||
        method == :external || debugJacobian || returnoperatingpoint ||
        returnsystem
    sys = HBSystem(Rbnm, invLnm, Gnm, Cnm, wmodesm, wmodes2m, bnm,
        Ljb, Ljbm, Lmean, Nbranches, freqindexmap, conjsourceindices,
        conjtargetindices, phimatrix, phimatrixtd, irfftplan, rfftplan,
        modelayout, realjacobianplan, complexjacobianplan, backend;
        realbackward = realrepresentation)

    # the residual and the complex (holomorphic) Jacobian, an approximation
    # to the exact Jacobian, for method == :quasinewton
    function fj!(F, Jx, x)
        setpoint!(sys, x)
        isnothing(F) || residual!(F, sys)
        isnothing(Jx) || jacobian!(Jx, sys)
        return nothing
    end

    # the residual and the exact Jacobian of the equivalent real system, for
    # method == :newton
    function fjreal!(Fr, Jr, xr)
        setpoint!(sys, xr)
        isnothing(Fr) || residual!(Fr, sys)
        isnothing(Jr) || jacobian!(Jr, sys)
        # write back the canonical real representation of the point. through
        # the plan rather than complex_to_real!, whose serial cursor over a
        # BitVector is host only; the two agree on CPU().
        applycomplextoreal!(xr, sys.nonlineartermplan, sys.x)
        return nothing
    end

    # diagnostics for each invocation of the nonlinear solver, stored in
    # the output so the solution process can be assessed after the run
    solverstages = IterationInfo[]

    # the source scale at which the source stepping detected a fold in
    # the branch of solutions, if any
    sourcefold = Ref(NaN)

    # use this for debugging purposes to return the residual and Jacobian
    # functions along with the ingredients from which the Jacobians are
    # assembled, so reference implementations can be constructed, eg. in the
    # tests.
    # The canonical layout, built once for whichever method uses it. Both
    # the matrix free path and the direct one solve the same system in the
    # same coordinates; only the way they apply the Jacobian differs.
    canonwork = if dcexplicit
        tr = dcexplicit ? transportrows(dcplan, bnmsource, Nmodes) : nothing
        Lc = tobackend(backend, compositelayout(modelayout, frequencies.modes;
            nvdc = isnothing(tr) ? 0 : nvoltages(tr)))
        # a block's zero frequency row is `i = 0` in the stamp; with an
        # average voltage to respond to it becomes the block's own relation,
        # which is what makes it visible at direct current
        br = (dcexplicit && Nauxscattering > 0) ?
            dcblockrows(stampedblocks, dcplan.componentof, Nmodes,
                dcplan.modeindex, Lc.nac, Nnodes - 1, Lmean) : nothing
        # the zero frequency block holds one entry per node and then one per
        # auxiliary unknown; only the nodal ones carry the coupling.
        #
        # Not named `w`: that is this function's drive frequency, and
        # assigning to it here rebound the argument, so every circuit which
        # injected direct current reported the work object as its drive and
        # `hbsolve` could not compute the mode frequencies from it.
        CanonicalWork(Lc, tobackend(backend,
                convert(Vector{precision}, xr));
            transport = tr, blockrows = br, nnodaldc = Nnodes - 1)
    else
        nothing
    end

    if returnsystem
        # Everything an external solver needs, without solving: the
        # evaluation object, the initial value, the real representation
        # layout, and the assembled real Jacobian only if it was asked for
        # (on a multi-tone problem its plan is the largest object in the
        # solve, and exactly what a matrix-free solver avoids).
        return (sys=sys, xr=xr, Fr=Fr, modelayout=modelayout,
            Jr=(assemblejacobian ? Jr : nothing), Nnodal=Nnodal,
            dcplan=dcplan, dcsol=dcsol, bnmsource=bnmsource,
            # with direct current active, `sys` alone is not the system:
            # the average voltages and the blocks' zero frequency rows are
            # in here, and a caller must apply them
            canonicalwork=canonwork, dcexplicit=dcexplicit,
            frequencies=frequencies,
            Amatrixindicesaliased=Amatrixindicesaliased,
            Amatrixconjindices=Amatrixconjindices, Ljb=Ljb, Lmean=Lmean,
            Rbnm=Rbnm, Nmodes=Nmodes, Nbranches=Nbranches, Nfreq=Nfreq,
            invLnm=invLnm, Gnm=Gnm, Cnm=Cnm, Amatrixmodes=Amatrixmodes)
    end
    if debugJacobian
        return (F=F, x=x, Fr=Fr, xr=xr, Jx=Jxb, Jr=Jr, fj=fj!, fjreal=fjreal!,
            canonicalwork=canonwork, dcexplicit=dcexplicit,
            sys=sys, Nnodal=Nnodal, mnaindices=mnaindices,
            gaugeindices=gaugeindices, floatingcomponents=floatingcomponents,
            coupledbranches=coupledbranches,
            complexjacobianplan=complexjacobianplan,
            realjacobianplan=realjacobianplan, phimatrix=phimatrix,
            cosphimatrix=(x -> (setpoint!(sys, x);
                JosephsonCircuits._updatecosphimatrix!(sys);
                sys.phimatrix)), modelayout=modelayout,
            Amatrixindices=Amatrixindices, Amatrixmodes=Amatrixmodes,
            modes=modes,
            Amatrixindicesaliased=Amatrixindicesaliased,
            Amatrixconjindices=Amatrixconjindices, Ljb=Ljb, Ljbm=Ljbm,
            Lmean=Lmean, Rbnm=Rbnm, invLnm=invLnm, Gnm=Gnm, Cnm=Cnm,
            wmodesm=wmodesm, wmodes2m=wmodes2m, Nmodes=Nmodes,
            Nbranches=Nbranches, Nfreq=Nfreq)
    end

    # The residual is a sum of terms whose size is set by the applied
    # source, so no iteration can push it below the rounding error of that
    # sum. The nondimensionalizing scale (see `calcsolverscale`) is read off
    # the port reference impedances, and a circuit whose interior sits at a
    # very different impedance -- a hundred ohm bridge between two ports
    # made nearly open at a gigaohm, say -- is left with a scaled source
    # many orders above one and a residual floor above the default absolute
    # `ftol`. Asking for less than the floor is asking the iteration to
    # converge on noise, and whether it reports success is then decided by
    # which way the last rounding happened to fall. Raise the tolerance to
    # the floor when the floor is the larger. The floor is the `sqrt(n)*eps`
    # accumulation over the terms of the sum with room to spare, since what
    # is wanted is a tolerance the arithmetic can reach and not a sharp
    # bound on the rounding. For a circuit driven near its characteristic
    # impedance it is orders of magnitude below `ftol` and nothing changes.
    ftol = max(ftol,
        16*sqrt(length(xr))*eps(real(eltype(xr)))*norm(bnmsource))

    # solve the nonlinear system
    # The canonical layout is wired into the matrix free path, which is
    # where the scatter and the gather sit in the inner loop and where the
    # cost of them has to be known. The assembled Jacobian methods would
    # additionally need `P J Pᵀ`, which stage 5 brings with the transport
    # rows; until then say so rather than ignore the request.
    info = if method == :quasinewton

        solveonbackend!(fj!, F, Jxb, x, backend; iterations = iterations,
            ftol = ftol, rtol = rtol, andersondepth = andersondepth,
            factorization = factorization)

    elseif method == :newton

        # solve the equivalent real nonlinear system with the exact real
        # Jacobian then convert back to complex
        info = if dcexplicit
            Lc = canonwork.layout
            nc = canonicaldim(Lc)
            # the pattern, from one assembly at the starting point
            fjreal!(nothing, Jr, xr)
            jplan = canonicaljacobianplan(Jr, canonwork)
            Jc = jplan.J
            uc = zeros(eltype(xr), nc)
            Fc = zeros(eltype(Fr), nc)
            gathercanonical!(uc, xr, Lc)
            out = solveonbackend!(
                canonicalfj(fjreal!, canonwork, Jr, jplan), Fc, Jc, uc, backend;
                iterations = iterations, ftol = ftol, rtol = rtol,
                andersondepth = andersondepth, factorization = factorization)
            scattercanonical!(xr, uc, Lc)
            scattercanonical!(Fr, Fc, Lc)
            if dcexplicit
                dccanonical = Array(uc)
                dcsol = dcsolutionfrom(dcplan,
                    Array(view(uc, voltagerange(Lc))))
            end
            out
        else
            solveonbackend!(fjreal!, Fr, Jr, xr, backend;
                iterations = iterations, ftol = ftol, rtol = rtol,
                andersondepth = andersondepth, factorization = factorization)
        end
        real_to_complex!(x,xr,modelayout.isreal)
        real_to_complex!(F,Fr,modelayout.isreal)
        info
    elseif method == :newtonkrylov

        # the matrix free real Jacobian
        jvpreal!(Jvr, vr) = jacobianvectorproduct!(Jvr, sys, vr)

        # The right preconditioner is the Jacobian materialized on a
        # coarser frequency grid: the mode coupling is kept in full only on
        # the modes where the linear circuit responds strongly (measured,
        # see modesoftness) and reduced to the mode diagonal elsewhere, so
        # its factorization is a fraction of the full Jacobian's. The setup
        # above was done once, on the user's harmonic selection; the coarse
        # grid exists only inside the preconditioner, so the residual, the
        # Jacobian-vector product and the solution are always those of the
        # requested truncation.
        # `krylovrecycle > 0` wraps the base preconditioner in a recycled
        # deflation subspace (see RecyclingPreconditioner). The intended
        # pairing is with `krylovcouplingmodes = :none`: a batch of small
        # per mode factorizations plus dense level 2 BLAS, which needs no
        # large sparse factorization and is the part of this solver that
        # ports to a GPU.
        #
        # The default is the mode block diagonal, whose factorization is a
        # batch of small independent per mode solves rather than one large
        # sparse factorization. On a strongly pumped line the block diagonal
        # alone stalls; what rescues it is escalation, which grows the base
        # only on repeated linear failures and in practice fires once or
        # twice.
        #
        # `krylovcouplingmodes = :band => p` restricts the retained coupling by
        # harmonic *offset* rather than by column (see `modebandmask`). That is
        # the restriction the Toeplitz structure of the nonlinear term asks for,
        # and at equal fill it is not close: on an eight mode chain driven to
        # max|phi| = 1.9 rad a bandwidth of one converges a Newton path in 118
        # GMRES iterations where a two column selection storing the same number
        # of nonzeros fails to converge in 1051. Its fill grows linearly in the
        # mode count where the full Jacobian's grows quadratically, so it
        # overtakes the full factorization once there are enough modes: measured
        # on a fixed circuit at max|phi| ~ 1.75 rad, the ratio of banded to full
        # total solve time is 1.28 at 8 modes, 0.75 at 16, 0.61 at 24 and 0.31
        # at 32, while the bandwidth that wins stays at one. It escalates by one
        # offset at a time rather than jumping to the full operator. Measured on a two tone line, this is the only configuration
        # whose standing improves with problem size: at 288 cells it is 3.72 s
        # against 8.19 s for mode selection and 5.48 s for the direct solve,
        # having been the slower of the two at 128 cells.
        #
        # `krylovrecycle > 0` additionally wraps the base in a recycled
        # deflation subspace. It also rescues the block diagonal, and needs no
        # sparse factorization at all, which makes it the natural fit for a
        # GPU; but it is a second answer to the same problem and measured a
        # net loss against escalation at 192 cells and above, so it is off by
        # default. `krylovcouplingmodes = 12, krylovrecycle = 0` recovers the
        # earlier frequency based mode selection, which is still the fastest
        # option on moderately sized lines.
        base = ModeCouplingPreconditioner(sys, Amatrixindicesaliased,
            Amatrixconjindices, Ljb, Lmean, Rbnm, Nmodes, Nbranches, Nfreq,
            invLnm, Gnm, Cnm, modelayout;
            couplingmodes = krylovcouplingmodes,
            factorization = factorization, precision = precision,
            Amatrixmodes = Amatrixmodes)

        # the Krylov workspaces are allocated `similar` to the vectors handed
        # in, so handing in device vectors is what puts the whole iteration on
        # the device. On CPU() tobackend adopts them and these are no-ops.
        xrb = tobackend(backend, convert(Vector{precision}, xr))
        Frb = tobackend(backend, convert(Vector{precision}, Fr))

        # built from the vector rather than from its length, so the deflation
        # subspace lives where the iteration does. The `n::Integer` method
        # allocates on the host, which on a device backend leaves the whole
        # recycling path handing host arrays to device kernels.
        pc = if krylovrecycle > 0
            RecyclingPreconditioner(base, jvpreal!, xrb;
                kmax = krylovrecycle, kharvest = krylovharvest)
        else
            base
        end

        # Defaults tuned for this preconditioner, overridable through
        # `krylovkwargs`. The preconditioner is refreshed eagerly, which the
        # generic default does not do: rebuilding the block diagonal costs a
        # fraction of a full factorization, while a preconditioner frozen at
        # zero flux is stale immediately, because at zero flux the Jacobian
        # has no harmonic coupling at all.
        #
        # The restart length cannot be set independently of the
        # preconditioner. A restricted but not diagonal preconditioner leaves
        # a small number of directions which the Krylov space must resolve
        # before a restart discards the progress made on them, and until it
        # can, the iteration count is pinned and does not respond to the
        # preconditioner at all. Measured on a two tone grid of 74 modes, a
        # banded preconditioner needs 1441 iterations without converging at a
        # restart of 30, 60 and 120, and converges in 216 at 240; the block
        # diagonal converges at none of them.
        #
        # the tuned solver defaults (long restart cycle, eager refresh)
        # live on `nlsolvekrylov!` itself, the single source of truth, so a
        # direct caller gets the same solver; `krylovkwargs` overrides any
        # of them
        info = if dcexplicit
            # The canonical layout groups the state by role rather than by
            # node, which is what stage 5 needs and what this measures. It
            # is a permutation here, so the iteration is the same one in a
            # rotated basis and must reach the same point; the benchmark in
            # `bench/compositelayout.jl` compares the two.
            work = canonwork
            L = work.layout
            ucb = similar(xrb, canonicaldim(L))
            Fcb = similar(Frb, canonicaldim(L))
            # the explicit voltages start at zero, which is the answer when
            # nothing injects direct current; `gathercanonical!` writes the
            # flux blocks only and would leave them undefined
            fill!(ucb, zero(eltype(ucb)))
            gathercanonical!(ucb, xrb, L)
            out = nlsolvekrylov!(canonicalresidual(fjreal!, work),
                canonicaljvp(jvpreal!, work), Fcb, ucb,
                CanonicalPreconditioner(pc, work);
                iterations = iterations, ftol = ftol, rtol = rtol,
                linearsolver = linearsolver,
                krylovkwargs...)
            scattercanonical!(xrb, ucb, L)
            scattercanonical!(Frb, Fcb, L)
            # report the voltages this path actually solved for, rather than
            # the eliminated ones. They agree while the direct current
            # devices are linear conductances, which the tests assert; they
            # will not once a block's port current is an unknown too.
            if dcexplicit
                dccanonical = Array(ucb)
                dcsol = dcsolutionfrom(dcplan,
                    Array(view(ucb, voltagerange(L))))
            end
            out
        else
            nlsolvekrylov!(fjreal!, jvpreal!, Frb, xrb, pc;
                iterations = iterations, ftol = ftol, rtol = rtol,
                linearsolver = linearsolver,
                krylovkwargs...)
        end
        # back to the host for the complex representation returned to the
        # caller, whose conversion walks a BitVector serially
        copyto!(xr, tohost(xrb))
        copyto!(Fr, tohost(Frb))
        real_to_complex!(x,xr,modelayout.isreal)
        real_to_complex!(F,Fr,modelayout.isreal)
        info
    elseif method == :external
        # hand the caller's root finder the system as a problem object; it
        # needs nothing from this scope beyond what the problem carries.
        # With direct current the problem carries the augmentation too, so
        # the caller solves the system which was posed rather than the
        # harmonic part of it, in the canonical unknowns.
        aug = dcexplicit ? DCAugmentation(canonwork, Jr) : nothing
        parts = (Amatrixindicesaliased = Amatrixindicesaliased,
             Amatrixconjindices = Amatrixconjindices, Ljb = Ljb,
             Lmean = Lmean, Rbnm = Rbnm, Nmodes = Nmodes,
             Nbranches = Nbranches, Nfreq = Nfreq, invLnm = invLnm,
             Gnm = Gnm, Cnm = Cnm, Amatrixmodes = Amatrixmodes)
        u0 = if dcexplicit
            v = zeros(Float64, canonicaldim(canonwork.layout))
            gathercanonical!(v, xr, canonwork.layout)
            v
        else
            copy(xr)
        end
        prob = HBNonlinearProblem(sys, modelayout, u0,
            dcexplicit ? (isnothing(aug.jplan) ? nothing : aug.jplan.J) : Jr,
            parts; augmentation = aug)
        u, extconverged = solverobject.f(prob, copy(u0))
        Fc = similar(u)
        hbresidual!(Fc, prob, u)
        if dcexplicit
            scattercanonical!(xr, u, canonwork.layout)
            scattercanonical!(Fr, Fc, canonwork.layout)
            dccanonical = Array(u)
            dcsol = dcsolutionfrom(dcplan,
                Array(view(u, voltagerange(canonwork.layout))))
        else
            copyto!(xr, u); copyto!(Fr, Fc)
        end
        real_to_complex!(x, xr, modelayout.isreal)
        real_to_complex!(F, Fr, modelayout.isreal)
        IterationInfo("external", 1.0, 0.0, extconverged, 0,
            [norm(Fc)], Float64[], Int[], Bool[], [])
    else
        throw(ArgumentError("Method $(method) is not defined."))
    end


    push!(solverstages, IterationInfo(info.label, 1.0,
        info.regularization, info.converged, info.iterations,
        info.normresidual, info.alpha, info.backtracks,
        info.andersonaccepted, info.krylov))

    if !info.converged
        @warn string(lazy"Solver did not converge.")
    end

    converged = info.converged

    # validate the original, ungauged Kirchhoff current law equations
    # directly by reconstructing their residuals from the augmented
    # residual and the state (the gauge fixing equations add x[g] to the
    # augmented residual of each gauge row g). the acceptance policy,
    # including its block-relative infinity-norm tolerance and its
    # rejection of non-finite residuals, is implemented and unit tested in
    # mnavalidatekcl.
    if converged && !isempty(gaugeindices)
        setpoint!(sys, x)
        # `sys` holds its state wherever the backend put it, while `F` and `x`
        # are host vectors, so the residual is evaluated into a buffer like the
        # system's own state and copied back for the host side validation.
        # Evaluating straight into `F` launches the backward term kernel with a
        # host array argument, which does not compile.
        Fgauge = similar(sys.x)
        residual!(Fgauge, sys)
        copyto!(F, Fgauge)
        # On the explicit path the system carries the applied source and the
        # resistor current arrives through the transport coupling, which
        # lives outside it. Adding it back, and correcting the source the
        # same way, hands this check the pair the elimination would have
        # produced, so it validates the physics rather than the path.
        bnmkcl = bnm
        if dcexplicit && !isnothing(dcsol)
            for p in 2:Nnodes
                F[(p-2)*Nmodes + dcplan.modeindex] +=
                    dcsol.scaledcurrent[p-1]
            end
            bnmkcl = applydcconductance(bnmsource, dcplan, dcsol, Nmodes)
        end
        kclok, normkcl, kcltol = mnavalidatekcl(F, x, gaugeindices,
            Nnodal, bnmkcl, ftol)
        if !kclok
            @warn "The original (ungauged) Kirchhoff current law equations are violated beyond the solver resolution at the returned solution, indicating that a gauge fixing equation absorbed an incompatibility, such as a net direct current injected into a floating subnetwork, into the arbitrary flux reference. Marking the solution as not converged." normkcl kcltol
            converged = false
        end
    end

    # The static flux partition treats a junction as a short at zero
    # frequency, which is a statement about the junction being in the zero
    # voltage state. The point is converged and the sine of the branch flux
    # is cached there, so the direct current each junction carries is a mean
    # over the reconstructed period and costs nothing to look at.
    if converged && !isempty(Ljb.nzind)
        setpoint!(sys, x)
        residual!(similar(sys.x), sys)
        branchnames = Dict{Int,Vector{String}}()
        for i in eachindex(componenttypes)
            componenttypes[i] === :Lj || continue
            key = (nodeindices[1,i], nodeindices[2,i])
            b = get(edge2indexdict, key, 0)
            iszero(b) || push!(get!(Vector{String}, branchnames, b),
                String(componentnames[i]))
        end
        checkjunctiondc(tohost(sys.sintd), Ljb.nzind, branchnames)
    end

    # remove the auxiliary variables so the output contains only the node
    # fluxes
    nodeflux = x[1:Nnodal]

    # the norm of the residual at the initial value, for the diagnostics
    normF0 = info.normresidual[1]
    # the norm of the residual at the returned solution, for the
    # diagnostics
    normFfinal = info.normresidual[end]


    # assemble the stored diagnostics of the solution process
    solverinfo = SolverInfo(solverstages, normF0, normFfinal, converged,
        sourcefold[])

    # calculate the scattering parameters for the pump
    Nports = length(portindices)
    S = zeros(Complex{Float64}, Nports*Nmodes, Nports*Nmodes)
    inputwave = zeros(Complex{Float64}, Nports*Nmodes)
    outputwave = zeros(Complex{Float64}, Nports*Nmodes)
    if !isempty(S)
        calcinputoutput!(inputwave, outputwave, nodeflux,
            bnm[1:Nnodal]/Lmean,
            portindices, portindices, portimpedances,
            portimpedances, nodeindices, componenttypes, wmodes, symfreqvar)
        calcscatteringmatrix!(S, inputwave, outputwave)
    end

    nodefluxout = if keyedarrays
        nodevariabletokeyed(nodeflux, modes, nodenames)
    else
        nodeflux
    end

    #
    Sout = if !isempty(S) && keyedarrays
        Stokeyed(S, modes, portnumbers, modes, portnumbers)
    else
        S
    end

    # the exact Jacobian of the equivalent real system at the converged
    # solution, and everything else needed to propagate a component
    # perturbation through the operating point with the implicit function
    # theorem. The Jacobian stored during the iteration is from the last
    # Newton step, so it is re-assembled here.
    operatingpoint = if returnoperatingpoint
        # retain the augmentation the solver assembled above (the promoted
        # resistor constitutive equations, the coupled inductor blocks, and
        # the gauge fixing rows, on the full augmented layout) rather than
        # reconstructing it, so the operating point cannot drift from the
        # system which was actually solved. calcresidualsensitivity masks
        # the constitutive equation rows of one promoted resistor out of it;
        # the coupled inductor and gauge entries lie outside that mask.
        # An operating point is a host object, whichever backend solved for
        # it. What differentiates it evaluates transforms of one direction at
        # a time, per component, so it needs a working system and not a
        # snapshot of arrays; a system on a backend would make every one of
        # those a transfer, spread across the callers, invisible to a test
        # suite which has no device. So a solve on a backend gets a host
        # twin, built once here from the same ingredients, and everything
        # downstream is the host code the tests already cover.
        #
        # The same reasoning gives the Jacobian: its values are on the
        # backend behind a structure held transposed, and what reads it are
        # sparse direct factorizations and a structural row mask, so it comes
        # home too. Both are the size of the state and are retained once per
        # pump solve rather than produced per signal frequency; the signal
        # sweep that follows still runs wherever it was asked to.
        opsys = if phimatrix isa Array
            setpoint!(sys, x)
            sys
        else
            hostphimatrix = zeros(eltype(phimatrix), size(phimatrix))
            hostphitd, hostirfftplan, hostrfftplan =
                plan_applynl(hostphimatrix, CPU())
            twin = HBSystem(Rbnm, invLnm, Gnm, Cnm, wmodesm, wmodes2m, bnm,
                Ljb, Ljbm, Lmean, Nbranches, freqindexmap, conjsourceindices,
                conjtargetindices, hostphimatrix, hostphitd, hostirfftplan,
                hostrfftplan, modelayout, nothing, nothing, CPU())
            setpoint!(twin, x)
            twin
        end
        # the Jacobian is the one the solver assembled, not the twin's: it is
        # the augmentation the solve actually used
        setpoint!(sys, x)
        jacobian!(Jr, sys)
        # with an explicit block the implicit function theorem applies to
        # the canonical system, so the point carries that too: the work, the
        # canonical state and the canonical Jacobian assembled there
        dcop = if dcexplicit && !isnothing(dccanonical)
            jp = canonicaljacobianplan(hostsparse(Jr), canonwork)
            canonicaljacobian!(jp, hostsparse(Jr))
            DCOperatingPoint(canonwork, dccanonical, copy(jp.J), jp)
        else
            nothing
        end
        HBOperatingPoint(opsys, copy(x), hostsparse(Jr), modelayout, Nnodal,
            Lmean, wmodes, Amna,
            mnaindices, coupledbranches, Nmodes, Nnodes, dcop)
    else
        nothing
    end

    # ground is dropped and the values keyed by node name, so that this
    # reads like `nodeflux` rather than one index out of step with it
    dcout = if isnothing(dcsol)
        nothing
    else
        v = dcsol.nodevoltage[2:end]
        keyedarrays ? nodevariabletokeyed(v, nodenames) : v
    end

    return NonlinearHB(w, frequencies, nodefluxout, Rbnmout, Ljb, Lb, Ljbm,
        Nmodes, Nbranches, nodenames, portnumbers, modes, Sout, solverinfo,
        operatingpoint, dcout)

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
            # check if the mode is in the dictionary of modes
            if haskey(modedict, mode)
                # if we find the mode and the port, set that branch in bbm
                # equal to the current scaled by the mean inductance and the
                # flux quantum.
                modeindex = modedict[mode]
                key = (nodeindices[1, portindex], nodeindices[2, portindex])
                bbm[(edge2indexdict[key]-1)*Nmodes+modeindex] += Lmean*current/phi0
            else
                throw(ArgumentError(lazy"Source mode $(mode) not found."))
            end
        else
            throw(ArgumentError(lazy"Source port $(port) not found."))
        end
    end

    return nothing
end

"""
    hbnlsolve(w, Nharmonics, sources, circuit::Circuit,
        circuitdefs = Dict{Symbol,Number}(); sorting = :name,
        keyword arguments...)

Nonlinear harmonic balance solution of a typed [`Circuit`](@ref); see
[`hbsolve`](@ref).
"""
function hbnlsolve(w::NTuple{N,Number}, Nharmonics::NTuple{N,Int}, sources,
        circuit::Circuit, circuitdefs::AbstractDict = Dict{Symbol,Number}();
        sorting::Symbol = :name, kwargs...) where N
    return hbnlsolve(w, Nharmonics, sources, elaborate(circuit),
        circuitdefs; sorting = sorting, kwargs...)
end
