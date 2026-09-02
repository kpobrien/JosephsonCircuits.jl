# =========================================================================
# The nonlinear harmonic balance solve.
# =========================================================================

"""
    hbnlsolve(w::NTuple{N,Number}, Nharmonics::NTuple{N,Int}, sources,
        circuit, circuitdefs; iterations = 1000, maxharmonics = Nharmonics,
        maxintermodorder = Inf, dc = false, odd = true, even = false,
        ftol = 1e-8, rtol = 0.0, method = :newtonkrylov,
        andersondepth = method == :quasinewton ? 5 : 0, x0 = nothing,
        symfreqvar = nothing, sorting = :number, keyedarrays = true,
        sensitivitynames = String[], returnoperatingpoint = false,
        krylovcouplingmodes = :none, krylovrecycle = 0, krylovharvest = 8,
        krylovkwargs = (;), stagedkwargs = (;), factorization = nothing,
        backend = CPU(), precision = Float64, debugJacobian = false,
        linearsolver = InternalGMRES(), returnsystem = false,
        assemblejacobian = true)

Solve the nonlinear harmonic balance problem of a circuit driven by any
number of strong tones (pumps) at any number of ports, including direct
current biases and flux pumping through a current source and a mutual
inductor. [`hblinsolve`](@ref) linearizes the circuit about the operating
point found here; [`hbsolve`](@ref) runs the two in sequence.

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
auxiliary variables have magnitudes comparable to the node fluxes. One
gauge fixing equation per floating inductive or Josephson subnetwork and
zero frequency mode makes circuits which are singular at direct current in
a purely nodal formulation (nodes or subnetworks with no inductive path to
ground) solvable without workaround inductors; if the net direct current
injected into such a subnetwork is nonzero no periodic solution exists and
an `ArgumentError` is thrown. The reported residual norms are those of the
augmented system, and the returned structure contains only the node fluxes
and the original incidence matrix. Commensurate drive frequencies whose
retained intermodulation products reach (numerically) zero frequency are
rejected with an `ArgumentError`. See `src/mna.jl`.

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
- `circuit`: a typed [`Circuit`](@ref), a legacy netlist of
    `(name, node1, node2, value)` tuples, or a [`CompiledCircuit`](@ref).
- `circuitdefs`: a dictionary from the symbols or symbolic variables used
    as component values to their numerical values. Optional when every
    component value is numeric.

# Keywords
- `iterations = 1000`: the maximum number of nonlinear solver iterations
    before it returns unconverged.
- `maxharmonics = Nharmonics`: an upper bound on the absolute harmonic
    index retained for each tone; see [`truncfreqs`](@ref).
- `maxintermodorder = Inf`: keep only the modes whose harmonic indices
    have an absolute sum of at most this order, a diamond truncation of the
    multi-tone Fourier space.
- `dc = false`: retain the zero frequency mode.
- `odd = true`: retain the odd harmonics, which four wave mixing couples
    through.
- `even = false`: retain the even harmonics, which three wave mixing
    couples through.
$(_DOC_FTOL)
$(_DOC_METHOD)
$(_DOC_STAGEDKWARGS)
$(_DOC_ANDERSON)
$(_DOC_NLKWARGS)
- `symfreqvar = nothing`: the symbolic frequency variable, such as `w`,
    when component values are expressions in the frequency.
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

# A fully numeric circuit needs no component definitions.
function hbnlsolve(w::NTuple{N,Number}, Nharmonics::NTuple{N,Int}, sources,
    circuit; kwargs...) where {N}
    return hbnlsolve(w, Nharmonics, sources, circuit, Dict{Any,Any}();
        kwargs...)
end

"""
    hbnlsolve(w::NTuple{N,Number}, sources, frequencies::Frequencies{N},
        indices::FourierIndices{N}, psc::CompiledCircuit, cg::CircuitGraph,
        nm::CircuitMatrices; kwargs...)

The nonlinear harmonic balance solve on an already compiled circuit `psc`
with its graph `cg` and matrices `nm`, at the mode set `frequencies` with
its Fourier indices `indices`. This is what the other methods call after
building those; it takes every keyword of the general method except the
ones which describe the mode set (`Nharmonics`, `maxharmonics`,
`maxintermodorder`, `dc`, `odd`, `even`, `sorting`, `stagedkwargs`).

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

See the general [`hbnlsolve`](@ref) docstring for the formulation and the
keywords.
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

    # A solver object (`NewtonKrylov()`, ...) carries its own options; a
    # symbol is the plain spelling. Resolve it first, since the setup below
    # reads `method` to decide whether to build the real representation.
    solverobject = method
    method = solvermethod(method)
    andersondepth = solverobject isa QuasiNewton ?
        solverobject.andersondepth : andersondepth
    krylovkwargs = merge(krylovkwargs, solverkwargs(solverobject))

    # the default factorization matches the backend: KLU on the host, cuDSS
    # on a device, where a host factorization could not be applied to the
    # device vectors of the Krylov iteration
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

    # Reject non-finite inputs up front: the frequency and source checks
    # below compare against bounds built from these values, and an infinite
    # or NaN value would pass or fail those comparisons silently.
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
    # The residual is computed with cyclic Fourier transforms, so its exact
    # Jacobian couples modes whose differences alias back onto the sampled
    # grid. `hbmatind` with `alias = true` includes those couplings (the
    # sum couplings of `hbconjmatind` always alias), and the exact real
    # Jacobian of `method = :newton` is assembled from them, so it is the
    # exact derivative for multi-tone problems too; it matches the
    # matrix-free products of `HBSystem` to machine precision. The complex
    # holomorphic Jacobian of `method = :quasinewton` keeps the truncated
    # (non-aliased) indices: it is an approximation either way, and the
    # aliased couplings would only make it denser to factorize.
    Amatrixindicesaliased = hbmatind(frequencies; alias = true)[2]

    # generate the frequencies of the modes
    Nmodes = length(modes)
    wmodes = calcmodefreqs(w,modes)

    # Only the all-zero mode tuple may have zero frequency: it is exactly
    # 0.0 by construction, which is what the `iszero` tests of the gauge
    # fixing and source checks rely on. Any other mode whose frequency
    # cancels to zero, or to within roundoff (commensurate drives), would
    # duplicate the direct current coordinate with different conjugacy
    # assumptions and give vanishing capacitor and resistor stamps, so it
    # is rejected here rather than allowed to form a singular system.
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

    # Nondimensionalize with the solver inductance scale Z0/w0 (see
    # `calcsolverscale`): the scaled entries are of order one for circuits
    # driven near their characteristic impedance and frequency, and `ftol`
    # is independent of the unit system. The scale multiplies rows only, so
    # the node fluxes and all physical outputs are unchanged. The local
    # name `Lmean` is what the matrix and plan interfaces call the scale.
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
    # resistors stay as node conductances, which is algebraically identical
    # and keeps the augmented system smaller. The (empty) list is the hook
    # for elements a conductance cannot express, such as voltage sources.
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
    # periodic flux carries no voltage, so a resistor is an open circuit
    # there; the missing coordinate is the average voltage, constant on each
    # static flux component because an inductor or a zero voltage junction
    # is a short at direct current. The average voltages are solved as
    # explicit unknowns with their own transport rows (see dcconductance.jl),
    # so the system is built on the applied source and nothing is
    # subtracted from it. `bnmsource` is a copy because `bnm` is handed to
    # the `HBSystem`, which writes into it when the drive is scaled, while
    # the direct current rows and the port outputs need the source as
    # applied.
    bnmsource = copy(bnm)
    dcplan = dcconductanceplan(floatingcomponents, Gnm, wmodes, Nmodes,
        Nnodes)

    # Most circuits inject no direct current, and then every average voltage
    # and block port current is zero: the explicit block would cost every
    # residual, product and preconditioner application for an answer of
    # zero. Without it, the `i = 0` zero frequency rows the scattering stamp
    # writes are exactly right.
    dcexplicit = !isnothing(dcplan) && dcinjected(dcplan, bnmsource, Nmodes)

    # With direct current active the average voltages, their transport rows
    # and the blocks' zero frequency rows live in a wrapper around the
    # `HBSystem`, not in it, so the raw system is not the one being solved.
    # Handing it out or differentiating it would define a different
    # operating point, so those requests are refused.
    if dcexplicit
    end
    # the canonical Jacobian is assembled on the host as a `SparseMatrixCSC`;
    # there is no device assembly
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
    # The scattering block contribution: the pump mode frequencies are
    # fixed, so the blocks' constitutive equations
    # im*w_m*Lmean*B(w_m)*phi - C(w_m)*i = 0 and the Kirchhoff current law
    # couplings of their auxiliary port currents form a constant matrix,
    # folded into the augmentation like the promoted resistor equations
    # (see `scatteringlinearterm`). The stamped blocks are kept because the
    # explicit direct current path needs each block's auxiliary base.
    stampedblocks = StampedScatteringBlock[]
    if Nauxscattering > 0
        # a constant matrix built once on the host; the augmented system it
        # produces is what every backend then solves
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

    # The Jacobian and its assembly plan for the chosen method. A plan maps
    # the Fourier coefficients of cos(phi(t)) and the linear term matrices
    # straight into the nonzeros of the Jacobian, so no branch matrices,
    # incidence products or real conversions are formed while iterating.
    # Only the chosen method's machinery is built, unless `debugJacobian`
    # asks for both to compare them.
    #
    # The complex Jacobian. On a backend the structure is stored transposed,
    # because a device factorization is compressed by rows, with its values
    # on the backend; the host `Jx` is built either way, because
    # `debugJacobian` compares against it and the sensitivity calculation
    # reads its structure.
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

    # `:newtonkrylov` is deliberately absent: its steps come from the
    # matrix-free Jacobian-vector product, and its preconditioner assembles
    # its own much sparser restricted plan, so the full Jacobian plan (the
    # largest object in a multi-tone solve) is never built.
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

    # The evaluation object of the nonlinear system: the residual, the
    # matrix-free Jacobian-vector and Hessian-vector products and the
    # assembled Jacobians are all computed through it, in the complex or the
    # equivalent real representation. The real form of the linear term is
    # needed exactly when the solve uses the real representation: every
    # method but `:quasinewton`, including an external solver and the
    # `returnsystem` path, since the residual is not complex differentiable.
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
        # write back the real representation of the point through the plan;
        # `complex_to_real!` walks a BitVector serially and is host only
        applycomplextoreal!(xr, sys.nonlineartermplan, sys.x)
        return nothing
    end

    # the diagnostics of each solver invocation, returned in the output
    solverstages = IterationInfo[]

    # reserved for source stepping continuation; always NaN
    sourcefold = Ref(NaN)

    # The canonical layout of the state with direct current active, built
    # once for whichever method uses it: the matrix-free and the direct
    # paths solve the same system in the same coordinates and differ only
    # in how they apply the Jacobian.
    canonwork = if dcexplicit
        tr = dcexplicit ? transportrows(dcplan, bnmsource, Nmodes) : nothing
        Lc = tobackend(backend, compositelayout(modelayout, frequencies.modes;
            nvdc = isnothing(tr) ? 0 : nvoltages(tr)))
        # a block's zero frequency row is `i = 0` in the stamp; with an
        # average voltage to respond to, it becomes the block's own relation
        br = (dcexplicit && Nauxscattering > 0) ?
            dcblockrows(stampedblocks, dcplan.componentof, Nmodes,
                dcplan.modeindex, Lc.nac, Nnodes - 1, Lmean) : nothing
        # the zero frequency block holds one entry per node and then one per
        # auxiliary unknown; only the nodal ones carry the coupling. (Not
        # named `w`, which is this function's drive frequency argument.)
        CanonicalWork(Lc, tobackend(backend,
                convert(Vector{precision}, xr));
            transport = tr, blockrows = br, nnodaldc = Nnodes - 1)
    else
        nothing
    end

    if returnsystem
        # Everything an external solver needs, without solving: the
        # evaluation object, the initial value, the real representation
        # layout, and the assembled real Jacobian if it was asked for.
        return (sys=sys, xr=xr, Fr=Fr, modelayout=modelayout,
            Jr=(assemblejacobian ? Jr : nothing), Nnodal=Nnodal,
            dcplan=dcplan, dcsol=dcsol, bnmsource=bnmsource,
            # with direct current active `sys` alone is not the system: the
            # average voltages and the blocks' zero frequency rows are here
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
    # sum. The scale (see `calcsolverscale`) comes from the port impedances,
    # and a circuit whose interior sits at a very different impedance is
    # left with a scaled source many orders above one and a rounding floor
    # above the default `ftol`. Asking for less than the floor is asking the
    # iteration to converge on noise, so the tolerance is raised to the
    # floor when the floor is larger. The floor is a `sqrt(n)*eps`
    # accumulation over the terms with room to spare; for a circuit driven
    # near its characteristic impedance it is far below `ftol`.
    ftol = max(ftol,
        16*sqrt(length(xr))*eps(real(eltype(xr)))*norm(bnmsource))

    # Solve the nonlinear system. The canonical layout is supported by the
    # matrix-free path; the assembled Jacobian methods would additionally
    # need the permuted Jacobian `P J P'`, which is not implemented, so the
    # request is refused rather than ignored.
    info = if method == :quasinewton

        solveonbackend!(fj!, F, Jxb, x, backend; iterations = iterations,
            ftol = ftol, rtol = rtol, andersondepth = andersondepth,
            factorization = factorization)

    elseif method == :newton

        # solve the equivalent real system with the exact real Jacobian,
        # then convert back to complex
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

        # The right preconditioner is the Jacobian with its mode coupling
        # restricted (see `ModeCouplingPreconditioner`), so its factorization
        # is a fraction of the full Jacobian's; the residual, the
        # Jacobian-vector product and the solution are always those of the
        # requested truncation.
        #
        # The default `:none` keeps the mode block diagonal, whose
        # factorization is a batch of small independent per mode solves. On
        # a strongly pumped line it stalls on its own; escalation, which
        # grows the coupling set on repeated linear failures, rescues it and
        # in practice fires once or twice. `:band => p` restricts the
        # coupling by harmonic offset, which is the restriction the Toeplitz
        # structure of the nonlinear term suggests; its fill grows linearly
        # rather than quadratically in the mode count, so it wins over the
        # full factorization once there are enough modes, and it escalates
        # one offset at a time. `krylovrecycle > 0` additionally wraps the
        # preconditioner in a recycled deflation subspace, which also rescues
        # the block diagonal and needs no sparse factorization at all, but
        # is off by default.
        base = ModeCouplingPreconditioner(sys, Amatrixindicesaliased,
            Amatrixconjindices, Ljb, Lmean, Rbnm, Nmodes, Nbranches, Nfreq,
            invLnm, Gnm, Cnm, modelayout;
            couplingmodes = krylovcouplingmodes,
            factorization = factorization, precision = precision,
            Amatrixmodes = Amatrixmodes)

        # the Krylov workspaces are allocated `similar` to the vectors handed
        # in, so device vectors put the whole iteration on the device; on
        # `CPU()` these are no-ops
        xrb = tobackend(backend, convert(Vector{precision}, xr))
        Frb = tobackend(backend, convert(Vector{precision}, Fr))

        # built from the vector rather than its length, so the deflation
        # subspace lives where the iteration does
        pc = if krylovrecycle > 0
            RecyclingPreconditioner(base, jvpreal!, xrb;
                kmax = krylovrecycle, kharvest = krylovharvest)
        else
            base
        end

        # The solver defaults (a long restart cycle and an eager preconditioner
        # refresh) are those of `nlsolvekrylov!` itself; `krylovkwargs`
        # overrides any of them. The refresh is eager because rebuilding the
        # block diagonal is cheap while a preconditioner frozen at zero flux,
        # where the Jacobian has no harmonic coupling, is stale at once. The
        # restart length must be long enough for the Krylov space to resolve
        # the directions a restricted preconditioner leaves, or a restart
        # discards the progress on them and the iteration count stops
        # responding to the preconditioner at all.
        info = if dcexplicit
            # The canonical layout groups the state by role rather than by
            # node. It is a permutation, so the iteration is the same one in
            # a rotated basis and must reach the same point.
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
            # report the voltages this path solved for; they agree with the
            # eliminated ones while the direct current devices are linear
            # conductances, which the tests assert
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
        # caller; the conversion walks a BitVector serially
        copyto!(xr, tohost(xrb))
        copyto!(Fr, tohost(Frb))
        real_to_complex!(x,xr,modelayout.isreal)
        real_to_complex!(F,Fr,modelayout.isreal)
        info
    elseif method == :external
        # hand the caller's root finder the system as a problem object. With
        # direct current the problem carries the augmentation too, so the
        # caller solves the posed system in the canonical unknowns.
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

    # Validate the original, ungauged Kirchhoff current law equations by
    # reconstructing their residuals from the augmented residual and the
    # state (a gauge row g adds x[g] to the augmented residual). The
    # acceptance policy is in `mnavalidatekcl`.
    if converged && !isempty(gaugeindices)
        setpoint!(sys, x)
        # `F` and `x` are host vectors while `sys` holds its state on the
        # backend, so the residual is evaluated into a backend buffer and
        # copied back for the host side validation
        Fgauge = similar(sys.x)
        residual!(Fgauge, sys)
        copyto!(F, Fgauge)
        # On the explicit path the resistor current arrives through the
        # transport coupling outside the system. Adding it back and correcting
        # the source the same way gives this check the pair the elimination
        # would have produced, so it validates the physics rather than the
        # path.
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
    # frequency, i.e. in its zero voltage state. The sine of the branch flux
    # is cached at the converged point, so the direct current each junction
    # carries costs nothing to check.
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

    # drop the auxiliary variables; the output holds only the node fluxes
    nodeflux = x[1:Nnodal]

    # the residual norms at the initial value and at the solution, for the
    # diagnostics
    normF0 = info.normresidual[1]
    normFfinal = info.normresidual[end]


    # the diagnostics of the solution process
    solverinfo = SolverInfo(solverstages, normF0, normFfinal, converged,
        sourcefold[])

    # the scattering parameters at the pump modes
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

    # The exact Jacobian of the equivalent real system at the converged
    # solution, with everything else needed to propagate a component
    # perturbation through the operating point by the implicit function
    # theorem. The Jacobian held during the iteration is from the last
    # Newton step, so it is reassembled here.
    operatingpoint = if returnoperatingpoint
        # The augmentation the solver assembled above (the coupled inductor
        # blocks and the gauge fixing rows) is kept rather than rebuilt, so
        # the operating point cannot drift from the system which was solved.
        #
        # An operating point is a host object whichever backend solved for
        # it: what differentiates it evaluates transforms one direction at a
        # time, per component, and doing that across a device transfer would
        # spread the cost over every caller. So a solve on a backend gets a
        # host twin of the system, built once here from the same ingredients,
        # and the Jacobian, read by sparse direct factorizations and a row
        # mask, comes back to the host as well. Both are retained once per
        # pump solve; the signal sweep still runs wherever it was asked to.
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
        # the Jacobian is the one the solver assembled, not the twin's
        setpoint!(sys, x)
        jacobian!(Jr, sys)
        # with an explicit direct current block the implicit function theorem
        # applies to the canonical system, so the point carries its work, its
        # state and its Jacobian too
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

    # ground is dropped and the values are keyed by node name, like
    # `nodeflux`
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

The source vector in the branch basis: for each source, the current at
its port and mode scaled by `Lmean/phi0`, at the index of that port's
branch and mode; zero elsewhere. See also [`addsources!`](@ref).

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

Fill `bbm` with the source vector of [`calcsources`](@ref). A source at a
port number or mode not in the circuit is ignored.
"""
function addsources!(bbm, modes, sources, portindices, portnumbers,
    nodeindices, edge2indexdict, Lmean, Nnodes, Nbranches, Nmodes)

    # zero the vector
    fill!(bbm,0)

    # port number -> position, and mode -> position
    portdict = Dict{eltype(portnumbers), eltype(portindices)}()
    for i in eachindex(portindices)
        portdict[portnumbers[i]] = portindices[i]
    end

    modedict = Dict{eltype(modes), Int}()
    for i in eachindex(modes)
        modedict[modes[i]] = i
    end

    for source in sources

        port = source[:port]
        mode = source[:mode]
        current = source[:current]


        if haskey(portdict, port)
            portindex = portdict[port]

            if haskey(modedict, mode)
                # the current at this port's branch and this mode, scaled by
                # the solver scale over the flux quantum
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
        circuitdefs = Dict{Symbol,Number}(); sorting = :name, kwargs...)

The nonlinear solve of a typed [`Circuit`](@ref), with every keyword of
the general method. `circuitdefs` is needed only when component values are
symbolic, and `sorting` defaults to `:name` because hierarchical net names
are not integers.
"""
function hbnlsolve(w::NTuple{N,Number}, Nharmonics::NTuple{N,Int}, sources,
        circuit::Circuit, circuitdefs::AbstractDict = Dict{Symbol,Number}();
        sorting::Symbol = :name, kwargs...) where N
    return hbnlsolve(w, Nharmonics, sources, elaborate(circuit),
        circuitdefs; sorting = sorting, kwargs...)
end
