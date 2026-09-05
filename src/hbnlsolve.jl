# the keyword pairs `kw` (a named tuple or the `Base.Pairs` a keyword
# splat carries) as a named tuple without the given keys

"""
    HBReuse()

What one solve builds and a later solve of the same circuit at new
component values takes over: the aliased mode coupling index, the padded
linear term with its augmentation ([`PaddedLinearTerm`](@ref)), the
[`HBSystem`](@ref) with its transforms and workspaces, the mode coupling
preconditioner with its structure and symbolic factorization, the
recycled deflation candidates of the last converged solve under a
[`Recycling`](@ref) or [`Floquet`](@ref) preconditioner, and the Krylov
vectors. Hand one to [`hbnlsolve`](@ref) as `reuse`; it is filled by the
first solve and rebound to the new values by every later one, so a sweep
pays for its transforms, index maps and factorization symbolics once, and
each solve starts from the deflation subspace the previous one harvested
rather than from nothing. [`hbcache`](@ref) carries one for its sweeps.

Only a [`NewtonKrylov`](@ref) method reads it. Every solve which shares it must be
of the same circuit topology at the same mode grid with the same options,
which is what a cache guarantees; a system whose sparse structure moved is
refused rather than rebound.
"""
mutable struct HBReuse
    indicesaliased::Any
    linear::Any        # the `PaddedLinearTerm`, refilled at each point
    sys::Any
    maps::Any          # the system's `ValueMaps`, found on first rebind
    preconditioner::Any
    # the `RecyclingState` or `FloquetState` of the last *converged* solve,
    # whose candidates the next solve inherits; `nothing` without recycling
    recycling::Any
    krylov::Base.RefValue{Any}
    krylovcanonical::Base.RefValue{Any}
end

HBReuse() = HBReuse(nothing, nothing, nothing, nothing, nothing, nothing,
    Ref{Any}(nothing), Ref{Any}(nothing))

# =========================================================================
# The nonlinear harmonic balance solve.
# =========================================================================

"""
    hbnlsolve(w::NTuple{N,Number}, Nharmonics::NTuple{N,Int}, sources,
        circuit, circuitdefs; iterations = 1000,
        Nevaluationharmonics = map(i -> 2i, Nharmonics),
        maxintermodorder = Inf, dc = false, odd = true, even = false,
        ftol = 1e-8, rtol = 0.0, method = NewtonKrylov(), x0 = nothing,
        symfreqvar = nothing, sorting = :number, keyedarrays = true,
        sensitivitynames = String[], returnoperatingpoint = false,
        frequencywindow = (0, Inf), backend = CPU(), debugJacobian = false,
        returnsystem = false, assemblejacobian = true)

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
- `Nharmonics::NTuple{N,Int}`: the largest absolute harmonic index of each
    tone retained as an unknown, so the modes of the returned solution. The
    length of the tuple must equal the number of non-commensurate tones. The
    nonlinearity is evaluated on the larger `Nevaluationharmonics` grid.
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
- `Nevaluationharmonics = map(i -> 2i, Nharmonics)`: the harmonics
    of each tone on the grid where the nonlinearity is sampled, at least
    `Nharmonics`. The nonlinear products of the retained modes reach
    higher orders, and a grid with `M` harmonics folds order `p` back to
    `p - (2M + 1)` (aliasing). Products of two retained modes reach twice
    `Nharmonics` and fold outside the retained set once the grid is half
    again as large, the three halves rule, but the leading nonlinearity of
    a junction is cubic, and products of three retained modes reach three
    times `Nharmonics`, which the three halves grid folds onto the tone
    itself. Twice the retained set is the default: it dealiases the cubic
    products and leaves the fifth order ones, whose folded contributions
    land outside the retained set. Only the transforms grow, the unknowns
    are the same. Measured against a grid four times the retained set on a
    64-junction line with three tones retaining (6,4,4): the unpadded grid
    is off by 1.5e-4 of the strongest mode in the modes of order four and
    1.5e-6 in the tones, three halves by 1.8e-7 and 4e-11, twice by 2.7e-11
    and 4.5e-15, at the same solve time and GMRES step count, which the
    padding took from 196 to 16.
- `maxharmonics`: deprecated and ignored with a warning; `Nharmonics` is
    the retained set and `Nevaluationharmonics` the sampling grid.
- `frequencywindow = (0, Inf)`: a lower and upper bound, in the units of
    `w`, on the absolute frequency `abs(dot(w, mode))` of the retained
    modes, a truncation by frequency beside the truncations by order; the
    zero frequency mode follows `dc`. See [`truncfreqs`](@ref) for why a
    floor matters on incommensurate tones.
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
    Nevaluationharmonics::NTuple{N,Int} = map(i -> 2i, Nharmonics),
    maxharmonics = nothing,
    frequencywindow = (0, Inf),
    maxintermodorder = Inf, dc::Bool = false, odd::Bool = true,
    even::Bool = false, x0 = nothing, ftol = 1e-8,
    switchofflinesearchtol = nothing, alphamin = nothing,
    method::AbstractHBNonlinearSolver = NewtonKrylov(),
    symfreqvar = nothing, sorting = :number, keyedarrays::Bool = true,
    sensitivitynames::Vector{String} = String[],
    returnoperatingpoint::Bool = false,
    backend = CPU(), debugJacobian = false,
    returnsystem::Bool = false, assemblejacobian::Bool = true,
    ) where {N}

    # deprecation warning for maxharmonics, whose role `Nharmonics` took
    # when the sampling grid became `Nevaluationharmonics`.
    if !isnothing(maxharmonics)
        Base.depwarn(lazy"The `maxharmonics` kwarg is deprecated and no longer used. `Nharmonics` is the retained set of modes and `Nevaluationharmonics` the grid on which the nonlinearity is sampled. Please remove it to avoid errors in future versions.", :hbnlsolve; force=true)
    end

    all(map(>=, Nevaluationharmonics, Nharmonics)) || throw(ArgumentError(
        lazy"`Nevaluationharmonics` = $(Nevaluationharmonics) must be at least `Nharmonics` = $(Nharmonics) in every tone."))

    if method isa Staged
        return stagedsolve(method, w, Nharmonics, sources, circuit, circuitdefs;
            iterations = iterations,
            Nevaluationharmonics = Nevaluationharmonics,
            frequencywindow = frequencywindow,
            maxintermodorder = maxintermodorder, dc = dc, odd = odd,
            even = even, ftol = ftol, symfreqvar = symfreqvar,
            sorting = sorting, keyedarrays = keyedarrays,
            sensitivitynames = sensitivitynames,
            returnoperatingpoint = returnoperatingpoint, backend = backend)
    end

    # calculate the frequency struct
    freq = removeconjfreqs(
        truncfreqs(
            calcfreqsrdft(Nevaluationharmonics); dc = dc, odd = odd,
            even = even, maxintermodorder = maxintermodorder,
            maxharmonics = Nharmonics, w = w,
            frequencywindow = frequencywindow,
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
        method = method,
        symfreqvar = symfreqvar, keyedarrays = keyedarrays,
        sensitivitynames = sensitivitynames,
        returnoperatingpoint = returnoperatingpoint,
        backend = backend, debugJacobian = debugJacobian,
        returnsystem = returnsystem,
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
ones which describe the mode set (`Nharmonics`, `Nevaluationharmonics`,
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
    method::AbstractHBNonlinearSolver = NewtonKrylov(),
    symfreqvar = nothing, keyedarrays::Bool = true,
    sensitivitynames::Vector{String} = String[],
    returnoperatingpoint::Bool = false,
    backend = CPU(), debugJacobian = false,
    returnsystem::Bool = false, assemblejacobian::Bool = true,
    reuse::Union{Nothing,HBReuse} = nothing,
    )

    method isa Staged && throw(ArgumentError(
        "a `Staged` method is solved by the general `hbnlsolve` method, which builds each stage's system."))
    precision = solverprecision(method)
    # the sparse factorization of the direct methods
    directfactorization = something(
        method isa Newton || method isa QuasiNewton ? method.factorization : nothing,
        KLUfactorization())

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
    # what a previous solve of this circuit built and this one takes over;
    # see `HBReuse`. Only the matrix free path is built to be rebound.
    reusing = !isnothing(reuse) && method isa NewtonKrylov
    Amatrixindicesaliased = if reusing && !isnothing(reuse.indicesaliased)
        reuse.indicesaliased
    else
        a = hbmatind(frequencies; alias = true)[2]
        reusing && (reuse.indicesaliased = a)
        a
    end

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
    # a reused system already holds this array, its time domain twin and
    # the plans, and they are read from it below rather than allocated
    reusesys = reusing && !isnothing(reuse.sys)
    if reusesys
        rs = reuse.sys
        (size(rs.phimatrix) == Nwtuple &&
         eltype(rs.phimatrix) == Complex{precision} &&
         KernelAbstractions.get_backend(rs.phimatrix) == backend) ||
            throw(ArgumentError("the system handed in for reuse was built for a different mode grid, precision or backend; hand in a fresh `HBReuse`."))
    end
    phimatrix, phimatrixtd, irfftplan, rfftplan = if reusesys
        reuse.sys.phimatrix, nothing, nothing, nothing
    else
        pm = tobackend(backend, zeros(Complex{precision}, Nwtuple))
        # create an array to hold the time domain data for the RFFT. also
        # generate the plans, which come from the package extension on a
        # device backend.
        (pm, plan_applynl(pm, backend)...)
    end

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
    # periodic flux carries no voltage, so a resistor is an open circuit at
    # DC. The missing coordinate is the average voltage, constant on each
    # static flux component because an inductor or zero-voltage junction is
    # a short there. The average voltages are unknowns with their own
    # transport rows, and the resistor and block currents they drive reach
    # the nodes through the coupling, so the system is built on the applied
    # source and nothing is subtracted from it. `bnmsource` is a copy and
    # not an alias: `bnm` is handed to the `HBSystem`, which owns it from
    # then on and writes into it when the drive is scaled, while the
    # transport rows read the drive as applied and the continuation
    # interface scales their constant term itself.
    bnmsource = copy(bnm)
    dcplan = dcconductanceplan(floatingcomponents, Gnm, wmodes, Nmodes,
        Nnodes)

    # The direct current block exists when there is something in it: an
    # average voltage to find, or a scattering block whose zero frequency
    # current the stamp's `i = 0` row may or may not have right. It is
    # classified in either case, below, once the system exists to classify
    # it against; it is carried through the solve only when it has
    # something to find. Nothing injects direct current in most circuits,
    # and then every average voltage and every determined block current is
    # zero: the explicit block would carry a subsystem whose answer is the
    # zero it starts at and charge for it on every residual, product and
    # preconditioner application. The periodic system as stamped is then
    # exact, and is solved alone.
    hasdc = !isnothing(dcplan) &&
        (!isempty(dcplan.components) || Nauxscattering > 0)
    dcexplicit = hasdc && dcinjected(dcplan, bnmsource, Nmodes)

    # The average voltages, their transport rows and the blocks' zero
    # frequency rows live in a wrapper around the system rather than in it,
    # so the raw `HBSystem` is not the system being solved: its scattering
    # rows still say `i = 0` and its resistors still carry no direct
    # current. Handing that object out, or differentiating it, would define
    # a second and different operating point. Refuse instead, until the
    # composite system exists to hand out.
    # the canonical Jacobian is assembled as a host `SparseMatrixCSC`; on a
    # device backend `Jr` is not one, and there is no device assembly yet
    if dcexplicit && method isa Newton && !(backend isa CPU)
        throw(ArgumentError("a circuit which injects direct current is solved with an assembled canonical Jacobian under `Newton()`, and that assembly is host only. Use `NewtonKrylov()` on this backend, which is matrix free, or run `Newton()` on the CPU."))
    end
    if dcexplicit && method isa QuasiNewton
        throw(ArgumentError("a circuit which injects direct current is solved with the average voltages as unknowns, which needs the real system; `QuasiNewton()` solves the complex holomorphic one. Use `NewtonKrylov()`, `Newton()` or an `ExternalSolver`."))
    end
    # the zero every average voltage sits at when none is injected, which
    # is the answer for a circuit with a zero frequency mode and no direct
    # current, and is replaced by the solved voltages when there is some
    dcsol = isnothing(dcplan) ? nothing :
        dcsolutionfrom(dcplan, zeros(Float64, length(dcplan.components)))
    # the converged canonical state, which is the point the whole system is
    # differentiated at when the direct current block is active
    dccanonical = nothing

    gaugeindices = calcdcgaugeindices(floatingcomponents, wmodes, Nmodes)
    # The padding, the augmentation and their union are structure; a sweep
    # assembles them once and refills the values (see PaddedLinearTerm).
    # The coupled inductor rows are the augmentation's only value
    # dependent entries and are reassembled, small, at every point.
    linear = reusing ? reuse.linear : nothing
    AmnaL = calcAmnaind(coupledbranches, Lb, Mb, cg.Rbn, Nmodes,
        Nnodal + Nauxr, Nnodal + Naux, Lmean)
    if isnothing(linear)
        Amna = calcAmna(mnaindices, nodeindices, vvn, gaugeindices, wmodes,
            Nmodes, Nnodes, Lmean)
        # pad with the coupled inductor auxiliary variables and add their
        # constitutive equations and Kirchhoff current law couplings, which are
        # real and frequency independent (see calcAmnaind).
        Amna = mnapad(Amna, length(coupledbranches)*Nmodes + Nauxscattering)
        Amna = spaddkeepzeros(Amna, AmnaL)
        # The scattering block contribution: the pump mode frequencies are
        # fixed, so the blocks' constitutive equations
        # im*w_m*Lmean*B(w_m)*phi - C(w_m)*i = 0 and the Kirchhoff current law
        # couplings of their auxiliary port currents form a constant matrix,
        # folded into the augmentation like the promoted resistor equations
        # (see `scatteringlinearterm`). The stamped blocks are kept because
        # the explicit direct current path needs each block's auxiliary base.
        stampedblocks = StampedScatteringBlock[]
        if Nauxscattering > 0
            # a constant matrix built once on the host; the augmented system
            # it produces is what every backend then solves
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
        invLnm0 = invLnm
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
        # promoted resistors would put values into calcAmna's rows, which
        # the refill does not follow; there are none today (see mnaindices)
        if reusing && isempty(mnaindices)
            reuse.linear = PaddedLinearTerm(Rbnm, invLnm, Gnm, Cnm, Amna,
                AmnaL, invLnm0, bnm, wmodesm, wmodes2m, stampedblocks)
        end
    else
        refill!(linear, invLnm, Gnm, Cnm,
            isempty(coupledbranches) ? nothing : AmnaL, bbm, Rbnm)
        Rbnm = linear.Rbnm
        invLnm = linear.invLnm
        Gnm = linear.Gnm
        Cnm = linear.Cnm
        Amna = linear.Amna
        bnm = linear.bnm
        wmodesm = linear.wmodesm
        wmodes2m = linear.wmodes2m
        stampedblocks = linear.stampedblocks
    end
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
    Jx, complexjacobianplan = if method isa QuasiNewton || debugJacobian
        plancomplexjacobian(Amatrixindices, Ljb, Lmean, Rbnm, Nmodes,
            Nbranches, Nfreq, invLnm, Gnm, Cnm)
    else
        nothing, nothing
    end
    devicex = method isa QuasiNewton && !(backend isa CPU)
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
    Jr, realjacobianplan = if method isa Newton || debugJacobian ||
            ((returnsystem || method isa ExternalSolver) && assemblejacobian) ||
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
    realrepresentation = method isa Newton || method isa NewtonKrylov ||
        method isa ExternalSolver || debugJacobian || returnoperatingpoint ||
        returnsystem
    sys = if reusesys
        # the same transforms, maps, kernels and workspaces at the new
        # values, written through maps found once so nothing is allocated
        # a conversion with no fixed map is remembered as such, so it is
        # probed once and not at every point
        isnothing(reuse.maps) &&
            (reuse.maps = something(valuemaps(reuse.sys), :none))
        s = rebind!(reuse.sys, invLnm, Gnm, Cnm, bnm, Ljb, Ljbm, Lmean;
            maps = reuse.maps === :none ? nothing : reuse.maps)
        reuse.sys = s
        s
    else
        s = HBSystem(Rbnm, invLnm, Gnm, Cnm, wmodesm, wmodes2m, bnm,
            Ljb, Ljbm, Lmean, Nbranches, freqindexmap, conjsourceindices,
            conjtargetindices, phimatrix, phimatrixtd, irfftplan, rfftplan,
            modelayout, realjacobianplan, complexjacobianplan, backend;
            realbackward = realrepresentation)
        reusing && (reuse.sys = s)
        s
    end

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
    # `xr` is read and not rewritten: the real representation has one slot
    # per real entry and no redundant ones, so the round trip through the
    # complex form is a copy and writing it back was a pass over the state,
    # and on a device a crossing of the bus, per residual for nothing
    function fjreal!(Fr, Jr, xr)
        setpoint!(sys, xr)
        isnothing(Fr) || residual!(Fr, sys)
        isnothing(Jr) || jacobian!(Jr, sys)
        return nothing
    end

    # the diagnostics of each solver invocation, returned in the output
    solverstages = IterationInfo[]

    # reserved for source stepping continuation; always NaN
    sourcefold = Ref(NaN)

    # use this for debugging purposes to return the residual and Jacobian
    # functions along with the ingredients from which the Jacobians are
    # assembled, so reference implementations can be constructed, eg. in the
    # tests.
    # The canonical layout, built once for whichever method uses it. Both
    # the matrix free path and the direct one solve the same system in the
    # same coordinates; only the way they apply the Jacobian differs.
    #
    # Built whenever the block exists, because building it is what
    # classifies it: the complete descriptor -- every transport row and
    # every block relation -- is what says whether the direct current is
    # determined, and a short or an ideal through in parallel with an
    # inductive path is refused here whether or not anything drives it.
    # When nothing does, the classification is all that is wanted of it and
    # the periodic system is solved as stamped, so the work is not kept.
    canonwork = if hasdc
        tr = transportrows(dcplan, bnmsource, Nmodes)
        Lc = tobackend(backend, compositelayout(modelayout, frequencies.modes;
            nvdc = nvoltages(tr)))
        # a block's zero frequency row is `i = 0` in the stamp; with an
        # average voltage to respond to it becomes the block's own relation,
        # which is what makes it visible at direct current. A block with no
        # zero frequency data has to give one only when it is asked to
        # carry direct current, and stays open otherwise.
        br = Nauxscattering > 0 ?
            dcblockrows(stampedblocks, dcplan.componentof, Nmodes,
                dcplan.modeindex, Nnodes - 1, Lmean;
                required = dcexplicit) : nothing
        # the zero frequency block holds one entry per node and then one per
        # auxiliary unknown; only the nodal ones carry the coupling.
        #
        # Not named `w`: that is this function's drive frequency, and
        # assigning to it here rebound the argument, so every circuit which
        # injected direct current reported the work object as its drive and
        # `hbsolve` could not compute the mode frequencies from it.
        work = CanonicalWork(Lc, tobackend(backend,
                convert(Vector{precision}, xr));
            transport = tr, blockrows = br, nnodaldc = Nnodes - 1)
        dcexplicit ? work : nothing
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
    info = if method isa QuasiNewton

        solveonbackend!(fj!, F, Jxb, x, backend; iterations = iterations,
            ftol = ftol, rtol = rtol, andersondepth = method.anderson,
            factorization = directfactorization)

    elseif method isa Newton

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
                andersondepth = 0, factorization = directfactorization)
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
                andersondepth = 0, factorization = directfactorization)
        end
        real_to_complex!(x,xr,modelayout.isreal)
        real_to_complex!(F,Fr,modelayout.isreal)
        info
    elseif method isa NewtonKrylov

        # the matrix free real Jacobian
        jvpreal!(Jvr, vr) = jacobianvectorproduct!(Jvr, sys, vr)

        # The right preconditioner is the Jacobian with its mode coupling
        # restricted (see `ModeCouplingPreconditioner`), so its factorization
        # is a fraction of the full Jacobian's; the residual, the
        # Jacobian-vector product and the solution are always those of the
        # requested truncation.
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
        # the preconditioner the method asks for: a mode coupling
        # preconditioner, possibly wrapped in a deflation below
        spec = method.preconditioner
        inner = spec isa AbstractModeCoupling ? spec : spec.inner
        base = if reusing && !isnothing(reuse.preconditioner)
            # the same structure, coupling set and symbolic factorization,
            # with the linear term and the junction coefficients refreshed
            rebind!(reuse.preconditioner, sys)
        else
            b = ModeCouplingPreconditioner(sys, Amatrixindicesaliased,
                Amatrixconjindices, Ljb, Lmean, Rbnm, Nmodes, Nbranches,
                Nfreq, invLnm, Gnm, Cnm, modelayout; spec = inner,
                precision = precision, Amatrixmodes = Amatrixmodes)
            reusing && (reuse.preconditioner = b)
            b
        end

        # the Krylov workspaces are allocated `similar` to the vectors handed
        # in, so device vectors put the whole iteration on the device; on
        # `CPU()` these are no-ops
        xrb = tobackend(backend, convert(Vector{precision}, xr))
        Frb = tobackend(backend, convert(Vector{precision}, Fr))

        # With direct current the canonical layout groups the state by role
        # rather than by node. It is a permutation, so the iteration is the
        # same one in a rotated basis and must reach the same point. The
        # recycling wrapper, when there is one, goes *outside* the canonical
        # wrapper: its candidates are then corrections of the whole canonical
        # state, and the Arnoldi factorization it harvests is that of the
        # operator GMRES actually ran, rather than a restriction of it.
        if dcexplicit
            work = canonwork
            L = work.layout
            ucb = similar(xrb, canonicaldim(L))
            Fcb = similar(Frb, canonicaldim(L))
            # the explicit voltages start at zero, which is the answer when
            # nothing injects direct current; `gathercanonical!` writes the
            # flux blocks only and would leave them undefined
            fill!(ucb, zero(eltype(ucb)))
            gathercanonical!(ucb, xrb, L)
            pcbase = CanonicalPreconditioner(base, work)
            jvsolve = canonicaljvp(jvpreal!, work)
            usolve, Fsolve = ucb, Fcb
        else
            pcbase = base
            jvsolve = jvpreal!
            usolve, Fsolve = xrb, Frb
        end

        # Built from the solve vector rather than its length, so the
        # deflation lives where the iteration does. Across a sweep the
        # wrapper is rebuilt, since its product closes over the system, but
        # the candidates are inherited: the directions the block diagonal
        # leaves are those of the circuit, which a nearby operating point
        # shares, and the pair is rebuilt against the rebound base at the
        # first refresh of the new solve. The inherited state is copied and
        # committed back only after a converged solve (below), so a failed
        # point cannot seed the next one.
        pc = if spec isa Floquet
            prev = reusing ? reuse.recycling : nothing
            # the residual-image form with physical candidates: `kharvest`
            # is the number of singular directions per harvest, the
            # harmonic Ritz ones coming on top of it
            state = if prev isa FloquetState && size(prev.X, 1) == length(usolve)
                copy(prev)
            else
                FloquetState(usolve)
            end
            FloquetPreconditioner(pcbase, jvsolve, usolve;
                kmax = spec.size, kharvest = spec.harvest, nritz = spec.ritz,
                kcandidate = spec.candidates, benefittol = spec.benefittol,
                cycleharvest = spec.cycleharvest, state = state,
                (isnothing(spec.ranktol) ? (;) : (; ranktol = spec.ranktol))...)
        elseif spec isa Recycling
            prev = reusing ? reuse.recycling : nothing
            state = if prev isa RecyclingState && size(prev.U, 1) == length(usolve)
                copy(prev)
            else
                RecyclingState(usolve)
            end
            RecyclingPreconditioner(pcbase, jvsolve, usolve;
                kmax = spec.size, kharvest = spec.harvest, form = spec.form,
                state = state)
        else
            pcbase
        end

        # the linear solver, the refresh policy and the escalation are the
        # method's; see `NewtonKrylov`
        info = if dcexplicit
            out = nlsolvekrylov!(canonicalresidual(fjreal!, canonwork),
                jvsolve, Fsolve, usolve, pc;
                iterations = iterations, ftol = ftol, rtol = rtol,
                linearsolver = method.linearsolver, refresh = method.refresh,
                escalate = method.escalate,
                workspace = reusing ? reuse.krylovcanonical : nothing)
            L = canonwork.layout
            scattercanonical!(xrb, usolve, L)
            scattercanonical!(Frb, Fsolve, L)
            # the voltages are read off the canonical state the solve
            # converged in; the block currents beside them stay in the
            # state, where the operating point and the sensitivities read
            # them
            dccanonical = Array(usolve)
            dcsol = dcsolutionfrom(dcplan,
                Array(view(usolve, voltagerange(L))))
            out
        else
            nlsolvekrylov!(fjreal!, jvsolve, Fsolve, usolve, pc;
                iterations = iterations, ftol = ftol, rtol = rtol,
                linearsolver = method.linearsolver, refresh = method.refresh,
                escalate = method.escalate,
                workspace = reusing ? reuse.krylov : nothing)
        end
        if reusing && !(spec isa AbstractModeCoupling) && info.converged
            reuse.recycling = pc.state
        end
        # back to the host for the complex representation returned to the
        # caller; the conversion walks a BitVector serially
        copyto!(xr, tohost(xrb))
        copyto!(Fr, tohost(Frb))
        real_to_complex!(x,xr,modelayout.isreal)
        real_to_complex!(F,Fr,modelayout.isreal)
        info
    elseif method isa ExternalSolver
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
        u, extconverged = method.f(prob, copy(u0))
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
            [norm(Fc)], Float64[], Int[], Bool[], [],
            extconverged ? :converged : :external)
    else
        throw(ArgumentError("Method $(method) is not defined."))
    end


    push!(solverstages, IterationInfo(info.label, 1.0,
        info.regularization, info.converged, info.iterations,
        info.normresidual, info.alpha, info.backtracks,
        info.andersonaccepted, info.krylov, info.reason))

    if !info.converged
        @warn lazy"Solver did not converge: $(stallmessage(info.reason))."
    end

    converged = info.converged

    # validate the original, ungauged Kirchhoff current law equations
    # directly by reconstructing their residuals from the augmented
    # residual and the state (the gauge fixing equations add x[g] to the
    # augmented residual of each gauge row g). the acceptance policy,
    # including its block-relative infinity-norm tolerance and its
    # rejection of non-finite residuals, is implemented and unit tested in
    # mnavalidatekcl.
    # whether the system has been set to the accepted point below, so the
    # junction check does not set it twice
    pointset = false
    if converged && !isempty(gaugeindices)
        setpoint!(sys, x)
        pointset = true
        # `sys` holds its state wherever the backend put it, while `F` and `x`
        # are host vectors, so the residual is evaluated into a buffer like the
        # system's own state and copied back for the host side validation.
        # Evaluating straight into `F` launches the backward term kernel with a
        # host array argument, which does not compile.
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
    # frequency, which is a statement about the junction being in the zero
    # voltage state. The direct current each junction carries is the mean of
    # the sine of its branch flux over the reconstructed period, which is
    # the pointwise sine the residual caches: setting the point is a
    # transform and the sine a pass over the grid, and neither the backward
    # transform nor the linear term is needed to read it.
    if converged && !isempty(Ljb.nzind)
        pointset || setpoint!(sys, x)
        _ensuresin!(sys)
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
                throw(ArgumentError(lazy"Source mode $(mode) is not among the retained modes; the truncation (`Nharmonics`, `maxintermodorder`, `dc`, `odd`, `even`, `frequencywindow`) removed it."))
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
