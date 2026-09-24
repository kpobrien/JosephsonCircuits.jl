# =====================================================================
# A reusable harmonic balance solver for parameter sweeps and optimizer
# loops, over a typed circuit whose values are written in terms of
# parameters.
#
# The expensive, reusable part of a solve is fixed by the topology and the
# harmonic selection: the compiled circuit, the mode grid and the Fourier
# index maps. The cheap part is fixed by the component values.
# `HBCache` holds the first and recomputes the second from the circuit's
# values at the point's definitions, so a loop over parameter values pays
# the compilation once.
#
# The cache also carries the one piece of state worth keeping between
# solves: the previously converged operating point, which warm starts the
# next one. A cold solve of a driven line can take many Newton iterations
# where a warm start from a nearby solution takes a few.
# =====================================================================

"""
    HBCache

A reusable harmonic balance solver over a typed circuit with parameters:
the compiled circuit and its definitions, the mode grid and its Fourier
index maps, the solver options, and the last converged operating point,
which [`hbsolve!`](@ref) uses to warm start the next solve.

Built by [`hbcache`](@ref). `converged` reports whether the last solve
succeeded, and a solve which does not converge also warns with the reason
it stopped. Check it: a solve which does not converge returns a state that
looks like a solution and is not one, and comparing timings or gradients
against it is meaningless.
"""
mutable struct HBCache{N,K,P}
    compiled::CompiledCircuit
    # the definitions of every parameter, which a point overrides, and the
    # keys of the definitions under each parameter name
    definitions::Dict{Any,Any}
    definitionkeys::Dict{Symbol,Vector{Any}}
    plan::P
    structure::Any
    frequencies::Frequencies{N}
    indices::FourierIndices{N}
    Nmodes::Int
    w::NTuple{N,Float64}
    sources::Vector{SourceTuple{N}}
    kwargs::K
    x::Union{Nothing,Vector{Complex{Float64}}}
    converged::Bool
    nsolves::Int
    # the system, the preconditioner and the Krylov vectors of the last
    # solve, rebound to each new point rather than rebuilt; see `HBReuse`
    reuse::HBReuse
    # the circuit matrices of the last point, refilled at the next
    nm::Union{Nothing,CircuitMatrices}
    matrixworkspace::Union{Nothing,CircuitMatrixWorkspace}
end

"""
    hbcache(w, Nharmonics, sources, circuit, circuitdefs = Dict{Symbol,Any}();
        dc = false, odd = true, even = false, maxintermodorder = Inf,
        Nevaluationharmonics = map(i -> 2i, Nharmonics),
        frequencywindow = (0, Inf), kwargs...)

A reusable nonlinear solver over a typed [`Circuit`](@ref), or the
[`CompiledCircuit`](@ref) of one, whose values are written in terms of
parameters (symbols, or the parameters of [`@params`](@ref)), with
`circuitdefs` giving every parameter a number. [`hbsolve!`](@ref) takes
the parameters to move as a named tuple keyed by their names, reads the
rest from `circuitdefs`, and evaluates the values at the point; a value
which is not a number at the definitions (a frequency dependent one) is
refused.

The harmonic selection keywords match [`hbnlsolve`](@ref); the remaining
keywords are stored and forwarded to every solve as keywords of
`hbnlsolve` on the compiled circuit, so they are the solver keywords
(`method`, `atol`, `rtol`, `iterations`, `backend`, ...) and are validated
here: a keyword the compiled circuit solve does not accept is an
`ArgumentError` at construction rather than a failure at the first solve,
as are `x0` and `reuse`, which the cache manages itself (the warm start
through `warmstart`, the reuse object internally), `keyedarrays = true`,
since the state is kept as plain vectors (`false` is accepted as what the
cache does anyway), and `method = Staged()`, which the cache does not
support since the continuation builds its own systems at its own
truncations.

# Examples
```julia
circuit = Circuit([(:P1, 1, 0, Port(1; Z0 = 50.0)), (:C1, 1, 2, Capacitor(:Cc)),
    (:Lj1, 2, 0, JosephsonJunction(:Lj)), (:C2, 2, 0, Capacitor(1000e-15))])
cache = hbcache((2*pi*4.75e9,), (8,),
    [(mode=(1,), port=1, current=1e-8)], circuit,
    Dict(:Lj => 1000e-12, :Cc => 100e-15))
for Lj in (900:25:1100)*1e-12
    sol = hbsolve!(cache, (Lj = Lj,))
    cache.converged || break
end
```
"""
function hbcache(w::NTuple{N,Number}, Nharmonics::NTuple{N,Int}, sources,
        circuit::CompilableCircuit,
        circuitdefs::AbstractDict = Dict{Symbol,Any}();
        Nevaluationharmonics::NTuple{N,Int} = map(i -> 2i, Nharmonics),
        maxintermodorder = Inf, frequencywindow = (0, Inf),
        dc::Bool = false, odd::Bool = true, even::Bool = false, kwargs...) where {N}
    all(map(>=, Nevaluationharmonics, Nharmonics)) || throw(ArgumentError(
        lazy"`Nevaluationharmonics` = $(Nevaluationharmonics) must be at least `Nharmonics` = $(Nharmonics) in every tone."))
    checkcachekwargs(kwargs)
    compiled = compile(circuit)
    definitions = definitiontable(circuitdefs)
    bound = bindvalues(compiled, cachevalues(compiled, definitions))
    # the inputs in their canonical forms, once, for every solve
    w = tonefrequencies(w)
    sources = sourcetable(sources, w)
    frequencies = removeconjfreqs(
        truncfreqs(calcfreqsrdft(Nevaluationharmonics); dc = dc, odd = odd,
            even = even, maxintermodorder = maxintermodorder,
            maxharmonics = Nharmonics, w = w,
            frequencywindow = frequencywindow))
    indices = fourierindices(frequencies)
    Nmodes = length(frequencies.modes)
    plan = circuitmatrixplan(compiled; Nmodes = Nmodes)
    return HBCache(compiled, definitions, definitionkeys(definitions), plan,
        structuralkey(bound),
        # as a named tuple: a keyword splat of mixed value types is a
        # `Pairs{Symbol,Any}`, and splatting that into every solve hands the
        # solver keywords of unknown type
        frequencies, indices, Nmodes, w, sources, NamedTuple(kwargs),
        nothing, false, 0, HBReuse(), nothing, nothing)
end

# The values of the circuit at its definitions, every one a number, as one
# real or complex vector, which is how the assembly takes them.
function cachevalues(compiled::CompiledCircuit, definitions)
    vvn = numericvalues(compiled, definitions)
    vals = Complex{Float64}[v for v in vvn]
    return all(v -> iszero(imag(v)), vals) ? real.(vals) : vals
end

# the keys of the definitions under each parameter name: a parameter may
# be defined under its parameter object, its symbol or its string, and a
# point moves every key of its name
function definitionkeys(definitions::AbstractDict)
    index = Dict{Symbol,Vector{Any}}()
    for key in keys(definitions)
        name = definitionname(key)
        isnothing(name) && continue
        push!(get!(Vector{Any}, index, name), key)
    end
    return index
end

# the definitions with the parameters of the point `p` moved, under every
# key of each name (see `definitionkeys`); a parameter the definitions do
# not hold is added under its symbol
function definitionsat(definitions::AbstractDict, index::AbstractDict,
        p::NamedTuple)
    d = copy(definitions)
    for (name, value) in zip(keys(p), values(p))
        moved = get(index, name, nothing)
        if isnothing(moved)
            d[name] = value
        else
            for key in moved
                d[key] = value
            end
        end
    end
    return d
end
# the keywords the compiled circuit solve accepts, read off its method so
# the check cannot drift from the signature
function compiledsolvekwargs()
    m = which(hbnlsolve, (NTuple{1,Float64}, Vector{SourceTuple{1}},
        Frequencies{1}, FourierIndices{1}, CompiledCircuit,
        CircuitMatrices))
    return Base.kwarg_decl(m)
end

"""
    checkcachekwargs(kwargs)

Validate the solver keywords an [`hbcache`](@ref) stores for every solve.
`x0` and `reuse` are the cache's own to manage and are refused, as is
`keyedarrays = true`, since the state is kept as plain vectors (`false` is
accepted); so is `method = Staged()`, which the compiled circuit solve does
not take; and so is any keyword that solve does not accept, which would
otherwise fail at the first [`hbsolve!`](@ref) with a method error.
"""
function checkcachekwargs(kwargs)
    for k in (:x0, :reuse)
        haskey(kwargs, k) && throw(ArgumentError(
            lazy"`$(k)` is managed by the cache and cannot be stored in it: the warm start is `hbsolve!`'s `warmstart` and the reuse object is the cache's own."))
    end
    if haskey(kwargs, :keyedarrays) && kwargs[:keyedarrays]
        throw(ArgumentError(
            "`keyedarrays = true` cannot be stored in the cache, whose state is kept as plain vectors for the warm start; index the returned arrays by position, or convert them."))
    end
    if haskey(kwargs, :method) && kwargs[:method] isa Staged
        throw(ArgumentError(
            "`method = Staged()` is not supported through the cache, since the continuation builds its own systems at its own truncations; solve with `hbnlsolve` instead."))
    end
    allowed = compiledsolvekwargs()
    for k in keys(kwargs)
        k in allowed || throw(ArgumentError(
            lazy"`$(k)` is not a keyword of the compiled circuit solve the cache runs; the solver keywords are $(allowed)."))
    end
    return nothing
end

"""
    componentvalues(cache::HBCache, p::NamedTuple)

The component values of the cache's circuit at `p`, in the compiled order:
the circuit's values at its definitions with the parameters of `p` moved.
"""
componentvalues(cache::HBCache, p::NamedTuple) =
    cachevalues(cache.compiled,
        definitionsat(cache.definitions, cache.definitionkeys, p))

"""
    reset!(cache::HBCache)

Discard the stored operating point, so the next [`hbsolve!`](@ref) starts
cold. Use this when the parameters move far enough that the previous
solution is a worse starting point than zero, or when crossing to a
different solution branch.
"""
function reset!(cache::HBCache)
    cache.x = nothing
    cache.converged = false
    return cache
end

"""
    hbsolve!(cache::HBCache, p::NamedTuple; warmstart = true)

Solve the nonlinear harmonic balance problem of `cache` at the design
parameters `p`, warm starting from the previously converged operating
point. Returns the [`NonlinearHB`](@ref) solution; `cache.converged`
reports whether it converged.

The compiled circuit and the mode grid are reused, and so are the system,
the preconditioner and the Krylov vectors of the previous solve, rebound to
the new component values (see [`HBReuse`](@ref)); only the numeric matrices
and the solve itself are recomputed. If the previous solve did not converge
its state is not used, because starting from a non-solution is usually
worse than starting cold; `warmstart = false` starts cold without
discarding the stored point, unlike [`reset!`](@ref). A component value
which crosses a structural boundary (an inductance open or shorted, a value
turned complex, a mutual coupling reaching one) invalidates the cached
sparsity patterns and is an `ArgumentError`; build a new cache for those
parameters.
"""
function hbsolve!(cache::HBCache, p::NamedTuple; warmstart::Bool = true)
    vvn = componentvalues(cache, p)
    # only the numbers moved, so the topology, the groups and the sparsity
    # patterns are reused and the matrices are refilled rather than rebuilt.
    # A value which crosses a structural boundary -- an inductance going
    # open or shorted, a capacitance going complex -- would change the
    # patterns, so it is refused rather than silently assembled against a
    # stale plan.
    bound = bindvalues(cache.compiled, vvn)
    if structuralkey(bound) != cache.structure
        throw(ArgumentError("a component value crossed a structural boundary (an inductance became open or shorted, a value became complex, or a mutual coupling reached one), so the cached sparsity patterns no longer apply. Build a new cache for these parameters."))
    end
    # into the storage of the previous point's matrices, once there are any
    nm = if isnothing(cache.nm)
        matrices = assemblematrices(cache.plan, bound)
        cache.matrixworkspace = CircuitMatrixWorkspace(cache.plan, matrices)
        matrices
    else
        matrices = assemblematrices!(cache.nm, cache.plan, bound, cache.matrixworkspace)
        if eltype(matrices.Mb) !== eltype(cache.nm.Mb)
            cache.matrixworkspace = CircuitMatrixWorkspace(cache.plan, matrices)
        end
        matrices
    end
    cache.nm = nm
    x0 = (warmstart && cache.converged) ? initialguess(cache.x) : ComplexF64[]
    # keyed arrays are a presentation convenience and pure overhead in a
    # loop; the stored state has to be a plain vector for the warm start
    nl = hbnlsolve(cache.w, cache.sources, cache.frequencies,
        cache.indices, cache.compiled, nm;
        x0 = x0, keyedarrays = false, reuse = cache.reuse, cache.kwargs...)
    cache.x = vec(collect(nl.nodeflux))
    cache.converged = nl.solverinfo.converged
    cache.nsolves += 1
    return nl
end
