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
# where a warm start from a nearby solution takes a few; a start the
# parameters jumped away from can lie outside Newton's basin, so a warm
# start which fails is retried once, cold.
# =====================================================================

"""
    HBCache

A reusable harmonic balance solver over a typed circuit with parameters:
the compiled circuit and a copy of its definitions, the mode grid and its
Fourier index maps, the solver options, and the last converged operating
point, which [`hbsolve!`](@ref) uses to warm start the next solve.

Built by [`hbcache`](@ref). `converged` reports whether the last point
converged, from its warm start or from the cold start a failed warm start
is retried from (see [`hbsolve!`](@ref)), and a point which does not
converge also warns with the reason it stopped and leaves the stored
operating point as it was. Check it: a solve which does not converge
returns a state that looks like a solution and is not one, and comparing
timings or gradients against it is meaningless. A cache is mutable and
every solve through it rewrites what it holds, so it belongs to one solve
at a time: do not share a cache between concurrent solves; build one per
task.
"""
mutable struct HBCache{N,K,P}
    compiled::CompiledCircuit
    # the definitions of every parameter, which a point overrides, and the
    # keys of the definitions under each parameter name
    definitions::Dict{Any,Any}
    definitionkeys::Dict{Symbol,Vector{Any}}
    plan::P
    frequencies::Frequencies{N}
    indices::FourierIndices{N}
    Nmodes::Int
    w::NTuple{N,Float64}
    sources::Vector{SourceTuple{N}}
    kwargs::K
    # the node fluxes of the last converged point, the warm start, and
    # whether the last solve converged
    x::Union{Nothing,Vector{Complex{Float64}}}
    converged::Bool
    nsolves::Int
    # the system of the last solve and what its method solves with,
    # rebound to each new point rather than rebuilt; see `HBReuse`
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
`circuitdefs` giving every parameter a number. The cache keeps a copy of
`circuitdefs`, so editing the dictionary afterwards does not move it.
[`hbsolve!`](@ref) takes the parameters to move as a named tuple keyed by
their names, reads the rest from the definitions, and evaluates the values
at the point; a name the definitions do not hold is refused, and so is a
value which is not a number at the definitions (a frequency dependent
one).

The harmonic selection keywords match [`hbnlsolve`](@ref); the remaining
keywords are stored and forwarded to every solve as keywords of
`hbnlsolve` on the compiled circuit, so they are the solver keywords
(`method`, `atol`, `rtol`, `iterations`, `backend`, ...) and are validated
here: a keyword the compiled circuit solve does not accept is an
`ArgumentError` at construction rather than a failure at the first solve,
as are `x0` and `reuse`, which the cache manages itself (the warm start
through `warmstart`, the reuse object internally), `keyedarrays = true`,
since the state is kept as plain vectors (`false` is accepted as what the
cache does anyway), `returnsystem = true` and `debugJacobian = true`,
which return the system rather than solve it, and `method = Staged()`,
which the cache does not support since the continuation builds its own
systems at its own truncations.

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
    checkcachekwargs(kwargs)
    compiled = compile(circuit)
    # a copy, so that the caller's dictionary does not move the cache
    definitions = Dict{Any,Any}(circuitdefs)
    # every value a number at the definitions, or refused naming it
    cachevalues(compiled, definitions)
    # the inputs in their canonical forms, once, for every solve
    w = tonefrequencies(w)
    sources = sourcetable(sources, w)
    frequencies, indices = pumpmodeset(w, Nharmonics, Nevaluationharmonics;
        dc = dc, odd = odd, even = even, maxintermodorder = maxintermodorder,
        frequencywindow = frequencywindow)
    Nmodes = length(frequencies.modes)
    plan = circuitmatrixplan(compiled; Nmodes = Nmodes)
    return HBCache(compiled, definitions, definitionkeys(definitions), plan,
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
# key of each name (see `definitionkeys`); a name the definitions do not
# hold is refused, since nothing in the circuit could read it
function definitionsat(definitions::AbstractDict, index::AbstractDict,
        p::NamedTuple)
    d = copy(definitions)
    for (name, value) in zip(keys(p), values(p))
        moved = get(index, name, nothing)
        isnothing(moved) && throwunknownparameter(name, keys(index))
        for key in moved
            d[key] = value
        end
    end
    return d
end

function throwunknownparameter(name, known)
    names = join(sort!(collect(known)), ", ")
    throw(ArgumentError(lazy"the point moves `$(name)`, which the definitions do not hold; the parameters are $(names)."))
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
accepted); so are `returnsystem = true` and `debugJacobian = true`, whose
solves return the system rather than a solution; so is `method =
Staged()`, which the compiled circuit solve does not take; and so is any
keyword that solve does not accept, which would otherwise fail at the
first [`hbsolve!`](@ref) with a method error.
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
    for k in (:returnsystem, :debugJacobian)
        haskey(kwargs, k) && kwargs[k] == true && throw(ArgumentError(
            lazy"`$(k) = true` returns the system rather than solving it, which a cache is for; build it with `hbnlsolve` or `hbnonlinearproblem`."))
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

# whether the matrices `nm` hold the element types the values of `b`
# assemble to, so that they can be refilled in place
function refillable(nm::CircuitMatrices, b::BoundCircuit)
    return eltype(nm.Cnm) === eltype(b.capacitors) &&
        eltype(nm.Gnm) === eltype(b.resistors) &&
        eltype(nm.Lb) === eltype(b.inductors) &&
        eltype(nm.Ljb) === eltype(b.junctions) &&
        eltype(nm.Mb) === promote_type(eltype(b.inductors),
            eltype(b.mutualinductors))
end

"""
    reset!(cache::HBCache)

Discard the stored operating point, so the next [`hbsolve!`](@ref) starts
cold, and with it what the previous solves taught the preconditioner: a
coupling set grown by escalation or by measurement, and the deflation
candidates of a [`Floquet`](@ref) preconditioner, so the next solve builds
the preconditioner its method asks for. A warm start which fails is
retried cold without discarding any of this (see [`hbsolve!`](@ref)); use
it when the parameters jump far enough that the previous solution is a
worse starting point than zero, which saves that failed attempt, when
crossing to a different solution branch, or to drop what the
preconditioner grew.
"""
function reset!(cache::HBCache)
    cache.x = nothing
    cache.converged = false
    cache.reuse.preconditioner = nothing
    cache.reuse.recycling = nothing
    return cache
end

"""
    hbsolve!(cache::HBCache, p::NamedTuple; warmstart = true)

Solve the nonlinear harmonic balance problem of `cache` at the design
parameters `p`, warm starting from the last converged operating point.
Returns the [`NonlinearHB`](@ref) solution; `cache.converged` reports
whether it converged. A name of `p` which the definitions of the
cache do not hold is an `ArgumentError`.

The compiled circuit and the mode grid are reused, and so are the system
of the previous solve and what its method solves with, rebound to the new
component values (see [`HBReuse`](@ref)): under `NewtonKrylov` the
preconditioner and the Krylov vectors, under `Newton` and `QuasiNewton`
the assembled Jacobian and its factorization, whose fill reducing ordering
and symbolic analysis are kept. Only the numeric matrices and the solve
itself are recomputed. The matrices are refilled on the
patterns of the compiled circuit, which do not depend on the values, and
assembled anew when a value changes the element type of its group (a
resistance or a capacitance turned complex).

A start the parameters jumped away from can lie outside Newton's basin,
so a warm started solve which does not converge is retried once from a
cold start, the solve `warmstart = false` makes, keeping what the cache
reuses. The warm attempt's messages reach the caller only when it is the
outcome, so a retried point warns as its cold solve does, and a point
costs at most the two solves, each within the solver's own budget. A
retried point returns the cold solve, with the warm attempt's record
ahead of its own in `solverinfo.stages`; where the circuit has several
operating points, it is the one a cold start reaches, which need not be
on the branch the sweep followed. A retry which converges becomes the
stored point, and a point which fails both ways leaves the stored point
as it was, so the next one starts from the last solution rather than from
a non-solution or from nothing. `warmstart = false` starts cold, is not
retried, and keeps the stored point, unlike
[`JosephsonCircuits.reset!`](@ref).
"""
function hbsolve!(cache::HBCache, p::NamedTuple; warmstart::Bool = true)
    vvn = componentvalues(cache, p)
    # only the numbers moved, so the topology, the groups and the sparsity
    # patterns are reused, and the matrices are refilled in the storage of
    # the previous point's when they hold the element types the values
    # assemble to, and assembled anew otherwise
    bound = bindvalues(cache.compiled, vvn)
    nm = if isnothing(cache.nm) || !refillable(cache.nm, bound)
        matrices = assemblematrices(cache.plan, bound)
        cache.matrixworkspace = CircuitMatrixWorkspace(cache.plan, matrices)
        matrices
    else
        assemblematrices!(cache.nm, cache.plan, bound, cache.matrixworkspace)
    end
    cache.nm = nm
    nl = if warmstart && !isnothing(cache.x)
        # a jump of the parameters can leave the warm start outside
        # Newton's basin, so a warm attempt which fails is retried once,
        # cold, with what the cache reuses. Its messages are held until it
        # is known to be the outcome, so that only the outcome's reach the
        # caller.
        held = HeldMessages(Base.CoreLogging.current_logger())
        warm = Base.CoreLogging.with_logger(held) do
            cachesolve(cache, nm, initialguess(cache.x))
        end
        if warm.solverinfo.converged
            release!(held)
            warm
        else
            # the record holds both attempts, the warm one first
            cold = cachesolve(cache, nm, ComplexF64[])
            prepend!(cold.solverinfo.stages, warm.solverinfo.stages)
            cold
        end
    else
        cachesolve(cache, nm, ComplexF64[])
    end
    cache.converged = nl.solverinfo.converged
    cache.converged && (cache.x = vec(collect(nl.nodeflux)))
    cache.nsolves += 1
    return nl
end

# A solve of the cache's problem with the matrices `nm` from the start
# `x0`, cold when it is empty. Keyed arrays are a presentation convenience
# and pure overhead in a loop, and the stored state has to be a plain
# vector for the warm start.
cachesolve(cache::HBCache, nm::CircuitMatrices, x0::Vector{ComplexF64}) =
    hbnlsolve(cache.w, cache.sources, cache.frequencies, cache.indices,
        cache.compiled, nm; x0 = x0, keyedarrays = false,
        reuse = cache.reuse, cache.kwargs...)

# The log messages of an attempt held until it is known to be the outcome
# of its point: released to the logger they were meant for if it is, and
# dropped if a retry replaces it. It logs what that logger would, and is
# read only when a message is logged, so the logger is held untyped.
struct HeldMessages <: Base.CoreLogging.AbstractLogger
    logger::Base.CoreLogging.AbstractLogger
    messages::Vector{Any}
end
HeldMessages(logger::Base.CoreLogging.AbstractLogger) =
    HeldMessages(logger, Any[])
Base.CoreLogging.min_enabled_level(h::HeldMessages) =
    Base.CoreLogging.min_enabled_level(h.logger)
Base.CoreLogging.shouldlog(h::HeldMessages, args...) =
    Base.CoreLogging.shouldlog(h.logger, args...)
Base.CoreLogging.catch_exceptions(h::HeldMessages) =
    Base.CoreLogging.catch_exceptions(h.logger)
function Base.CoreLogging.handle_message(h::HeldMessages, args...; kwargs...)
    push!(h.messages, (args, kwargs))
    return nothing
end
function release!(h::HeldMessages)
    for (args, kwargs) in h.messages
        Base.CoreLogging.handle_message(h.logger, args...; kwargs...)
    end
    return nothing
end
