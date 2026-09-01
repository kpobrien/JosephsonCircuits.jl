# =====================================================================
# A reusable harmonic balance solver for parameter sweeps and optimizer
# loops over circuit builders.
#
# The expensive, reusable part of a solve is fixed by the topology and the
# harmonic selection: the parse, the graph, the mode grid and the Fourier
# index maps. The cheap part is fixed by the component values. `HBCache`
# holds the first and recomputes the second from the builder, so a loop
# over parameter values pays the parse once.
#
# The cache also carries the one piece of state worth keeping between
# solves: the previously converged operating point, used to warm start the
# next one. That is where the time actually goes -- a cold solve of a
# driven line can take many Newton iterations, while warm started from a
# nearby solution it takes a few.
# =====================================================================

"""
    HBCache

A reusable harmonic balance solver over a circuit builder: the parsed
sorted circuit, the circuit graph, the mode grid and its Fourier index
maps, the solver options, and the last converged operating point, which
[`hbsolve!`](@ref) uses to warm start the next solve.

Built by [`hbcache`](@ref). `converged` reports whether the last solve
succeeded. Check it: a solve which does not converge returns a state that
looks like a solution and is not one, and comparing timings or gradients
against it is meaningless.
"""
mutable struct HBCache{N,K,P}
    builder::Any
    compiled::CompiledCircuit
    plan::P
    structure::Any
    valueorder::Vector{Int}
    cg::CircuitGraph
    frequencies::Frequencies{N}
    indices::FourierIndices{N}
    Nmodes::Int
    w::NTuple{N,Number}
    sources::Vector
    kwargs::K
    x::Union{Nothing,Vector{Complex{Float64}}}
    converged::Bool
    nsolves::Int
end

"""
    hbcache(w, Nharmonics, sources, builder, p::NamedTuple;
        dc = false, odd = true, even = false, maxintermodorder = Inf,
        maxharmonics = Nharmonics, sorting = :number, kwargs...)

A reusable nonlinear solver over the circuit builder `builder`, parsed
once at the parameter point `p`. The harmonic selection keywords match
[`hbnlsolve`](@ref); the remaining keywords are stored and forwarded to
every solve.

The builder must keep the circuit topology fixed as the parameters vary;
[`hbsolve!`](@ref) checks the component names of each new evaluation
against the parse.

# Examples
```julia
make(; Lj, Cc) = [("P1","1","0",1), ("R1","1","0",50.0),
    ("C1","1","2",Cc), ("Lj1","2","0",Lj), ("C2","2","0",1000e-15)]
cache = hbcache((2*pi*4.75e9,), (8,),
    [(mode=(1,), port=1, current=1e-8)], make,
    (Lj = 1000e-12, Cc = 100e-15))
for Lj in (900:25:1100)*1e-12
    sol = hbsolve!(cache, (Lj = Lj, Cc = 100e-15))
    cache.converged || break
end
```
"""
function hbcache(w::NTuple{N,Number}, Nharmonics::NTuple{N,Int}, sources,
        builder, p::NamedTuple;
        maxharmonics::NTuple{N,Int} = Nharmonics, maxintermodorder = Inf,
        dc::Bool = false, odd::Bool = true, even::Bool = false,
        sorting = :number, kwargs...) where {N}

    frequencies = removeconjfreqs(
        truncfreqs(calcfreqsrdft(Nharmonics); dc = dc, odd = odd,
            even = even, maxintermodorder = maxintermodorder,
            maxharmonics = maxharmonics))
    indices = fourierindices(frequencies)
    Nmodes = length(frequencies.modes)

    # A netlist of tuples, not a typed `Circuit`: the value table below is
    # built by walking the builder's output entry by entry, which a
    # `Circuit` does not support, and a typed circuit's compiled table can
    # hold generated entries the builder never returned. Rebinding a typed
    # circuit needs the compiler's own value plan; until it has one, this
    # takes what it can consume.
    circuit0 = builder(; p...)
    circuit0 isa AbstractVector || throw(ArgumentError(lazy"hbcache needs a builder returning a netlist of (name, node1, node2, value) tuples; this one returned a $(typeof(circuit0)). A typed `Circuit` is not supported here yet, because the cache rebinds values by walking the builder's output."))
    compiled = compile(circuit0; sorting = sorting)
    # the loop enumeration is quadratic in the number of inductive
    # loops and nothing here reads it
    cg = calccircuitgraph(compiled; loops = false)
    bound = bind(compiled)
    plan = circuitmatrixplan(compiled, cg, bound; Nmodes = Nmodes)

    # where each component the builder returns lands in the value table.
    # The topology is fixed as the parameters vary, so this is fixed too, and
    # looking it up once turns a string hash and a dictionary probe per
    # component per solve into an array read.
    valueorder = [get(compiled.componentnamedict, String(first(c)), 0)
                  for c in circuit0]

    return HBCache(builder, compiled, plan, structuralkey(bound), valueorder,
        cg, frequencies, indices, Nmodes, w, collect(sources), kwargs,
        nothing, false, 0)
end

"""
    componentvalues(cache::HBCache, p::NamedTuple)

The component values of `cache.builder` at `p`, in the parsed sorted
order, without re-parsing. The builder output must have the same
component names as the parse and fully numeric values.
"""
function componentvalues(cache::HBCache, p::NamedTuple)
    # A function barrier. The builder is stored as `Any`, so calling it
    # yields a value of unknown type and everything downstream of it in the
    # same function is dynamically dispatched -- once per component, per
    # solve. Handing the result to a function which specializes on its
    # concrete type costs one dispatch instead of thousands.
    return gathercomponentvalues(cache.builder(; p...), cache.valueorder,
        cache.compiled.componentnames, cache.compiled.componentnamedict)
end

function gathercomponentvalues(circuit, order::Vector{Int},
        names::Vector{String}, namedict::Dict{String,Int})
    length(circuit) == length(names) ||
        throw(ArgumentError(lazy"the builder returned $(length(circuit)) components where the parse has $(length(names)); the circuit topology must be fixed as the parameters vary."))
    vals = Vector{Complex{Float64}}(undef, length(names))
    @inbounds for (k, c) in enumerate(circuit)
        name = first(c)
        # the cached position, confirmed by comparing the name rather than
        # hashing it. A builder which reorders its components between calls
        # still lands in the right place, it just pays for the lookup.
        i = (k <= length(order) && !iszero(order[k]) &&
             isequal(name, names[order[k]])) ? order[k] :
            get(namedict, String(name), 0)
        iszero(i) && throw(ArgumentError(lazy"the builder returned the component $(name), which is not in the parsed circuit; the circuit topology must be fixed as the parameters vary."))
        v = c[4]
        v isa Number || throw(ArgumentError(lazy"the component $(name) has the non-numeric value $(v); the cached solver requires a fully numeric builder output."))
        vals[i] = Complex{Float64}(v)
    end
    # keep purely real value vectors real, which the assembly prefers
    return all(v -> iszero(imag(v)), vals) ? real.(vals) : vals
end

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

The parse, the graph and the mode grid are reused; only the component
values, the numeric matrices and the solve itself are recomputed. If the
previous solve did not converge its state is not used, because starting
from a non-solution is usually worse than starting cold.
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
    nm = assemblematrices(cache.plan, bound)
    x0 = (warmstart && cache.converged) ? cache.x : nothing
    # keyed arrays are a presentation convenience and pure overhead in a
    # loop; the stored state has to be a plain vector for the warm start
    nl = hbnlsolve(cache.w, cache.sources, cache.frequencies,
        cache.indices, cache.compiled, cache.cg, nm;
        x0 = x0, keyedarrays = false, cache.kwargs...)
    cache.x = vec(collect(nl.nodeflux))
    cache.converged = nl.solverinfo.converged
    cache.nsolves += 1
    return nl
end
