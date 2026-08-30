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
mutable struct HBCache{N,K}
    builder::Any
    psc::ParsedSortedCircuit
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

    psc = parsesortcircuit(builder(; p...), sorting = sorting)
    cg = calccircuitgraph(psc)

    return HBCache(builder, psc, cg, frequencies, indices, Nmodes, w,
        collect(sources), kwargs, nothing, false, 0)
end

"""
    componentvalues(cache::HBCache, p::NamedTuple)

The component values of `cache.builder` at `p`, in the parsed sorted
order, without re-parsing. The builder output must have the same
component names as the parse and fully numeric values.
"""
function componentvalues(cache::HBCache, p::NamedTuple)
    circuit = cache.builder(; p...)
    length(circuit) == length(cache.psc.componentvalues) ||
        throw(ArgumentError(lazy"the builder returned $(length(circuit)) components where the parse has $(length(cache.psc.componentvalues)); the circuit topology must be fixed as the parameters vary."))
    vals = Vector{Complex{Float64}}(undef, length(circuit))
    for c in circuit
        name = String(first(c))
        i = get(cache.psc.componentnamedict, name, 0)
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
    nm = numericmatrices(cache.psc, cache.cg, vvn; Nmodes = cache.Nmodes)
    x0 = (warmstart && cache.converged) ? cache.x : nothing
    # keyed arrays are a presentation convenience and pure overhead in a
    # loop; the stored state has to be a plain vector for the warm start
    nl = hbnlsolve(cache.w, cache.sources, cache.frequencies,
        cache.indices, cache.psc, cache.cg, nm;
        x0 = x0, keyedarrays = false, cache.kwargs...)
    cache.x = vec(collect(nl.nodeflux))
    cache.converged = nl.solverinfo.converged
    cache.nsolves += 1
    return nl
end
