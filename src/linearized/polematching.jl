"""
    JosephsonCircuits.matchpoles(reference, candidates;
        indices = eachindex(reference.poles),
        candidateindices = eachindex(candidates.poles),
        minoverlap = 0.9, maxdistance = Inf, aliases = true)

Match selected `reference` poles to distinct `candidates` by maximizing the
total physical profile overlap, using node voltages and junction fluxes.
Works across different harmonic truncations: absent harmonics are zero
padded, with their tails retained in the normalization. Node and junction
axes must agree. Voltages on both sides use the reference frequency scale.

Returns one named tuple per selected reference pole: `index` (zero when
unmatched), `overlap`, `alias` and `distance`. With `aliases = true` and
equal nonzero pump frequencies, the candidate's nearest Floquet alias is
compared: `candidate.poles[index] - im*alias*pumpfrequency`. Its harmonic
coefficients are shifted by `alias` at the same time. `distance` is the
absolute difference of these physical rates, in inverse seconds.

`minoverlap` rejects weak matches and `maxdistance` optionally limits rate
motion. `candidateindices` selects acceptable candidates without copying a
result; returned indices still refer to the original result. Unmatched poles
are reported explicitly; a parameter continuation
should reduce its step or expand its search. An overlap is not a proof of
branch identity at a degeneracy. Remove duplicate discoveries and boundary
artifacts from the candidate set before interpreting matches. Modes with
no observable node voltage or junction flux cannot be matched this way.
No new circuit solve is launched by this function.

# Examples
```jldoctest
julia> c = Circuit([(:r, 1, 0, Resistor(2.0)), (:c, 1, 0, Capacitor(0.5))]);

julia> p = hbstability(c);

julia> matched = only(JosephsonCircuits.matchpoles(p, p; indices = [1]));

julia> matched.index == 1 && matched.overlap ≈ 1 && matched.distance == 0
true
```
"""
function matchpoles(reference::HBStabilityResult, candidates::HBStabilityResult;
        indices = eachindex(reference.poles),
        candidateindices = eachindex(candidates.poles), minoverlap::Real = 0.9,
        maxdistance::Real = Inf, aliases::Bool = true)
    reference.nodes == candidates.nodes && reference.junctionbranches == candidates.junctionbranches ||
        throw(ArgumentError("pole profiles must have identical node and junction axes."))
    0 <= minoverlap <= 1 && isfinite(minoverlap) && maxdistance >= 0 ||
        throw(ArgumentError("require 0 <= minoverlap <= 1 and maxdistance >= 0."))
    ids = collect(Int, indices)
    all(i -> 1 <= i <= length(reference.poles), ids) || throw(BoundsError(reference.poles, ids))
    cids = collect(Int, candidateindices)
    all(i -> 1 <= i <= length(candidates.poles), cids) || throw(BoundsError(candidates.poles, cids))
    allunique(cids) || throw(ArgumentError("candidateindices must be distinct."))
    nr, nc = length(ids), length(cids)
    nr == 0 && return NamedTuple{(:index,:overlap,:alias,:distance),Tuple{Int,Float64,Int,Float64}}[]
    overlap = zeros(nr, nc)
    shifts = zeros(Int, nr, nc)
    distances = fill(Inf, nr, nc)
    scale = reference.frequencyscale
    samepump = reference.pumpfrequency > 0 &&
        samefrequency(reference.pumpfrequency, candidates.pumpfrequency, reference.pumpfrequency)
    amodes = Dict(only(mode) => k for (k, mode) in enumerate(reference.modes))
    for (a, i) in enumerate(ids), (b, j) in enumerate(cids)
        shift = aliases && samepump ? round(Int, imag(candidates.poles[j]-reference.poles[i])/reference.pumpfrequency) : 0
        shifts[a,b] = shift
        distances[a,b] = abs(candidates.poles[j] - im*shift*reference.pumpfrequency - reference.poles[i])
        av, bv = view(reference.nodevoltage, :, :, i)./scale, view(candidates.nodevoltage, :, :, j)./scale
        af, bf = view(reference.junctionflux, :, :, i), view(candidates.junctionflux, :, :, j)
        denominator = sqrt((sum(abs2,av)+sum(abs2,af))*(sum(abs2,bv)+sum(abs2,bf)))
        (reference.internalonly[i] || candidates.internalonly[j] || denominator == 0) && continue
        product = 0.0+0im
        for (b, mode) in enumerate(candidates.modes)
            k = get(amodes, only(mode)+shift, 0)
            k == 0 && continue
            product += dot(view(av,k,:), view(bv,b,:)) + dot(view(af,k,:), view(bf,b,:))
        end
        overlap[a,b] = min(1.0, abs(product)/denominator)
    end
    # Dummy columns allow every reference to remain unmatched. Solve a
    # rectangular assignment rather than taking independent nearest modes.
    cost = zeros(nr, nc+nr)
    for a in 1:nr, j in 1:nc
        cost[a,j] = overlap[a,j] >= minoverlap && distances[a,j] <= maxdistance ? -overlap[a,j] : 1.0
    end
    assignment = poleassignment(cost)
    return [j <= nc && cost[a,j] < 0 ?
        (index = cids[j], overlap = overlap[a,j], alias = shifts[a,j], distance = distances[a,j]) :
        (index = 0, overlap = 0.0, alias = 0, distance = Inf)
        for (a,j) in enumerate(assignment)]
end

# Shortest augmenting path assignment, with a sentinel column at index 1.
function poleassignment(cost)
    n, m = size(cost)
    u, v = zeros(n), zeros(m+1)
    p, way = zeros(Int,m+1), zeros(Int,m+1)
    for i in 1:n
        p[1] = i
        j0 = 1
        minv, used = fill(Inf,m+1), falses(m+1)
        while true
            used[j0] = true
            i0, delta, j1 = p[j0], Inf, 0
            for j in 2:m+1
                used[j] && continue
                cur = cost[i0,j-1]-u[i0]-v[j]
                if cur < minv[j]
                    minv[j], way[j] = cur, j0
                end
                if minv[j] < delta
                    delta, j1 = minv[j], j
                end
            end
            for j in 1:m+1
                if used[j]
                    u[p[j]] += delta
                    v[j] -= delta
                else
                    minv[j] -= delta
                end
            end
            j0 = j1
            p[j0] == 0 && break
        end
        while true
            j1 = way[j0]
            p[j0] = p[j1]
            j0 = j1
            j0 == 1 && break
        end
    end
    assignment = zeros(Int,n)
    for j in 2:m+1
        p[j] > 0 && (assignment[p[j]] = j-1)
    end
    return assignment
end
