# Beyn's block-moment contour reduction (LAA 436, 2012, 3839–3863,
# doi:10.1016/j.laa.2011.03.030). Moments use the unit circle coordinate
# (s-center)/radius, avoiding powers of a GHz-sized physical frequency.
# One circuit-wide CSC pattern, including every local response entry even
# when it vanishes at a sampled frequency. All destinations are computed
# once; contour evaluation changes only nzval. A workspace belongs to one
# search, so independent searches can safely use the same pole system.
struct PoleMatrixWorkspace
    matrix::SparseMatrixCSC{ComplexF64,Int}
    q0::Vector{ComplexF64}
    q1::Vector{ComplexF64}
    q2::Vector{ComplexF64}
    destinations::Vector{Matrix{Int}}
end

function PoleMatrixWorkspace(sys)
    n = size(sys.Q0, 1)
    rows, cols = Int[], Int[]
    for Q in (sys.Q0, sys.Q1, sys.Q2), col in 1:n, k in nzrange(Q, col)
        push!(rows, rowvals(Q)[k]); push!(cols, col)
    end
    for term in sys.terms, col in term.cols, row in term.rows
        push!(rows, row); push!(cols, col)
    end
    Q = sparse(rows, cols, zeros(ComplexF64, length(rows)), n, n)
    function destination(row, col)
        range = nzrange(Q, col)
        return first(range) - 1 + searchsortedfirst(view(rowvals(Q), range), row)
    end
    function align(A)
        values = zeros(ComplexF64, nnz(Q))
        for col in 1:n, k in nzrange(A, col)
            values[destination(rowvals(A)[k], col)] = nonzeros(A)[k]
        end
        return values
    end
    destinations = [[destination(row, col) for row in term.rows, col in term.cols]
        for term in sys.terms]
    return PoleMatrixWorkspace(Q, align(sys.Q0), align(sys.Q1), align(sys.Q2), destinations)
end

function polematrix!(work::PoleMatrixWorkspace, sys, z)
    Q = work.matrix
    nonzeros(Q) .= work.q0 .+ z .* work.q1 .+ z^2 .* work.q2
    for (term, destinations) in zip(sys.terms, work.destinations)
        localmatrix = poleterm(term, z, sys.scale)
        for (k, value) in zip(destinations, localmatrix)
            nonzeros(Q)[k] += value
        end
    end
    all(isfinite, nonzeros(Q)) || throw(ArgumentError("the pole operator is nonfinite at s = $(z*sys.scale)."))
    return Q
end

polematrix(sys, z) = polematrix!(PoleMatrixWorkspace(sys), sys, z)

# The contributions of the nodes `nodes` of the unit circle to the moment
# sums `sums[k]`, the sum over the nodes `t` of `t^(k-1) X_t`, `X_t` the
# solve of the probes `V` at the point `center + radius*t` times
# `radius*t`, and to `bound`, the sum of their norms, returned. The sums
# hold every node a rule has reached, so a doubled rule adds only the
# nodes between, and the rule of `N` nodes is the sums over `N`; every
# point refactors the one pattern of `work` in `cache`.
function addmoments!(sums, bound, work, cache, factorization, sys, center, radius, nodes, V)
    for t in nodes
        tryfactorize!(cache, factorization, polematrix!(work, sys, center + radius*t))
        X = (cache.factorization \ V) .* (radius*t)
        bound += norm(X)
        power = one(ComplexF64)
        for k in eachindex(sums)
            sums[k] .+= power .* X
            power *= t
        end
    end
    return bound
end

# The poles inside the circle from the trapezoidal rule of `N` nodes: the
# block Hankel matrices of the moments, the first's rank cut at `ranktol`
# of its size or of the rule's bound, and the eigenvalues of the reduced
# pencil inside the unit circle, those outside being artifacts of the
# quadrature
function polereduce(sums, bound, N, K, ranktol, center, radius)
    A = [s ./ N for s in sums]
    n, l = size(first(A))
    H0, H1 = zeros(ComplexF64, n*K, l*K), zeros(ComplexF64, n*K, l*K)
    for i in 1:K, j in 1:K
        H0[(i-1)*n+1:i*n, (j-1)*l+1:j*l] = A[i+j-1]
        H1[(i-1)*n+1:i*n, (j-1)*l+1:j*l] = A[i+j]
    end
    S = svd(H0)
    rank = count(>(ranktol*max(first(S.S), bound/N, floatmin(Float64))), S.S)
    capacity = min(n*K, l*K)
    rank == 0 && return (; values = ComplexF64[], vectors = zeros(ComplexF64, n, 0),
        rank, capacity, moments = A, bound = bound/N)
    U, W = S.U[:, 1:rank], S.V[:, 1:rank]
    eig = eigen((U'*H1*W)*Diagonal(inv.(S.S[1:rank])))
    inside = findall(z -> isfinite(z) && abs(z) < 1, eig.values)
    return (; values = center .+ radius .* eig.values[inside],
        vectors = U[1:n, :]*eig.vectors[:, inside], rank, capacity, moments = A,
        bound = bound/N)
end

# the poles inside a circle, by the rule of `quadrature` nodes doubled
# until two rules agree (see ContourIntegral)
function polesolve(method::ContourIntegral, sys, factorization)
    n = size(sys.Q0, 1)
    l = min(n, method.probes)
    # every unit vector for a small reference; otherwise probes from a
    # fixed sequence, which leaves the caller's random stream alone
    V = l == n ? Matrix{ComplexF64}(I, n, n) :
        ComplexF64[sin(i*j*sqrt(2)) + im*cos(i*j*sqrt(3)) for i in 1:n, j in 1:l]/sqrt(n)
    center, radius = method.center/sys.scale, method.radius/sys.scale
    K, N = method.moments, method.quadrature
    work, cache = PoleMatrixWorkspace(sys), FactorizationCache()
    sums = [zeros(ComplexF64, n, l) for _ in 1:2K]
    bound = addmoments!(sums, 0.0, work, cache, factorization, sys, center, radius,
        (cis(2pi*j/N) for j in 0:N-1), V)
    previous = polereduce(sums, bound, N, K, method.ranktol, center, radius)
    distance(a, b) = isempty(a) ? 0.0 : isempty(b) ? Inf :
        maximum(x -> minimum(abs.(x .- b)), a)/radius
    result = nothing
    for _ in 1:method.refinements
        bound = addmoments!(sums, bound, work, cache, factorization, sys, center, radius,
            (cis(pi*(2j + 1)/N) for j in 0:N-1), V)
        N *= 2
        current = polereduce(sums, bound, N, K, method.ranktol, center, radius)
        change = max(distance(previous.values, current.values), distance(current.values, previous.values))
        momentchange = maximum(norm(a - b) for (a, b) in zip(current.moments, previous.moments))/
            max(current.bound, previous.bound, floatmin(Float64))
        saturated = current.rank == current.capacity
        samecount = current.rank == previous.rank && length(current.values) == length(previous.values)
        settled = samecount && change <= method.changetol && momentchange <= method.changetol
        converged = settled && length(current.values) == current.rank && !saturated
        boundarydistance = isempty(current.values) ? Inf :
            minimum(radius .- abs.(current.values .- center))*sys.scale
        diagnostics = (; quadrature = N, rank = current.rank, capacity = current.capacity,
            saturated, change, momentchange, boundarydistance)
        result = poleresult(sys, current.values, current.vectors,
            zeros(Int, length(current.values)), NamedTuple[diagnostics], converged,
            nothing, method.tol, method)
        result.converged && return result
        # more nodes cannot raise the rank's capacity: report it at once,
        # for more probes or moments
        settled && saturated && return result
        previous = current
    end
    return result
end
