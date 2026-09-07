# The fit of a rational scattering block to sampled scattering data, by
# vector fitting: a set of common poles refined by the relaxed pole
# relocation of Gustavsen and Semlyen, the residues of every entry and the
# feedthrough by least squares at the converged poles, a real state space
# realization from the poles and residues, and a passivity enforcement by
# the smallest perturbation of the residues that brings every singular
# value crossing one back under it, until the bounded real lemma's test
# passes. Nothing here needs a package beyond linear algebra.

"""
    RationalScattering(block::ScatteringParameters, npoles;
        frequencies = nothing, iterations = 30, passivity = true, atol = 1e-8)

A [`RationalScattering`](@ref) block fitted to the scattering data of
`block`, a tabulated or Touchstone block or any other, at `npoles` common
poles by vector fitting: the poles start as complex pairs spread over the
band, are relocated by the relaxed iteration of Gustavsen and Semlyen
until they settle, and the residues of every entry and the constant term
follow by least squares. The data is sampled at `frequencies` in Hz, by
default a tabulated block's own. The poles the data does not need drift
out of the band or settle on one another and are dropped, so ask for as
many as the data might need; too few settle on a poor fit that cannot be
made passive, which is an error. The fit is a real state space realization
with the port count states per pole kept, stable by construction, and is
made passive where the fit strays above unit singular value by the
smallest change of its residues that brings the crossings back, checked by
the same test as any rational block; `passivity = false` returns the raw
fit and rejects one that is not passive. The fitted block keeps the
reference impedances, the grounding and the noise model of `block`.
A delay is not a rational function: a cable is a [`TransmissionLine`](@ref)
of its delay in cascade with a fit of the data with that delay removed.
"""
function RationalScattering(block::ScatteringParameters, npoles::Integer; frequencies = nothing,
        iterations::Integer = 30, passivity::Bool = true, atol::Real = 1e-8)
    npoles >= 1 || throw(ArgumentError("fit at least one pole."))
    fs = if isnothing(frequencies)
        block.provider isa TabulatedMatrixProvider || throw(ArgumentError(
            "give the frequencies in Hz to sample the block at; only a tabulated block has its own."))
        copy(block.provider.frequencies) ./ (2pi)
    else
        Float64.(collect(frequencies))
    end
    all(f -> isfinite(f) && f > 0, fs) && issorted(fs) || throw(ArgumentError("the sample frequencies must be positive, finite and increasing."))
    n = block.nports
    S = zeros(ComplexF64, n, n, length(fs))
    evaluatescattering!(S, block, 2pi .* fs)
    poles, residues, D = vectorfit(S, 2pi .* fs, Int(npoles), Int(iterations))
    A, B, C = realization(poles, residues, n)
    if passivity
        A, B, C, D = enforcepassivity(A, B, C, D, 2pi .* fs; atol = atol)
    end
    fitted = RationalScattering(A, B, C, D; zref = block.zref, grounded = block.grounded, noise = block.noise, atol = atol)
    return fitted
end

# The relaxed vector fit: with the poles `a`, the unknowns of every entry
# are its residues and constant, and the shared unknowns the residues
# and the constant of the weight `sigma(s) = d + sum r_p/(s - a_p)` whose
# zeros are the next poles, with the relaxed normalization that the real
# part of `sigma` averages to one over the samples so that the trivial
# solution is excluded; complex poles come in conjugate pairs and enter
# through the real basis, so every unknown is real. Returns the poles,
# the residue matrices `(n, n, npoles)` and the constant `(n, n)`.
function vectorfit(S::AbstractArray{<:Complex,3}, ws::AbstractVector, npoles::Int, iterations::Int)
    n, K = size(S, 1), length(ws)
    wmin, wmax = first(ws), last(ws)
    # the starting poles: complex pairs along the band, a real one if odd
    poles = ComplexF64[]
    npairs = npoles ÷ 2
    spread = npairs == 1 ? [(wmin + wmax)/2] : collect(range(wmin, wmax; length = npairs))
    for w in spread
        push!(poles, complex(-0.01w, w))
        push!(poles, complex(-0.01w, -w))
    end
    isodd(npoles) && push!(poles, complex(-(wmin + wmax)/2, 0.0))
    poles = converge(S, ws, poles, iterations)
    poles = prunepoles(S, ws, poles, iterations)
    residues, D = fitresidues(S, ws, poles)
    return poles, residues, D
end

# the relocation iterated until the poles settle to a part in a billion,
# or until the fit at them is at the roundoff of the least squares, where
# the weight's columns are combinations of the entries' and the
# relocation would be arbitrary
roundoff(S) = 1e-12*maximum(abs, S)
function converge(S::AbstractArray{<:Complex,3}, ws::AbstractVector, poles::Vector{ComplexF64}, iterations::Int)
    wmin = first(ws)
    order(x) = (real(x), imag(x))
    for iteration in 1:iterations
        fiterror(S, ws, poles) <= roundoff(S) && break
        newpoles = relocate(S, ws, poles)
        change = length(newpoles) == length(poles) ?
            maximum(abs.(sort(newpoles; by = order) .- sort(poles; by = order)) ./ max.(abs.(poles), wmin)) : Inf
        poles = newpoles
        change < 1e-9 && break
    end
    return poles
end

# The poles the fit does not need, once the relocation has settled: a
# pole the data does not need drifts out of the band, far above it where
# it fits a constant that the constant term holds, or below it where the
# data cannot see it, carrying next to nothing; or it settles in a
# cluster with the pole it copies, split by less than the data can
# resolve, the members carrying large residues of opposite sign. The
# poles are tried without, least contribution first, and every cluster
# of poles closer than the sample spacing at its frequency is replaced by
# its mean and relocated again from there, as the mean of a cluster is
# not the pole it stands for; each simplification is kept only if the
# residues refitted reproduce the data as closely as before, within a
# factor of two above the roundoff of the least squares, so nothing the
# data needs is removed.
function fiterror(S::AbstractArray{<:Complex,3}, ws::AbstractVector, poles::Vector{ComplexF64})
    residues, D = fitresidues(S, ws, poles)
    err = 0.0
    for (k, w) in enumerate(ws)
        F = D .+ sum(residues[:, :, p] ./ (im*w - poles[p]) for p in eachindex(poles); init = zeros(size(D)))
        err = max(err, maximum(abs, F .- S[:, :, k]))
    end
    return err
end
function prunepoles(S::AbstractArray{<:Complex,3}, ws::AbstractVector, poles::Vector{ComplexF64}, iterations::Int)
    acceptable(err, before) = err <= 2before + roundoff(S)
    before = fiterror(S, ws, poles)
    while true
        n0 = length(poles)
        poles, before = dropneedless(S, ws, poles, before, acceptable)
        poles, before = mergeclusters(S, ws, poles, iterations, before, acceptable)
        length(poles) == n0 && break
    end
    return poles
end
# a pole and its conjugate are dropped when the fit without them is as
# close, the pole contributing least to the fit tried first
function dropneedless(S, ws, poles, before, acceptable)
    while true
        residues, _ = fitresidues(S, ws, poles)
        contribution(p) = maximum(opnorm(residues[:, :, p])/abs(im*w - poles[p]) for w in ws)
        candidates = sort(filter(p -> imag(poles[p]) >= 0, eachindex(poles)); by = contribution)
        dropped = false
        for p in candidates
            keep = [q for q in eachindex(poles) if q != p && !(imag(poles[p]) != 0 && poles[q] == conj(poles[p]))]
            trial = poles[keep]
            err = fiterror(S, ws, trial)
            if acceptable(err, before)
                # the last pole is needless only if a constant reproduces
                # the data, which is then no rational block
                isempty(trial) && throw(ArgumentError("no pole is needed: a constant reproduces the data, so give it as a ScatteringParameters matrix."))
                poles, before, dropped = trial, err, true
                break
            end
        end
        dropped || return poles, before
    end
end
function mergeclusters(S, ws, poles, iterations, before, acceptable)
    # the resolution of the data at a pole: the spacing of the samples
    # around its frequency, the first spacing below the band
    resolution(a) = begin
        k = searchsortedfirst(ws, abs(imag(a)))
        k <= 1 ? ws[2] - ws[1] : k > length(ws) ? ws[end] - ws[end - 1] : ws[k] - ws[k - 1]
    end
    # the representatives, real or in the upper half plane, a pole within
    # the resolution of the real axis being real
    reps = ComplexF64[]
    for a in poles
        imag(a) < -resolution(a) && continue
        push!(reps, abs(imag(a)) <= resolution(a) ? complex(real(a), 0.0) : a)
    end
    # the clusters, by the nearest representative within the resolution
    cluster = collect(1:length(reps))
    root(k) = (while cluster[k] != k; k = cluster[k]; end; k)
    for k in eachindex(reps), l in 1:k - 1
        abs(reps[k] - reps[l]) <= resolution(reps[k]) && (cluster[root(k)] = root(l))
    end
    merged = ComplexF64[]
    for r in unique(root.(eachindex(reps)))
        members = reps[root.(eachindex(reps)) .== r]
        a = sum(members)/length(members)
        if imag(a) == 0
            push!(merged, a)
        else
            push!(merged, a)
            push!(merged, conj(a))
        end
    end
    length(merged) == length(poles) && return poles, before
    merged = converge(S, ws, merged, iterations)
    err = fiterror(S, ws, merged)
    return acceptable(err, before) ? (merged, err) : (poles, before)
end

# the real basis of a pole set: columns of the matrix `Phi(s)` such that a
# real coefficient vector `c` gives `sum_p r_p/(s - a_p)` with conjugate
# residues on conjugate poles; and the map back from coefficients to the
# complex residues
function realbasis(ws::AbstractVector, poles::Vector{ComplexF64})
    K, N = length(ws), length(poles)
    Phi = zeros(ComplexF64, K, N)
    p = 1
    while p <= N
        a = poles[p]
        if imag(a) == 0 || abs(imag(a)) <= 1e-6*abs(a)
            Phi[:, p] .= 1 ./ (im .* ws .- real(a))
            p += 1
        else
            Phi[:, p] .= 1 ./ (im .* ws .- a) .+ 1 ./ (im .* ws .- conj(a))
            Phi[:, p + 1] .= im ./ (im .* ws .- a) .- im ./ (im .* ws .- conj(a))
            p += 2
        end
    end
    return Phi
end
function complexresidues(c::AbstractVector, poles::Vector{ComplexF64})
    N = length(poles)
    r = zeros(ComplexF64, N)
    p = 1
    while p <= N
        a = poles[p]
        if imag(a) == 0 || abs(imag(a)) <= 1e-6*abs(a)
            r[p] = c[p]
            p += 1
        else
            r[p] = complex(c[p], c[p + 1])
            r[p + 1] = conj(r[p])
            p += 2
        end
    end
    return r
end

# one relocation of the poles: the least squares of every entry with the
# shared weight, the weight's residues, and the zeros of the weight
function relocate(S::AbstractArray{<:Complex,3}, ws::AbstractVector, poles::Vector{ComplexF64})
    n, K, N = size(S, 1), length(ws), length(poles)
    Phi = realbasis(ws, poles)
    ne = n*n
    # The unknowns: per entry N residues and a constant, then the weight's
    # N residues and its constant; the equations: per entry K complex
    # samples with a zero right hand side, then the relaxation, all as
    # real equations. An entry's own unknowns enter only its own rows, so
    # they are eliminated entry by entry: the QR of the entry's block
    # `[A_e B_e]`, `A_e` the residue and constant columns, shared by every
    # entry, and `B_e` the weight's columns, leaves the trailing block of
    # its triangle as the entry's contribution to the weight's system,
    # `(N + 1)` rows against `(N + 1)` unknowns, and the weight is solved
    # from all of them and the relaxation at once. Nothing the size of
    # every entry's every sample is ever formed (Deschrijver's fast fit).
    K >= N + 1 || throw(ArgumentError(lazy"fitting $(N) poles needs at least $(N + 1) sample frequencies."))
    Ae = [real.(Phi) ones(K); imag.(Phi) zeros(K)]
    # the weight's columns scaled to their size over every entry
    scale = zeros(N + 1)
    for e in 1:ne
        i, j = (e - 1) % n + 1, (e - 1) ÷ n + 1
        Se = view(S, i, j, :)
        for p in 1:N
            scale[p] += sum(abs2, Se .* view(Phi, :, p))
        end
        scale[N + 1] += sum(abs2, Se)
    end
    scale = max.(sqrt.(scale), floatmin(Float64))
    reduced = zeros(ne*(N + 1) + 1, N + 1)
    rhs = zeros(ne*(N + 1) + 1)
    Be = zeros(2K, N + 1)
    for e in 1:ne
        i, j = (e - 1) % n + 1, (e - 1) ÷ n + 1
        Se = view(S, i, j, :)
        for p in 1:N
            v = Se .* view(Phi, :, p)
            Be[1:K, p] .= -real.(v) ./ scale[p]
            Be[K + 1:2K, p] .= -imag.(v) ./ scale[p]
        end
        Be[1:K, N + 1] .= -real.(Se) ./ scale[N + 1]
        Be[K + 1:2K, N + 1] .= -imag.(Se) ./ scale[N + 1]
        R = qr(hcat(Ae, Be)).R
        reduced[(e - 1)*(N + 1) + 1:e*(N + 1), :] .= R[N + 2:2N + 2, N + 2:2N + 2]
    end
    # the relaxation: the real part of the weight averages to one; scaled
    # to the size of the data so that it weighs as one sample of it
    weight = sqrt(sum(abs2, S)/K)
    for p in 1:N
        reduced[end, p] = weight*sum(real, view(Phi, :, p))/K/scale[p]
    end
    reduced[end, N + 1] = weight/scale[N + 1]
    rhs[end] = weight*1.0
    x = (reduced \ rhs) ./ scale
    csigma = x[1:N]
    dsigma = x[N + 1]
    # a weight whose constant vanished has its zeros at its poles: it is
    # held away from zero, as Gustavsen does
    abs(dsigma) < 1e-8 && (dsigma = copysign(1e-8, dsigma))
    # the zeros of the weight: the eigenvalues of A - b r'/d in the real
    # form of the pole set
    Ar, br = realpolematrix(poles)
    cr = realresiduerow(csigma, poles)
    zeros_ = eigvals(Ar .- br*transpose(cr) ./ dsigma)
    newpoles = ComplexF64[]
    for z in zeros_
        real(z) > 0 && (z = complex(-real(z), imag(z)))
        push!(newpoles, z)
    end
    # conjugate pairs paired, real poles real
    sort!(newpoles; by = x -> (real(x), imag(x)))
    out = ComplexF64[]
    used = falses(length(newpoles))
    for (k, z) in enumerate(newpoles)
        used[k] && continue
        if abs(imag(z)) <= 1e-6*abs(z)
            push!(out, complex(real(z), 0.0))
            used[k] = true
        else
            partner = findfirst(l -> !used[l] && l != k && abs(newpoles[l] - conj(z)) < 1e-8*abs(z), eachindex(newpoles))
            if isnothing(partner)
                push!(out, complex(real(z), 0.0))
                used[k] = true
            else
                push!(out, complex(real(z), abs(imag(z))))
                push!(out, complex(real(z), -abs(imag(z))))
                used[k] = true
                used[partner] = true
            end
        end
    end
    return out
end

# the real block diagonal matrix of a pole set, its input column of ones
# in the real form, and the coefficient row of the weight's residues in
# the same form, so that `A - b c'` has the weight's zeros as eigenvalues
function realpolematrix(poles::Vector{ComplexF64})
    N = length(poles)
    A = zeros(N, N)
    b = zeros(N)
    p = 1
    while p <= N
        a = poles[p]
        if imag(a) == 0 || abs(imag(a)) <= 1e-6*abs(a)
            A[p, p] = real(a)
            b[p] = 1.0
            p += 1
        else
            A[p, p] = real(a); A[p, p + 1] = imag(a)
            A[p + 1, p] = -imag(a); A[p + 1, p + 1] = real(a)
            b[p] = 2.0; b[p + 1] = 0.0
            p += 2
        end
    end
    return A, b
end
function realresiduerow(c::AbstractVector, poles::Vector{ComplexF64})
    N = length(poles)
    row = zeros(N)
    p = 1
    while p <= N
        a = poles[p]
        if imag(a) == 0 || abs(imag(a)) <= 1e-6*abs(a)
            row[p] = c[p]
            p += 1
        else
            row[p] = c[p]; row[p + 1] = c[p + 1]
            p += 2
        end
    end
    return row
end

# the residues of every entry and the constant at fixed poles, by least
# squares in the real basis
function fitresidues(S::AbstractArray{<:Complex,3}, ws::AbstractVector, poles::Vector{ComplexF64})
    n, K, N = size(S, 1), length(ws), length(poles)
    Phi = realbasis(ws, poles)
    M = zeros(2K, N + 1)
    M[1:K, 1:N] .= real.(Phi)
    M[K + 1:2K, 1:N] .= imag.(Phi)
    M[1:K, N + 1] .= 1.0
    residues = zeros(ComplexF64, n, n, N)
    D = zeros(n, n)
    F = qr(M)
    for i in 1:n, j in 1:n
        x = F \ vcat(real.(view(S, i, j, :)), imag.(view(S, i, j, :)))
        residues[i, j, :] .= complexresidues(x[1:N], poles)
        D[i, j] = x[N + 1]
    end
    return residues, D
end

# the real state space realization of poles and residue matrices: every
# real pole `a` with residue `R` is `n` states `A = a I, B = I, C = R`, and
# every conjugate pair `a, conj(a)` with residues `R, conj(R)` the `2n`
# states `A = [Re a I  -Im a I; Im a I  Re a I]`, `B = [I; 0]`,
# `C = [2 Re R  -2 Im R]`
function realization(poles::Vector{ComplexF64}, residues::AbstractArray{<:Complex,3}, n::Int)
    N = length(poles)
    blocks = Tuple{Matrix{Float64},Matrix{Float64},Matrix{Float64}}[]
    p = 1
    while p <= N
        a = poles[p]
        R = residues[:, :, p]
        # the states of a pole are the rank of its residue, `R = U S V'`
        # split as `C = U sqrt(S)`, `B = sqrt(S) V'`, so that the
        # realization is minimal and balanced; a real pole with a
        # residue of rank `r` takes `r` states, a pair `2r`
        F = svd(R)
        r = count(s -> s > 1e-12*max(F.S[1], floatmin(Float64)), F.S)
        Bc = Diagonal(sqrt.(F.S[1:r]))*F.Vt[1:r, :]
        Cc = F.U[:, 1:r]*Diagonal(sqrt.(F.S[1:r]))
        if imag(a) == 0 || abs(imag(a)) <= 1e-6*abs(a)
            r > 0 && push!(blocks, (real(a) .* Matrix(1.0I, r, r), real.(Bc), real.(Cc)))
            p += 1
        else
            # the complex state `z' = a z + Bc u`, `y = 2 Re(Cc z)`, as its
            # real and imaginary parts
            α, β = real(a), imag(a)
            r > 0 && push!(blocks, ([α .* Matrix(1.0I, r, r)  -β .* Matrix(1.0I, r, r); β .* Matrix(1.0I, r, r)  α .* Matrix(1.0I, r, r)],
                [real.(Bc); imag.(Bc)], [2 .* real.(Cc)  -2 .* imag.(Cc)]))
            p += 2
        end
    end
    isempty(blocks) && return zeros(0, 0), zeros(0, n), zeros(n, 0)
    nz = sum(size(b[1], 1) for b in blocks)
    A, B, C = zeros(nz, nz), zeros(nz, n), zeros(n, nz)
    offset = 0
    for (Ab, Bb, Cb) in blocks
        k = size(Ab, 1)
        A[offset + 1:offset + k, offset + 1:offset + k] .= Ab
        B[offset + 1:offset + k, :] .= Bb
        C[:, offset + 1:offset + k] .= Cb
        offset += k
    end
    return A, B, C
end

# The passivity enforcement. Where the largest singular value of the fit
# crosses one, found from the imaginary eigenvalues of the Hamiltonian
# matrix, which bound the violating bands, and the worst point of each
# band, the residues and the constant are perturbed by the least change,
# in the norm of the fit over the samples, that brings the singular value
# at each worst point to one less a margin, a linear constraint on the
# perturbation through the singular vectors; repeated until the test
# passes, with a bound on the rounds.
function enforcepassivity(A, B, C, D, ws; atol = 1e-8, rounds::Int = 20)
    n, nz = size(D, 1), size(A, 1)
    for round in 1:rounds
        bands = violations(A, B, C, D; atol = atol)
        isempty(bands) && return A, B, C, D
        # the worst point of each band
        points = Float64[]
        for (w1, w2) in bands
            grid = range(w1, w2; length = 51)
            k = argmax([opnorm(D .+ C*((im*w*I - A) \ B)) for w in grid])
            push!(points, grid[k])
        end
        # the perturbation of C and D, `delta S(i w) = delta C (i w I - A)^(-1) B + delta D`,
        # as a real vector; its effect at the samples for the norm and at
        # the points for the constraints
        nunk = n*nz + n*n
        function perturbation(w)
            X = (im*w*I - A) \ B
            # delta S[i, j] = sum_k delta C[i, k] X[k, j] + delta D[i, j]
            M = zeros(ComplexF64, n*n, nunk)
            for i in 1:n, j in 1:n
                row = (i - 1)*n + j
                for k in 1:nz
                    M[row, (k - 1)*n + i] = X[k, j]
                end
                M[row, n*nz + (i - 1)*n + j] = 1.0
            end
            return M
        end
        # the norm: the fit over the samples, as the normal matrix
        # accumulated sample by sample, never the stacked samples
        H = zeros(nunk, nunk)
        for w in ws
            M = perturbation(w)
            H .+= real.(M'*M)
        end
        # the constraints: Re(u' delta S v) = 1 - margin - sigma at each point
        margin = 1e-6
        Ceq = zeros(length(points), nunk)
        ceq = zeros(length(points))
        for (q, w) in enumerate(points)
            Sw = D .+ C*((im*w*I - A) \ B)
            F = svd(Sw)
            u, v, σ = F.U[:, 1], F.V[:, 1], F.S[1]
            M = perturbation(w)
            # u' delta S v = sum_ij conj(u_i) deltaS_ij v_j
            row = zeros(ComplexF64, nunk)
            for i in 1:n, j in 1:n
                row .+= conj(u[i])*v[j] .* view(M, (i - 1)*n + j, :)
            end
            Ceq[q, :] .= real.(row)
            ceq[q] = 1 - margin - σ
        end
        # the least norm change under the constraints: the KKT system
        H .+= 1e-12*I
        KKT = [H transpose(Ceq); Ceq zeros(length(points), length(points))]
        sol = KKT \ vcat(zeros(nunk), ceq)
        δ = sol[1:nunk]
        # the unknowns were laid out as delta C[i, k] at (k - 1) n + i and
        # delta D[i, j] at (i - 1) n + j
        C = C .+ reshape(δ[1:n*nz], n, nz)
        D = D .+ transpose(reshape(δ[n*nz + 1:end], n, n))
    end
    worst, _ = hinfnorm(A, B, C, D)
    throw(ArgumentError(lazy"the fit could not be made passive in $(rounds) rounds: its largest singular value is still $(worst); a fit this far from passive is a poor fit, so fit with more poles, or over a narrower band."))
end

# the bands of angular frequency where the largest singular value exceeds
# one: between the crossings of one, found by the pencil of the passivity
# test, the largest singular value is tested at a point of each interval,
# and adjacent intervals above one are one band; where the pencil is
# singular the band of the poles is sampled instead
function violations(A, B, C, D; atol = 1e-8)
    # the crossings of the level `1 + atol/2`, so that a fit lossless to
    # within the tolerance has none and the pencil is regular
    level = 1 + atol/2
    An, Bn, Cn, wscale = balancedrealization(A, B, C)
    crossings = pencilcrossings(An, Bn, Cn ./ level, D ./ level; singulartest = false)[1] .* wscale
    worst, w = hinfnorm(A, B, C, D)
    isempty(crossings) && return worst > 1 + atol ? [(max(w - wscale/100, 0.0), w + wscale/100)] : Tuple{Float64,Float64}[]
    edges = vcat(0.0, crossings, 2last(crossings) + wscale)
    bands = Tuple{Float64,Float64}[]
    for k in 1:length(edges) - 1
        w = (edges[k] + edges[k + 1])/2
        if opnorm(D .+ C*((im*w*I - A) \ B)) > level
            if !isempty(bands) && last(bands)[2] == edges[k]
                bands[end] = (last(bands)[1], edges[k + 1])
            else
                push!(bands, (edges[k], edges[k + 1]))
            end
        end
    end
    return bands
end
