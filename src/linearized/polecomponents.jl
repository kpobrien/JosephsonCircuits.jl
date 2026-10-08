# One analytic contribution L*f(scale*z + i*offset)*(R0 + z*R1).
# The frequency argument belongs to the INPUT harmonic of a converting
# device. Its output rows may belong to another harmonic.
# Matrices use compact coordinates; rows/cols embed the local contribution
# in the circuit. No term stores circuit-width sparse column pointers.
struct PoleResponseTerm
    response::Any
    offset::Float64
    rows::Vector{Int}
    cols::Vector{Int}
    left::Matrix{ComplexF64}
    right0::Matrix{ComplexF64}
    right1::Matrix{ComplexF64}
end

function poleterm(t::PoleResponseTerm, z, scale)
    value = t.response(scale*z + im*t.offset)
    all(isfinite, value) || throw(ArgumentError("a Laplace response is nonfinite at s = $(scale*z + im*t.offset)."))
    size(value) == (size(t.left,2), size(t.right0,1)) ||
        throw(DimensionMismatch("the Laplace callback returned a matrix of the wrong size."))
    return t.left*ComplexF64.(value)*(t.right0 + z*t.right1)
end

# Gather only the coordinates touched by a component. The transposed CSC
# maps give direct access to each physical node's entries without scanning
# the full circuit width for every response. Incidences contain one list
# of (physical node row, coefficient) pairs per port.
function polelocalmaps(fluxrows, commonrows, incidences, extracols = Int[])
    cols = copy(extracols)
    for entries in incidences, (row, _) in entries, map in (fluxrows, commonrows)
        append!(cols, view(rowvals(map), nzrange(map, row)))
    end
    sort!(unique!(cols))
    index = Dict(col => j for (j, col) in enumerate(cols))
    F = zeros(ComplexF64, length(incidences), length(cols))
    V = similar(F)
    fill!(V, 0)
    for (q, entries) in enumerate(incidences), (row, coefficient) in entries
        for (map, localmap) in ((fluxrows, F), (commonrows, V)), k in nzrange(map, row)
            localmap[q, index[rowvals(map)[k]]] += coefficient*nonzeros(map)[k]
        end
    end
    return cols, F, V
end

# A rational filter is kept in state space, including unobservable states.
# Couplings are (output harmonic minus input harmonic, output coefficient).
struct PoleFilter
    provider::Any
    couplings::Vector{Tuple{Int,ComplexF64}}
    statebase::Int
end

struct PoleBlock
    definition::Any
    signal::Vector{Int}
    ref::Vector{Int}
    currentbase::Int
    filters::Vector{PoleFilter}
end

polenstates(p::RationalScatteringProvider) = size(p.A, 1)
polenstates(p) = 0

polefrequency(p::RationalScatteringProvider) = maximum(abs, [p.groups.poles; p.states.poles]; init = 0.0)
polefrequency(p::ModulatedRationalProvider) = max(polefrequency(p.cosine), polefrequency(p.sine))
polefrequency(p::TransmissionLineProvider) = p.delay > 0 ? inv(p.delay) : 0.0
polefrequency(p) = 0.0
polefrequency(d::ScatteringParameters) = polefrequency(d.provider)
polefrequency(d::LinearizedScattering) = maximum(polefrequency, d.providers; init = 0.0)

# A diagonal state similarity balances B/scale against C. It preserves
# even hidden states, unlike a minimal transfer-function realization, and
# avoids a GHz pole represented by B=1,C=1e10 losing accuracy in the pencil.
function polebalanced(p::RationalScatteringProvider, scale)
    t = statebalance(p.B, p.C, scale)
    return (p.A ./ t).*transpose(t), p.B ./ t, p.C .* transpose(t)
end

# the factors `t` of that similarity, the states `t .* z` of the balanced
# realization's `z`, for the inputs `B` and the outputs `C`
function statebalance(B, C, scale)
    t = Float64[]
    for k in axes(B,1)
        b, c = norm(view(B,k,:))/scale, norm(view(C,:,k))
        push!(t, b > 0 && c > 0 ? sqrt(b)/sqrt(c) : c > 0 ? inv(c) : b > 0 ? b : 1.0)
    end
    all(x -> isfinite(x) && x > 0, t) || throw(ArgumentError("the rational state scaling is nonfinite; rescale the supplied realization."))
    return t
end

polezero(p::ConstantMatrixProvider) = p.A
polezero(p::RationalScatteringProvider) = p.D - p.C*(p.A\p.B)
polezero(p::TransmissionLineProvider) = [0.0 1.0; 1.0 0.0]
polezero(p::CallableMatrixProvider) = poleproviderresponse(p, 0.0+0im)

function poleproviderresponse(p::CallableMatrixProvider, s)
    p.form == :matrix && return p.f.f(s)
    value = zeros(ComplexF64, p.n, p.n)
    p.f.f(value, s)
    return value
end

function poleproviderresponse(p::TransmissionLineProvider, s)
    through = exp(-s*p.delay)
    return [0 through; through 0]
end

function checkpoledc(d, p, path)
    d.dcmodel isa ScatteringLimit && return nothing
    stated = dcscatteringmatrix(d.dcmodel, d.nports)
    actual = polezero(p)
    isapprox(stated, actual; atol = d.atol, rtol = d.atol) ||
        throw(ArgumentError("$path has a DC override different from its analytic zero-frequency limit; supply a consistent causal model for pole analysis."))
    return nothing
end

function checkpoleprovider(p, path; realresponse = true)
    if p isa ConstantMatrixProvider
        all(isfinite, p.A) && (!realresponse || all(isreal, p.A)) ||
            throw(ArgumentError("$path requires a finite real constant response or an explicit causal model."))
    elseif p isa CallableMatrixProvider
        p.f isa LaplaceResponse || throw(ArgumentError(
            "$path needs a LaplaceResponse or a RationalScattering fit for pole analysis."))
    elseif !(p isa Union{RationalScatteringProvider,TransmissionLineProvider})
        throw(ArgumentError("$path has frequency data without a Laplace realization; fit it with RationalScattering."))
    end
    return nothing
end

function poleblocks(psc, omega, m, base)
    blocks = PoleBlock[]
    for cb in psc.scatteringblocks
        d = cb.definition
        currentbase = base
        base += d.nports*m
        filters = PoleFilter[]
        if d isa LinearizedScattering
            # The envelope is a transient startup option; HB and this
            # analysis describe the fully-on periodic operating point.
            ratio = omega > 0 ? d.wp/omega : 0.0
            # a harmonic of the analysis's pump as harmonic balance tells one
            (!any(>(0), d.harmonics) || (omega > 0 && isharmonic(d.wp, round(Int, ratio), omega))) ||
                throw(ArgumentError("$(cb.path) needs a pumpfrequency commensurate with its pump $(d.wp)."))
            for (k, p) in zip(d.harmonics, d.providers)
                h = k*round(Int, ratio)
                rot = cis(k*d.phase)
                if p isa ModulatedRationalProvider
                    k > 0 || throw(ArgumentError("H_0 must be an unmodulated response."))
                    for (part, weight) in ((p.cosine, rot), (p.sine, im*rot))
                        push!(filters, PoleFilter(part, [(h, weight), (-h, conj(weight))], base))
                        base += polenstates(part)*m
                    end
                else
                    checkpoleprovider(p, cb.path; realresponse = k == 0)
                    k == 0 && checkpoledc(d, p, cb.path)
                    # Native real causal providers are their own analytic
                    # reflection, including a line's exact delay response.
                    if k == 0 || p isa Union{RationalScatteringProvider,TransmissionLineProvider}
                        pairs = k == 0 ? [(0, 1.0+0im)] : [(h, rot), (-h, conj(rot))]
                        push!(filters, PoleFilter(p, pairs, base))
                        base += polenstates(p)*m
                    else
                        # Complex converted responses need the analytic
                        # reflection conj(H(conj(s))), not conj(H(s)).
                        push!(filters, PoleFilter(p, [(h, rot)], base))
                        reflected = p isa ConstantMatrixProvider ? ConstantMatrixProvider(conj.(p.A)) :
                            CallableMatrixProvider(LaplaceResponse(s -> conj.(poleproviderresponse(p, conj(s)))), d.nports)
                        push!(filters, PoleFilter(reflected, [(-h, conj(rot))], base))
                    end
                end
            end
        else
            checkpoleprovider(d.provider, cb.path)
            checkpoledc(d, d.provider, cb.path)
            push!(filters, PoleFilter(d.provider, [(0, 1.0+0im)], base))
            base += polenstates(d.provider)*m
        end
        push!(blocks, PoleBlock(d, cb.signalnodes .- 1, cb.refnodes .- 1,
            currentbase, filters))
    end
    return blocks, base
end

# Native (unscaled-frequency) polynomial stamps, added to K and G, which
# are returned. Columns are node flux, Lscale*port current and
# Lscale*filter state, with harmonic fastest. The stamps are gathered as
# triplets and summed into both at once, since an entry inserted into a
# sparse matrix moves every entry after it.
function stampoleblocks(K, G, blocks, modes, m, Lscale, scale)
    Kstamps = (Int[], Int[], ComplexF64[])
    Gstamps = (Int[], Int[], ComplexF64[])
    modeindex = Dict(only(mode) => i for (i, mode) in enumerate(modes))
    for b in blocks
        d, base = b.definition, b.currentbase
        for q in 1:d.nports, a in 1:m
            r = base + (q-1)*m + a
            root = sqrt(d.zref[q])
            stamp!(Kstamps, r, r, -root)
            for (node, sign) in ((b.signal[q], 1), (b.ref[q], -1))
                node == 0 && continue
                col = (node-1)*m + a
                stamp!(Kstamps, col, r, sign)
                stamp!(Gstamps, r, col, sign*Lscale/root)
            end
        end
        for filter in b.filters
            p, zb = filter.provider, filter.statebase
            rational = p isa RationalScatteringProvider
            p isa Union{ConstantMatrixProvider,RationalScatteringProvider} || continue
            D = rational ? p.D : p.A
            nz = polenstates(p)
            if rational
                Ab, Bb, Cb = polebalanced(p, scale)
                for a in 1:m, j in 1:nz
                    r = zb + (j-1)*m + a
                    stamp!(Gstamps, r, r, 1)
                    for k in 1:nz
                        stamp!(Kstamps, r, zb + (k-1)*m + a, -Ab[j,k])
                    end
                    for q in 1:d.nports
                        root = sqrt(d.zref[q])
                        stamp!(Kstamps, r, base + (q-1)*m + a, -(Bb[j,q]*root/2))
                        for (node, sign) in ((b.signal[q], 1), (b.ref[q], -1))
                            node == 0 && continue
                            stamp!(Gstamps, r, (node-1)*m + a, -(sign*Lscale*Bb[j,q]/(2root)))
                        end
                    end
                end
            end
            for (h, weight) in filter.couplings, a in 1:m
                out = get(modeindex, only(modes[a])+h, 0)
                out == 0 && continue
                for rport in 1:d.nports
                    r = base + (rport-1)*m + out
                    for q in 1:d.nports
                        root = sqrt(d.zref[q])
                        stamp!(Kstamps, r, base + (q-1)*m + a, -(weight*D[rport,q]*root))
                        for (node, sign) in ((b.signal[q], 1), (b.ref[q], -1))
                            node == 0 && continue
                            stamp!(Gstamps, r, (node-1)*m + a, -(sign*weight*D[rport,q]*Lscale/root))
                        end
                    end
                    if rational
                        for k in 1:nz
                            stamp!(Kstamps, r, zb + (k-1)*m + a, -(2weight*Cb[rport,k]))
                        end
                    end
                end
            end
        end
    end
    # a circuit without blocks has no stamps and keeps its matrices
    n = size(K, 1)
    Ks = isempty(first(Kstamps)) ? K : K + sparse(Kstamps..., n, n)
    Gs = isempty(first(Gstamps)) ? G : G + sparse(Gstamps..., n, n)
    return Ks, Gs
end

# one term of a sum of stamps held as triplets
function stamp!((rows, cols, vals), i, j, v)
    push!(rows, i); push!(cols, j); push!(vals, v)
    return nothing
end

function poleblockterms(blocks, modes, offsets, fluxrows, commonrows, scale, Lscale)
    terms = PoleResponseTerm[]
    m = length(modes)
    modeindex = Dict(only(mode) => i for (i, mode) in enumerate(modes))
    for b in blocks, filter in b.filters
        p = filter.provider
        p isa Union{ConstantMatrixProvider,RationalScatteringProvider} && continue
        d = b.definition
        response = s -> poleproviderresponse(p, s)
        for a in 1:m
            incidences = [[((node-1)*m+a, sign/sqrt(d.zref[q]))
                for (node, sign) in ((b.signal[q], 1), (b.ref[q], -1)) if node > 0]
                for q in 1:d.nports]
            currents = [b.currentbase + (q-1)*m+a for q in 1:d.nports]
            cols, F, V = polelocalmaps(fluxrows, commonrows, incidences, currents)
            R0 = Lscale*(scale*V + im*offsets[a]*F)
            for q in 1:d.nports
                R0[q, searchsortedfirst(cols, currents[q])] += sqrt(d.zref[q])
            end
            R1 = Lscale*scale*F
            for (h, weight) in filter.couplings
                out = get(modeindex, only(modes[a])+h, 0)
                out == 0 && continue
                rows = [b.currentbase + (q-1)*m+out for q in 1:d.nports]
                left = Matrix(Diagonal(fill(-weight, d.nports)))
                push!(terms, PoleResponseTerm(response, offsets[a], rows, cols, left, R0, R1))
            end
        end
    end
    return terms
end

function polelumpedvalues(psc, values)
    static = Any[values...]
    dynamic = Int[]
    for (k, value) in enumerate(values)
        kind = psc.componenttypes[k]
        # Independent current sources have zero perturbation, whatever the
        # frequency dependence of their specified drive.
        if kind == :I
            static[k] = 0.0
        elseif value isa Number && isreal(value) &&
                (isfinite(value) || (kind == :R && value == Inf))
            continue
        elseif kind in (:R, :C, :L) && poleanalytic(value) && !(value isa Number)
            if kind == :L
                dc = laplacevalue(value, 0.0+0im)
                isfinite(dc) && isreal(dc) && !iszero(dc) || throw(ArgumentError(
                    "$(psc.componentnames[k]) needs finite nonzero DC inductance, as in harmonic balance; use a state realization for other constitutive laws."))
            end
            push!(dynamic, k)
            static[k] = kind == :R ? Inf : kind == :C ? 0.0 : real(laplacevalue(value, 0.0+0im))
        else
            throw(ArgumentError("$(psc.componentnames[k]) needs a real constant or an explicit LaplaceResponse for pole analysis."))
        end
    end
    # HB itself requires mutually coupled inductances to be constants.
    for (_, l1, l2) in psc.couplings
        (l1 in dynamic || l2 in dynamic) && throw(ArgumentError("mutually coupled inductances must be frequency independent."))
    end
    return static, dynamic
end

function polelumpedterms(psc, values, dynamic, m, offsets, fluxrows, commonrows, scale, Lscale)
    terms = PoleResponseTerm[]
    for k in dynamic, a in 1:m
        value, kind = values[k], psc.componenttypes[k]
        incidence = Dict{Int,Float64}()
        for (node, sign) in ((psc.nodeindices[1,k]-1, 1), (psc.nodeindices[2,k]-1, -1))
            if node > 0
                row = (node-1)*m+a
                incidence[row] = get(incidence, row, 0.0) + sign
            end
        end
        rows = sort!([row for (row, coefficient) in incidence if !iszero(coefficient)])
        isempty(rows) && continue
        entries = [(row, incidence[row]) for row in rows]
        left = reshape(ComplexF64[Lscale*incidence[row] for row in rows], :, 1)
        cols, F, V = polelocalmaps(fluxrows, commonrows, [entries])
        if kind == :L
            # Replace the DC inductance used for topology and scaling with
            # the actual Laplace law, without changing the static graph.
            invdc = inv(laplacevalue(value, 0.0+0im))
            response = s -> fill(inv(laplacevalue(value, s))-invdc, 1, 1)
            R0, R1 = F, zeros(ComplexF64, size(F))
        else
            response = kind == :R ? (s -> fill(inv(laplacevalue(value, s)), 1, 1)) :
                (s -> fill(s*laplacevalue(value, s), 1, 1))
            R0, R1 = scale*V + im*offsets[a]*F, scale*F
        end
        push!(terms, PoleResponseTerm(response, offsets[a], rows, cols, left, R0, R1))
    end
    return terms
end
