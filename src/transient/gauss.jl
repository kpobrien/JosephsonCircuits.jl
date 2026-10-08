# The two stage Gauss-Legendre collocation on the second order system: two
# stage fluxes per step, solved together by the shared Newton engine on the
# true stage equations, with the correction from one complex factorization.
# The stage matrix of the tableau has the conjugate eigenvalues
# `mu = 3 ± i sqrt(3)`, so with one frozen junction stiffness the two
# stages decouple into the complex system `(mu/h)^2 C + (mu/h) G + L + J*`
# and its conjugate, of which one is solved; that is a simplified Newton,
# converging linearly on the true residual. The complex matrix is the
# harmonic balance real Jacobian plan's assembly of the real part on its
# pattern, plus a constant imaginary part on the same pattern. The tangent
# and the adjoint, which are linear, solve the linearized stage equations
# exactly instead: the real matrix of both stages, each at its own junction
# stiffness, factorized at each step (see `StagePlan`).

"""
    GaussLegendre()

The two stage Gauss-Legendre collocation on the flux and its rate, the
default of [`transientsolve`](@ref): fourth order, A-stable, symplectic,
and free of numerical damping, so a lossless LC oscillation keeps its
energy and a resonator's angular frequency `w` is warped by
`(w dt)^4/720` rather than the trapezoidal rule's `(w dt)^2/12`. Each
step solves the two stage equations together, by Newton on one
complex factorization of the stage matrix at a frozen junction
stiffness, refreshed as the trapezoidal rule's is. Not L-stable: an
unresolved fast mode is not damped, as the trapezoidal rule does not damp
it.

Along a direction of the state without capacitance the equations are
algebraic. Where a resistor acts along it the rate is what the
constraint determines, read from it at the endpoint, and converges at
fourth order as the rest of the state does. Where none does the
constraint is on the flux alone: a junction on the direction or a source
driving it makes it nonlinear or moving, and each step projects its
endpoint onto it, while a linear constraint no source drives is an
invariant of the rule and holds by itself. The rate along every such
direction is not the rule's, whose update carries the rounding of every
stage solve forward along it, but is read from the differentiated
constraint wherever the state is reported, as the trapezoidal rule's is;
where a scattering block is on a direction, whose states move the
constraint at a rate the reading would have to solve for with the
block's port currents, the rate along the projected directions is the
derivative of the cubic through the state, the two stages and the
projected endpoint, third order. The tangent and the adjoint of a
Gauss-Legendre solve differentiate the full stage equations, the
projection and the readings included, and are exact for the recorded
steps: each step's linearized stage equations, the real matrix of both
stages at the step's recorded phases with the blocks' coupling at the
stages' weights, are factorized at the step and solved once for every
direction, the adjoint's by the transpose of the same matrix.
"""
struct GaussLegendre <: AbstractTransientIntegrator end

# The tableau, and what a step reads of it: the abscissae `c`, the
# tableau `A`, which the rational blocks' stages read, its inverse
# `A^{-1}` and the square of that as column major tuples, `A^{-1} 1`, the endpoint
# weights `b' A^{-1}` of the flux and `b' A^{-2}` of the rate, the
# predictor of the next stages from the last increments, and the
# eigenvalue `mu` of `A^{-1}` in the upper half plane with the entries of
# its eigenvector matrix `T = [t conj(t)]` and of the inverse the stage
# transform reads, and the derivative at the endpoint of the cubic
# through the state, the two stages and the endpoint, as weights on the
# four, for the rate along an algebraic direction a block is on.
struct GaussCoefficients
    c::NTuple{2,Float64}
    a::NTuple{4,Float64}
    ainv::NTuple{4,Float64}
    ainv2::NTuple{4,Float64}
    ainvone::NTuple{2,Float64}
    ex::NTuple{2,Float64}
    ev::NTuple{2,Float64}
    predict::NTuple{4,Float64}
    endrate::NTuple{4,Float64}
    mu::ComplexF64
    t11::ComplexF64
    t21::ComplexF64
    tinv11::ComplexF64
    tinv12::ComplexF64
end

function gausscoefficients()
    s3 = sqrt(3.0)
    A = [1/4 1/4-s3/6; 1/4+s3/6 1/4]
    c = (1/2 - s3/6, 1/2 + s3/6)
    b = [1/2, 1/2]
    Ainv = inv(A)
    Ainv2 = Ainv*Ainv
    mu = 3 + im*s3
    t = [-3 + 2s3, im*s3]
    T = [t conj(t)]
    Tinv = inv(T)
    # the eigen decomposition of A^{-1} the transform relies on
    isapprox(Ainv*t, mu .* t; rtol = 1e-12) || error("the Gauss-Legendre transform is inconsistent.")
    ex = transpose(b)*Ainv
    ev = transpose(b)*Ainv2
    # the start of the next step's stages: the collocation polynomial of
    # the step just taken, the quadratic through the state and the two
    # stage increments, extrapolated to the new stage times and taken
    # relative to the new state, as fixed coefficients on the increments
    M = [c[1] c[1]^2; c[2] c[2]^2]
    P = [c[1] 2c[1]+c[1]^2; c[2] 2c[2]+c[2]^2]*inv(M)
    # the derivative at 1 of the Lagrange basis on the nodes 0, c1, c2, 1
    nodes = (0.0, c[1], c[2], 1.0)
    endrate = ntuple(4) do j
        sum(prod((1 - nodes[l])/(nodes[j] - nodes[l]) for l in 1:4 if l != j && l != m; init = 1.0)/(nodes[j] - nodes[m])
            for m in 1:4 if m != j)
    end
    isapprox(sum(endrate[j]*nodes[j]^3 for j in 1:4), 3.0; atol = 1e-12) || error("the endpoint derivative is inconsistent.")
    return GaussCoefficients(c, (A[1, 1], A[2, 1], A[1, 2], A[2, 2]), (Ainv[1, 1], Ainv[2, 1], Ainv[1, 2], Ainv[2, 2]),
        (Ainv2[1, 1], Ainv2[2, 1], Ainv2[1, 2], Ainv2[2, 2]), (sum(Ainv[1, :]), sum(Ainv[2, :])),
        (ex[1], ex[2]), (ev[1], ev[2]), (P[1, 1], P[2, 1], P[1, 2], P[2, 2]), endrate, mu, T[1, 1], T[2, 1], Tinv[1, 1], Tinv[1, 2])
end

# the entries of a column major 2 by 2 tuple
@inline entry(m::NTuple{4,Float64}, i, j) = m[i + 2(j - 1)]

# The stage algebra of a rational block. Its states at the two stages are
# linear in the state at the start and the incident waves at the stages,
# `Z = Zz z + Zu U` with `Z = M^(-1) (1 (x) I)` and `Zu = M^(-1) h (A_tab (x) B)`
# for `M = I - h A_tab (x) A`, the reflected waves are `b = C Z + D U`,
# and the state at the end is `z' = Ez z + Eu U` with the Gauss weights;
# all constant matrices, factorized once per step size. In the complex
# stage basis the same algebra is the block's scattering matrix at the
# stage frequency `mu/h`, which the frozen operator carries exactly.
struct RationalStage
    block::Int
    Zz::Matrix{Float64}
    Zu::Matrix{Float64}
    Ez::Matrix{Float64}
    Eu::Matrix{Float64}
end


# The coupling of every rational block's stage algebra, grouped: sparse
# operators on the stage stacked unknowns `[stage 1; stage 2]` of a
# batch, so that the reflected waves of all the blocks at both stages,
# and their states after the step, are a few products on the backend
# with no per-block or per-condition work on the host. With `G` the
# gather of the rates across the ports and the port currents, `S` the
# scatter onto the blocks' rows, `W_d`, `W_x` the incident waves from
# the gathered increments and values, `Z_u`, `Z_z` the stacked states at
# the stages from the incident waves and the states, and `E_z`, `E_u`
# the states' update: `P_d = Z_u W_d G`, `P_x = Z_u W_x G`, `P_z = Z_z`
# give the stacked states, the reflected waves at a stage are
# `sum_j w_j(t) C_j` times the states of that stage, `C_0` the output
# matrix of every block's unconverted response with the weight one and
# the rest the modulated outputs of the pumped blocks with the weights
# of the stage's time (see [`BlockModulation`](@ref)), scattered by `S`;
# `E_d = E_u W_d G`, `E_x = E_u W_x G`; the transposes the adjoint reads;
# and on the host the output matrices `C_j` for the resting waves. A circuit
# without a pumped block has the one term `S C_0`, which is the coupling
# of a block which does not convert.
struct RationalCoupling{M}
    nstates::Int
    Pd::M
    Px::M
    Pz::M
    Ez::M
    Ed::M
    Ex::M
    Pxt::M
    Pzt::M
    Ezt::M
    Edt::M
    Ext::M
    # the scatter of the output of each term onto the blocks' rows,
    # `(n, nstates)`, and its transpose
    SC::Vector{M}
    SCt::Vector{M}
    # which block and which modulated output each term after the first
    # is, for its weight at a time
    terms::Vector{Tuple{Int,Int}}
    # the output matrices of every term on the host, `(nports, nstates)`,
    # for the resting waves
    Cblk::Vector{SparseMatrixCSC{Float64,Int}}
    # on the host, the scatter onto the blocks' rows and the stacked
    # states from the stage unknowns, `P_d + P_x`, for the exact stage
    # solve of a pumped block (see [`StageCorrection`](@ref)) and the two
    # stages' matrix (see `StagePlan`)
    Shost::SparseMatrixCSC{Float64,Int}
    Phost::SparseMatrixCSC{Float64,Int}
    # the port rows with a modulated output term, the only rows the
    # stage correction acts on, and the modulated outputs of every term
    # after the first at each stage on those rows from the stacked stage
    # unknowns, `W_ij = (C_j P_i)[modulated, :]` with `P_i` the stage's
    # rows of `P_d + P_x`, stage by stage and term by term, on the backend
    modulated::Vector{Int}
    W::Vector{M}
    # the output scatter entrywise in magnitude, which bounds its rounding
    # against the states themselves rather than against their largest:
    # the scatter and the states both run over the decades the poles of a
    # fit span, and a norm of each pairs the largest of one with the
    # largest of the other whatever rows they are in; and the largest row
    # sum of each, which gives the cheap bound that decides whether the
    # entrywise one is needed
    SCabs::Vector{M}
    SCrows::Vector{Float64}
end

function rationalcoupling(p::TransientProblem, stages::Vector{RationalStage}, gc::GaussCoefficients, h, Lscale,
        blockgather::SparseMatrixCSC, blockscatter::SparseMatrixCSC, backend)
    nports = size(blockscatter, 2)
    nz = blockstates(p)
    # the incident waves at both stages, `(2 nports)` stage major, from
    # the gathered rates and currents at both stages, `(4 nports)`
    di, dj, dv = Int[], Int[], Float64[]
    xi, xj, xv = Int[], Int[], Float64[]
    offset, prow = 0, 0
    for b in p.blocks
        np = length(b.signal)
        for i in 1:2, q in 1:np
            row = (i - 1)*nports + prow + q
            for l in 1:2
                push!(di, row); push!(dj, (l - 1)*2nports + offset + q)
                push!(dv, entry(gc.ainv, i, l)/h*phi0/(2sqrt(b.R[q])))
            end
            push!(xi, row); push!(xj, (i - 1)*2nports + offset + np + q)
            push!(xv, sqrt(b.R[q])*phi0/(2Lscale))
        end
        offset += 2np
        prow += np
    end
    Wd = sparse(di, dj, dv, 2nports, 4nports)
    Wx = sparse(xi, xj, xv, 2nports, 4nports)
    # qualified: the package's quantum optics library defines a dense
    # `blockdiag` of its own
    Gs = SparseArrays.blockdiag(blockgather, blockgather)
    # the stacked states at the stages from the incident waves and from
    # the states, and the states' update
    ui, uj, uv = Int[], Int[], Float64[]
    zi, zj, zv = Int[], Int[], Float64[]
    ei, ej, ev = Int[], Int[], Float64[]
    fi, fj, fv = Int[], Int[], Float64[]
    prows = cumsum([0; [length(b.signal) for b in p.blocks]])
    for rs in stages
        b = p.blocks[rs.block]
        np, nzb = length(b.signal), size(b.A, 1)
        pr, zb = prows[rs.block], b.zbase
        ucol = (i, q) -> (i - 1)*nports + pr + q
        zrow = (i, k) -> (i - 1)*nz + zb + k
        # the stage algebra of a block whose realization is several
        # filters side by side, a pumped block's, is block diagonal over
        # them, and its structural zeros are not stored
        for i in 1:2, k in 1:nzb
            for l in 1:2, r in 1:np
                v = rs.Zu[(i - 1)*nzb + k, (l - 1)*np + r]
                iszero(v) || (push!(ui, zrow(i, k)); push!(uj, ucol(l, r)); push!(uv, v))
            end
            for l in 1:nzb
                v = rs.Zz[(i - 1)*nzb + k, l]
                iszero(v) || (push!(zi, zrow(i, k)); push!(zj, zb + l); push!(zv, v))
            end
        end
        for k in 1:nzb, l in 1:nzb
            v = rs.Ez[k, l]
            iszero(v) || (push!(ei, zb + k); push!(ej, zb + l); push!(ev, v))
        end
        for k in 1:nzb, i in 1:2, q in 1:np
            v = rs.Eu[k, (i - 1)*np + q]
            iszero(v) || (push!(fi, zb + k); push!(fj, ucol(i, q)); push!(fv, v))
        end
    end
    Zu = sparse(ui, uj, uv, 2nz, 2nports)
    Zz = sparse(zi, zj, zv, 2nz, nz)
    Ezb = sparse(ei, ej, ev, nz, nz)
    Eub = sparse(fi, fj, fv, nz, 2nports)
    Pd = Zu*Wd*Gs
    Px = Zu*Wx*Gs
    Ed = Eub*Wd*Gs
    Ex = Eub*Wx*Gs
    # the output terms: the unconverted response of every block, then
    # each modulated output of the pumped blocks
    Cblk = SparseMatrixCSC{Float64,Int}[]
    terms = Tuple{Int,Int}[]
    ci, cj, cv = Int[], Int[], Float64[]
    for rs in stages
        b = p.blocks[rs.block]
        np, nzb = length(b.signal), size(b.A, 1)
        for q in 1:np, k in 1:nzb
            iszero(b.C[q, k]) || (push!(ci, prows[rs.block] + q); push!(cj, b.zbase + k); push!(cv, b.C[q, k]))
        end
    end
    push!(Cblk, sparse(ci, cj, cv, nports, nz))
    push!(terms, (0, 0))
    for rs in stages
        b = p.blocks[rs.block]
        np, nzb = length(b.signal), size(b.A, 1)
        for (mi, m) in enumerate(b.modulations)
            ci, cj, cv = Int[], Int[], Float64[]
            for q in 1:np, k in 1:nzb
                iszero(m.C[q, k]) || (push!(ci, prows[rs.block] + q); push!(cj, b.zbase + k); push!(cv, m.C[q, k]))
            end
            push!(Cblk, sparse(ci, cj, cv, nports, nz))
            push!(terms, (rs.block, mi))
        end
    end
    t = A -> sparse(transpose(A))
    d = A -> devicesparse(sparse(A), backend)
    SC = [blockscatter*C for C in Cblk]
    Phost = sparse(Pd + Px)
    modulated = sort!(unique!(reduce(vcat, [rowvals(C) for C in Cblk[2:end]]; init = Int[])))
    W = [d((Cblk[j]*Phost[(i - 1)*nz + 1:i*nz, :])[modulated, :]) for i in 1:2 for j in 2:length(Cblk)]
    return RationalCoupling(nz, d(Pd), d(Px), d(Zz), d(Ezb), d(Ed), d(Ex),
        d(t(Px)), d(t(Zz)), d(t(Ezb)), d(t(Ed)), d(t(Ex)),
        [d(M) for M in SC], [d(t(M)) for M in SC], terms, Cblk,
        sparse(blockscatter), Phost, modulated, W,
        [d(abs.(M)) for M in SC], [opnorm(M, Inf) for M in SC])
end

# the weight of every output term at the time `t`: one for the
# unconverted response, and the modulation of each converted output
function modulationweights!(w::AbstractVector, p::TransientProblem, cp::RationalCoupling, t)
    for (j, (bi, mi)) in enumerate(cp.terms)
        w[j] = j == 1 ? 1.0 : modulationweight(p.blocks[bi], p.blocks[bi].modulations[mi], t)
    end
    return w
end


function rationalstages(p::TransientProblem, gc::GaussCoefficients, h)
    tableau = [entry(gc.a, i, j) for i in 1:2, j in 1:2]
    stages = RationalStage[]
    for (k, b) in enumerate(p.blocks)
        nz = size(b.A, 1)
        nz == 0 && continue
        M = Matrix{Float64}(I, 2nz, 2nz) - h .* kron(tableau, b.A)
        F = lu(M)
        Zz = F \ kron(ones(2, 1), Matrix{Float64}(I, nz, nz))
        Zu = F \ (h .* kron(tableau, b.B))
        Asum = hcat(b.A, b.A)
        Ez = Matrix{Float64}(I, nz, nz) .+ (h/2) .* (Asum*Zz)
        Eu = (h/2) .* (Asum*Zu .+ hcat(b.B, b.B))
        push!(stages, RationalStage(k, Zz, Zu, Ez, Eu))
    end
    return stages
end

# the rational blocks' frequency dependent hybrid entries at the complex
# frequency `s` as a sparse matrix on the state: their unconverted
# response, which the frozen stage operator carries at the stage
# frequency `mu/h` (a pumped block's converted outputs are the stage
# correction's, see [`StageCorrection`](@ref))
function rationalmatrix(p::TransientProblem, s, Lscale, n)
    rows, cols, vals = Int[], Int[], ComplexF64[]
    for b in p.blocks
        size(b.A, 1) == 0 && continue
        np = length(b.signal)
        Sr = b.C*((s*I - b.A) \ b.B)
        Bb = -s*Lscale .* Sr .* transpose(1 ./ sqrt.(b.R))
        Cb = -Sr .* transpose(sqrt.(b.R))
        for q in 1:np, r in 1:np
            push!(rows, b.auxbase + q); push!(cols, b.auxbase + r); push!(vals, Cb[q, r])
            b.signal[r] > 0 && (push!(rows, b.auxbase + q); push!(cols, b.signal[r]); push!(vals, Bb[q, r]))
            b.ref[r] > 0 && (push!(rows, b.auxbase + q); push!(cols, b.ref[r]); push!(vals, -Bb[q, r]))
        end
    end
    return sparse(rows, cols, vals, n, n)
end

# The values of one orientation of the two stages' matrix (see
# `StagePlan`), the matrix or its transpose, in the compressed column order
# of that orientation, with its pattern `colptr` and `rowval`: the constant
# part `base`, and gathered per stored entry in a fixed order, so that a
# matrix is the same whatever assembles it, the junction terms, `jcoef`
# times the stiffness of the junction and stage `jrow` (the junctions at the
# first stage, then at the second), and the converted outputs, `bval` times
# the weight `brow` of the term and stage (the terms at the first stage,
# then at the second).
struct StageGather{VF, VI}
    colptr::VI
    rowval::VI
    base::VF
    jptr::VI
    jrow::VI
    jcoef::VF
    bptr::VI
    brow::VI
    bval::VF
end

"""
    StagePlan

The real matrix of a Gauss-Legendre step's two stages, which the tangent
and the adjoint solve. On the stage stacked unknowns `[stage 1; stage 2]`
its block `(i, l)` is

    A2_il C/h^2 + A1_il G/h + delta_il (L + J_i) - B_il,

with `A1` and `A2` the entries of the inverse tableau and of its square,
`J_i` the junction stamp at stage `i`'s phases, so that each stage has its
own stiffness, and `B_il` the rational blocks' coupling of the stages at
stage `i`'s weights, the converted outputs of a pumped block among them.
Its pattern is the step operator's with the junction pairs and the blocks'
coupling, in every block. The plan holds the pattern on the host, the
values' gather in the matrix's own order (see `StageGather`), which a
host factorization and the adjoint's device factorization read, and on a
device in its transpose's order, which the device reads the matrix in for
the tangent; the counts of the junctions and of the output terms; whether
the matrix is constant, without junctions and converted outputs; and the
host's fill reducing ordering of the pattern, each node's two stages
together in the step operator's order, rows and columns, which differ
where the step operator's diagonal has gaps.
"""
struct StagePlan{G}
    n::Int
    nj::Int
    nterms::Int
    constant::Bool
    pattern::SparseMatrixCSC{Float64,Int}
    natural::G
    transposed::Union{Nothing, G}
    ordering::FactorizationCache
    # the only constructor, and it takes the parameter, which a host's
    # plan leaves to its plain field alone
    StagePlan{G}(n, nj, nterms, constant, pattern, natural, transposed, ordering) where {G} =
        new{G}(n, nj, nterms, constant, pattern, natural, transposed, ordering)
end

# the plan of a system's two stages' matrix from its scaled matrices `C`,
# `G`, `L`, the step operator `K`, all on the host, the junction incidence
# `RJ` and the blocks' coupling or nothing, on `backend`
function stageplan(C, G, L, K, RJ, gc::GaussCoefficients, h, coupling, backend)
    n, nj = size(K, 1), size(RJ, 1)
    # the blocks' coupling of stage `l` into stage `i` through each output
    # term: the unconverted response at the weight one, then the converted
    # outputs at the stage's weights
    nz = isnothing(coupling) ? 0 : coupling.nstates
    nterms = isnothing(coupling) ? 0 : length(coupling.Cblk)
    blocks = [sparse(coupling.Shost*coupling.Cblk[j]*coupling.Phost[(i - 1)*nz + 1:i*nz, (l - 1)*n + 1:l*n])
        for j in 1:nterms, i in 1:2, l in 1:2]
    # every entry as row, column and value in the matrix's coordinates:
    # the constant part, the junction terms with their stiffness's row, and
    # the converted outputs with their weight's
    rb, cb, vb = Int[], Int[], Float64[]
    rj, cj, jrow, jcoef = Int[], Int[], Int[], Float64[]
    rt, ct, brow, bval = Int[], Int[], Int[], Float64[]
    place! = (M, i, l, s) -> begin
        I, J, V = SparseArrays.findnz(M)
        append!(rb, I .+ (i - 1)*n); append!(cb, J .+ (l - 1)*n); append!(vb, s .* V)
        nothing
    end
    for i in 1:2, l in 1:2
        place!(C, i, l, entry(gc.ainv2, i, l)/h^2)
        place!(G, i, l, entry(gc.ainv, i, l)/h)
        i == l && place!(L, i, l, 1.0)
        nterms > 0 && place!(blocks[1, i, l], i, l, -1.0)
    end
    RJt = sparse(transpose(RJ))
    for i in 1:2, k in 1:nj, pa in nzrange(RJt, k), pb in nzrange(RJt, k)
        push!(rj, rowvals(RJt)[pa] + (i - 1)*n); push!(cj, rowvals(RJt)[pb] + (i - 1)*n)
        push!(jrow, k + (i - 1)*nj); push!(jcoef, nonzeros(RJt)[pa]*nonzeros(RJt)[pb])
    end
    for j in 2:nterms, i in 1:2, l in 1:2
        I, J, V = SparseArrays.findnz(blocks[j, i, l])
        append!(rt, I .+ (i - 1)*n); append!(ct, J .+ (l - 1)*n)
        append!(brow, fill(j + (i - 1)*nterms, length(I))); append!(bval, .-V)
    end
    # the pattern: the step operator's, the junction pairs' and the
    # blocks' coupling's, in every block
    units = (I, J) -> sparse(mod1.(I, n), mod1.(J, n), zeros(length(I)), n, n)
    P = spaddkeepzeros(spaddkeepzeros(spaddkeepzeros(SparseMatrixCSC(n, n, copy(SparseArrays.getcolptr(K)), copy(rowvals(K)),
        zeros(nnz(K))), units(rb, cb)), units(rj, cj)), units(rt, ct))
    pattern = sparse([P P; P P])
    gather = transposed -> stagegather(pattern, transposed, (rb, cb, vb), (rj, cj, jrow, jcoef), (rt, ct, brow, bval), backend)
    natural = gather(false)
    return StagePlan{typeof(natural)}(n, nj, nterms, nj == 0 && isempty(bval), pattern, natural,
        backend isa CPU ? nothing : gather(true), FactorizationCache())
end

# the gather of the values of the matrix of pattern `pattern`, or of its
# transpose, from the entries of its constant part, its junction terms and
# its converted outputs (see `StageGather`)
function stagegather(pattern::SparseMatrixCSC, transposed::Bool, constantpart, junctions, converted, backend)
    S = transposed ? sparse(transpose(pattern)) : pattern
    colptr, rowval = SparseArrays.getcolptr(S), rowvals(S)
    position = (r, c) -> transposed ? storedposition(colptr, rowval, c, r) : storedposition(colptr, rowval, r, c)
    rb, cb, vb = constantpart
    base = zeros(nnz(S))
    for k in eachindex(rb)
        base[position(rb[k], cb[k])] += vb[k]
    end
    # each list grouped by its entries, in the order it was formed
    grouped = (rows, cols) -> begin
        at = [position(rows[k], cols[k]) for k in eachindex(rows)]
        order = sortperm(at; alg = MergeSort)
        ptr = zeros(Int, nnz(S) + 1)
        ptr[1] = 1
        for q in at
            ptr[q + 1] += 1
        end
        cumsum!(ptr, ptr)
        ptr, order
    end
    rj, cj, jrow, jcoef = junctions
    jptr, jorder = grouped(rj, cj)
    rt, ct, brow, bval = converted
    bptr, border = grouped(rt, ct)
    host = backend isa CPU
    index = x -> host ? Vector{Int}(x) : tobackend(backend, Vector{Int32}(x))
    value = x -> tobackend(backend, Vector{Float64}(x))
    return StageGather(index(colptr), index(rowval), value(base), index(jptr), index(jrow[jorder]), value(jcoef[jorder]),
        index(bptr), index(brow[border]), value(bval[border]))
end

# the fill reducing ordering of the two stages' pattern from the step
# operator's, `nothing` for KLU's own: each node's two stages one after the
# other, the nodes in the step operator's order, which eliminates the
# pattern's two by two blocks as that order eliminates its entries. Where
# the step operator's rows are matched to other columns, its diagonal
# having gaps (see `fillordering`), each stage's rows are matched to that
# stage's columns alike, and the fill is the pattern's so permuted.
function stageordering(ordering, pattern::SparseMatrixCSC, n::Int)
    isnothing(ordering) && return nothing
    perm = orderingpermutation(ordering)
    stages = p -> vec(transpose(hcat(p, p .+ n)))
    if ordering.rows == perm
        pairs = stages(perm)
        fill, _ = symbolicfill(_symmetricpattern(pattern), pairs)
        return FillOrdering(pairs, fill)
    end
    # the column matched to each row, in either stage
    cols = similar(perm)
    cols[ordering.rows] = perm
    rows = stages(ordering.rows)
    fill, _ = symbolicfill(_symmetricpattern(pattern[:, vcat(cols, cols .+ n)]), rows)
    return FillOrdering(stages(perm), fill, rows)
end

# What a Gauss-Legendre system holds beyond the trapezoidal one: the
# coefficients, the imaginary part of the stage matrix on the Jacobian's
# pattern, the complex matrix the factorization reads, the stage algebra
# of the rational blocks, and the plan of the two stages' matrix the
# tangent and the adjoint solve.
struct GaussStage{V, M, R, SM, P}
    coefficients::GaussCoefficients
    imvals::V
    cjacobian::M
    # the complex entries of the rational blocks' unconverted responses
    # at the stage frequency on the pattern, and the blocks' stage
    # algebra: the grouped coupling of the blocks, or nothing without
    # any, as a union within the backend's sparse matrix type, so that
    # the stage's type, and with it the system's, the stepper's and every
    # response's, is the backend's alone, a circuit with a block runs on
    # the code compiled for one without, and a presence check is a branch
    # rather than a specialization; and whether any block is pumped, which
    # is when the step's solves carry a correction
    rationalvals::R
    coupling::Union{Nothing, RationalCoupling{SM}}
    pumped::Bool
    # The fill reducing ordering of the stage matrix's pattern for a host
    # factorization, chosen once for the system: every fresh factorization
    # of a condition, a chunk or a checkpoint window on the system takes
    # it, and only reads it, so the chunks of a batch share it across
    # threads. Empty on a device, which orders its batch itself.
    ordering::FactorizationCache
    # the plan of the two stages' matrix, whose host ordering derives from
    # the stage matrix's
    stages::P
    # The only constructor, and it takes the parameters: `SM` appears in
    # the union field alone, so a circuit without a coupling passes
    # `nothing` and leaves it with nothing to infer from. `gaussstage`
    # reads it off the backend.
    GaussStage{V, M, R, SM, P}(coefficients, imvals, cjacobian, rationalvals, coupling, pumped, ordering,
        stages) where {V, M, R, SM, P} =
        new{V, M, R, SM, P}(coefficients, imvals, cjacobian, rationalvals, coupling, pumped, ordering, stages)
end

# the stage with its union field's type taken from the backend, and on
# the host the ordering `factorization` chooses for its pattern, with the
# two stages' ordering from it
function gaussstage(gc, imvals, cjacobian, rationalvals, coupling, backend, pumped, factorization, stages::StagePlan)
    SM = typeof(devicesparse(sparse(zeros(1, 1)), backend))
    ordering = FactorizationCache()
    if backend isa CPU
        chosen = fillordering(factorization, cjacobian)
        seedordering!(ordering, cjacobian, chosen)
        seedordering!(stages.ordering, stages.pattern, stageordering(chosen, stages.pattern, stages.n))
    end
    return GaussStage{typeof(imvals), typeof(cjacobian), typeof(rationalvals), SM, typeof(stages)}(gc, imvals, cjacobian,
        rationalvals, coupling, pumped, ordering, stages)
end

# The entries of a sparse matrix `A` placed on the pattern of the real
# Jacobian, which holds every entry of the step matrix and of the blocks'
# rows, as the values of the pattern; on a device the pattern is stored
# transposed, as the column structure of the transpose, and `A` is
# placed transposed on it.
patterncolumns(Jrs::SparseMatrixCSC) = SparseArrays.getcolptr(Jrs), rowvals(Jrs)
patterncolumns(Jrs::DeviceSparsePattern) = Array(Jrs.colptr), Array(Jrs.rowval)
function patternvalues(A::SparseMatrixCSC{T}, Jrs, transposed::Bool) where {T}
    colptr, rowval = patterncolumns(Jrs)
    n = length(colptr) - 1
    pattern = SparseMatrixCSC(n, n, Vector{Int}(colptr), Vector{Int}(rowval), ones(Bool, length(rowval)))
    B = dropzeros(transposed ? sparse(transpose(A)) : A)
    vals = zeros(T, length(rowval))
    vals[sparseaddmap(pattern, B)] .= nonzeros(B)
    return vals
end

# a stage of a stage pair: the column of an `(n, 2)` array, or the
# contiguous matrix of an `(n, ndir, 2)` array of directions
stage(a::AbstractMatrix, i) = view(a, :, i)
stage(a::AbstractArray{<:Any,3}, i) = view(a, :, :, i)
# The cubic Lagrange stencil that reads a current at a stage time of the
# step from recorded time `k` to `k + 1` off the recorded grid: the four
# grid indices and their weights, one sided at the ends of the record, and
# linear on a record too short for four points, whose last two weights
# are zero. Fourth order, so a smooth current keeps the method's order.
# Tuples, which a step reads without allocating.
function gaussstencil(k::Int, nt::Int, c::Float64)
    if nt < 4
        return (k, k + 1, k, k), (1 - c, c, 0.0, 0.0)
    end
    first = k == 1 ? 1 : k == nt - 1 ? nt - 3 : k - 1
    s = c + (k - first)
    weight = j -> prod((s - m)/(j - m) for m in 0:3 if m != j)
    return (first, first + 1, first + 2, first + 3), (weight(0), weight(1), weight(2), weight(3))
end
