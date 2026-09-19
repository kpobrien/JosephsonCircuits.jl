# The two stage Gauss-Legendre collocation on the second order system: two
# stage fluxes per step, solved together by the shared Newton engine on the
# true stage equations, with the correction from one complex factorization.
# The stage matrix of the tableau has the conjugate eigenvalues
# `mu = 3 ± i sqrt(3)`, so with one frozen junction stiffness the two
# stages decouple into the complex system `(mu/h)^2 C + (mu/h) G + L + J*`
# and its conjugate, of which one is solved; that is a simplified Newton,
# converging linearly on the true residual, and the same operator
# preconditions the exact solves of the tangent and the adjoint of the
# stage equations. The complex matrix is the harmonic balance real
# Jacobian plan's assembly of the real part on its pattern, plus a constant
# imaginary part on the same pattern.

"""
    GaussLegendre()

The two stage Gauss-Legendre collocation on the flux and its rate, the
default of [`transientsolve`](@ref): fourth order, A-stable, symplectic,
and free of numerical damping, so a lossless LC oscillation keeps its
energy and a resonator's frequency is warped by `(2 pi f dt)^4/720`
rather than the trapezoidal rule's `(2 pi f dt)^2/12`. Each step solves
the two stage equations together, by Newton on one
complex factorization of the stage matrix at a frozen junction
stiffness, refreshed as the trapezoidal rule's is. Not L-stable: an
unresolved fast mode is not damped, as the trapezoidal rule does not damp
it.

Along a direction of the state without capacitance the equations are
algebraic. Where a resistor acts along it the rate is what the
constraint determines and converges at second order. Where none does the
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
steps.
"""
struct GaussLegendre <: AbstractTransientIntegrator end

# The tableau, and what a step reads of it: the abscissae `c`, the inverse
# `A^{-1}` and its square as column major tuples, `A^{-1} 1`, the endpoint
# weights `b' A^{-1}` of the flux and `b' A^{-2}` of the rate, the
# predictor of the next stages from the last increments, and the
# eigenvalue `mu` of `A^{-1}` in the upper half plane with the entries of
# its eigenvector matrix `T = [t conj(t)]` and of the inverse the stage
# transform reads, and the derivative at the endpoint of the cubic
# through the state, the two stages and the endpoint, as weights on the
# four, for the rate along an algebraic direction a block is on.
struct GaussCoefficients
    c::NTuple{2,Float64}
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
    return GaussCoefficients(c, (Ainv[1, 1], Ainv[2, 1], Ainv[1, 2], Ainv[2, 2]),
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
# `E_d = E_u W_d G`, `E_x = E_u W_x G`; their transposes; and on the
# host the output matrices `C_j` for the resting waves. A circuit
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
    Pdt::M
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
    # solve of a pumped block (see [`StageCorrection`](@ref))
    Shost::SparseMatrixCSC{Float64,Int}
    Phost::SparseMatrixCSC{Float64,Int}
    # the port rows with a modulated output term, the only rows the
    # stage correction acts on
    modulated::Vector{Int}
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
    Gs = blockdiag(blockgather, blockgather)
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
    return RationalCoupling(nz, d(Pd), d(Px), d(Zz), d(Ezb), d(Ed), d(Ex),
        d(t(Pd)), d(t(Px)), d(t(Zz)), d(t(Ezb)), d(t(Ed)), d(t(Ex)),
        [d(M) for M in SC], [d(t(M)) for M in SC], terms, Cblk,
        sparse(blockscatter), sparse(Pd + Px),
        sort!(unique!(reduce(vcat, [rowvals(C) for C in Cblk[2:end]]; init = Int[]))))
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
    tableau = [1/4 1/4-sqrt(3)/6; 1/4+sqrt(3)/6 1/4]
    stages = RationalStage[]
    for (k, b) in enumerate(p.blocks)
        nz, np = size(b.A, 1), length(b.signal)
        nz == 0 && continue
        I2 = Matrix{Float64}(I, 2, 2)
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
# frequency `s` as a sparse matrix on the state
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

# The transfer of every rational block at the complex stage frequency
# `s`, `C (s I - A)^(-1) B` per output term: the first the unconverted
# response, then each modulated output of a pumped block, in the order
# of the coupling's terms, so that the stage operator can carry the
# converted coupling at a step's weights (see [`rationalvalues!`](@ref)).
# `nothing` for a block without states.
function rationalstageterms(p::TransientProblem, s)
    terms = Vector{Union{Nothing,Vector{Matrix{ComplexF64}}}}(undef, length(p.blocks))
    for (bi, b) in enumerate(p.blocks)
        if size(b.A, 1) == 0
            terms[bi] = nothing
            continue
        end
        X = (s*I - b.A) \ b.B
        terms[bi] = vcat([b.C*X], [m.C*X for m in b.modulations])
    end
    return terms
end

"""
    rationalvalues!(vals, p, s, Lscale, Jrs, transposed, terms, weights)

Overwrite `vals`, the entries of the rational blocks on the Jacobian's
pattern at the complex stage frequency `s`, from the blocks' stage
transfers `terms` (see `rationalstageterms`) weighted: the
unconverted response with one, and each modulated output of a pumped
block with its entry of `weights`, in the order of the coupling's terms
after the first, or with nothing of the modulated outputs when `weights`
is `nothing`. A pumped block's stage operator is refreshed this way at
every step with the mean of its two stages' weights, since its
converted coupling can be as large as its unconverted one, which the
frozen operator of the simplified Newton would not converge without.
"""
function rationalvalues!(vals, p::TransientProblem, s, Lscale, Jrs, transposed, terms, weights)
    fill!(vals, 0)
    colptr, rowval = patterncolumns(Jrs)
    place = (r, c, val) -> begin
        pos = 0
        rr, cc = transposed ? (c, r) : (r, c)
        for k in colptr[cc]:colptr[cc + 1] - 1
            rowval[k] == rr && (pos = k; break)
        end
        pos > 0 || error("an entry of a rational block is absent from the Jacobian's pattern.")
        vals[pos] += val
    end
    j = 1
    for (bi, b) in enumerate(p.blocks)
        isnothing(terms[bi]) && continue
        n = length(b.signal)
        Sr = copy(terms[bi][1])
        for (mi, m) in enumerate(b.modulations)
            j += 1
            isnothing(weights) && continue
            Sr .+= weights[j] .* terms[bi][mi + 1]
        end
        Bb = -s*Lscale .* Sr .* transpose(1 ./ sqrt.(b.R))
        Cb = -Sr .* transpose(sqrt.(b.R))
        for q in 1:n
            row = b.auxbase + q
            for r in 1:n
                place(row, b.auxbase + r, Cb[q, r])
                b.signal[r] > 0 && place(row, b.signal[r], Bb[q, r])
                b.ref[r] > 0 && place(row, b.ref[r], -Bb[q, r])
            end
        end
    end
    return vals
end

# What a Gauss-Legendre system holds beyond the trapezoidal one: the
# coefficients, the imaginary part of the stage matrix on the Jacobian's
# pattern, the complex matrix the factorization reads, and the stage
# algebra of the rational blocks.
struct GaussStage{V, M, R, SM}
    coefficients::GaussCoefficients
    imvals::V
    cjacobian::M
    # the complex entries of the rational blocks at the stage frequency
    # on the pattern, and the blocks' stage algebra: the grouped coupling
    # of the blocks, or nothing without any, as a union within the
    # backend's sparse matrix type, so that the stage's type, and with it
    # the system's, the stepper's and every response's, is the backend's
    # alone, a circuit with a block runs on the code compiled for one
    # without, and a presence check is a branch rather than a
    # specialization
    rationalvals::R
    rational::Vector{RationalStage}
    coupling::Union{Nothing, RationalCoupling{SM}}
    # the stage transfers of the blocks' output terms with the host
    # pattern of the Jacobian, the host copy of the values and whether
    # they are laid out transposed, for the refresh of a pumped block's
    # coupling at every step; and whether any block is pumped, which is
    # when the refresh happens
    stageterms::Any
    hostvals::Vector{ComplexF64}
    transposedvals::Bool
    pumped::Bool
    # The only constructor, and it takes the parameters: `SM` appears in
    # the union field alone, so a circuit without a coupling passes
    # `nothing` and leaves it with nothing to infer from. `gaussstage`
    # reads it off the backend.
    GaussStage{V, M, R, SM}(coefficients, imvals, cjacobian,
            rationalvals, rational, coupling, stageterms, hostvals, transposedvals, pumped) where {V, M, R, SM} =
        new{V, M, R, SM}(coefficients, imvals, cjacobian,
            rationalvals, rational, coupling, stageterms, hostvals, transposedvals, pumped)
end

# the stage with its union field's type taken from the backend
function gaussstage(gc, imvals, cjacobian, rationalvals, stages, coupling, backend,
        stageterms, hostvals, transposedvals, pumped)
    SM = typeof(devicesparse(sparse(zeros(1, 1)), backend))
    return GaussStage{typeof(imvals), typeof(cjacobian), typeof(rationalvals), SM}(gc, imvals, cjacobian,
        rationalvals, stages, coupling, stageterms, hostvals, transposedvals, pumped)
end

# the imaginary part of the stage matrix, `Im((mu/h)^2) C + Im(mu/h) G`,
# placed on the pattern of the real Jacobian, which holds every entry of
# `C` and `G`; on a device the pattern is stored transposed, as the
# column structure of the transpose
patterncolumns(Jrs::SparseMatrixCSC) = SparseArrays.getcolptr(Jrs), rowvals(Jrs)
patterncolumns(Jrs::DeviceSparsePattern) = Array(Jrs.colptr), Array(Jrs.rowval)
function gaussimaginary(gc::GaussCoefficients, h, C, G, Jrs, transposed)
    ci, gi = imag((gc.mu/h)^2), imag(gc.mu/h)
    Kim = ci .* C .+ gi .* G
    imvals = zeros(nnz(Jrs))
    colptr, rowval = patterncolumns(Jrs)
    for j in axes(Kim, 2), q in nzrange(Kim, j)
        val = nonzeros(Kim)[q]
        iszero(val) && continue
        i = rowvals(Kim)[q]
        r, c = transposed ? (j, i) : (i, j)
        pos = 0
        for k in colptr[c]:colptr[c + 1] - 1
            rowval[k] == r && (pos = k; break)
        end
        pos > 0 || error("an entry of the stage matrix is absent from the Jacobian's pattern.")
        imvals[pos] = val
    end
    return imvals
end

# a stage of a stage pair: the column of an `(n, 2)` array, or the
# contiguous matrix of an `(n, ndir, 2)` array of directions
stage(a::AbstractMatrix, i) = view(a, :, i)
stage(a::AbstractArray{<:Any,3}, i) = view(a, :, :, i)
# The cubic Lagrange stencil that reads a current at a stage time of the
# step from recorded time `k` to `k + 1` off the recorded grid: the four
# grid indices and their weights, one sided at the ends of the record, and
# linear on a record too short for four points. Fourth order, so a smooth
# current keeps the method's order.
function gaussstencil(k::Int, nt::Int, c::Float64)
    if nt < 4
        return [k, k + 1], [1 - c, c]
    end
    first = k == 1 ? 1 : k == nt - 1 ? nt - 3 : k - 1
    indices = [first, first + 1, first + 2, first + 3]
    s = c + (k - first)
    weights = [prod((s - m)/(j - m) for m in 0:3 if m != j) for j in 0:3]
    return indices, weights
end
