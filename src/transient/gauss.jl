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

The two stage Gauss-Legendre collocation on the flux and its rate: fourth
order, A-stable, symplectic, and free of numerical damping, so a lossless
LC oscillation keeps its energy and a resonator's frequency is warped by
`(2 pi f dt)^4/720` rather than the trapezoidal rule's `(2 pi f dt)^2/12`.
Each step solves the two stage equations together, by Newton on one
complex factorization of the stage matrix at a frozen junction
stiffness, refreshed as the trapezoidal rule's is. Not L-stable: an
unresolved fast mode is not damped, as the trapezoidal rule does not damp
it.

Along a direction of the state without capacitance the equations are
algebraic. Where a resistor acts along it the rate is what the
constraint determines and converges at second order; where none does
and the constraint is nonlinear or driven, a junction or a source on a
node without capacitor or resistor to ground, the plain rule's rate does
not converge at all, so each step projects its endpoint onto the
constraint along those directions and takes the rate along them from
the derivative of the cubic through the state, the two stages and the
projected endpoint, all on the constraint, which is third order. A
linear constraint no source drives is an invariant of the rule and needs
neither. The tangent and the adjoint of a Gauss-Legendre solve
differentiate the full stage equations, the projection included, and are
exact for the recorded steps.
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
# four, for the rate along an algebraic direction.
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
# the gathered increments and values, `C Z_u`, `C Z_z` the reflected
# waves from the incident waves and the states, and `E_z`, `E_u` the
# states' update: `M_d = S C Z_u W_d G`, `M_x = S C Z_u W_x G`,
# `M_s = S C Z_z`, `E_d = E_u W_d G`, `E_x = E_u W_x G`, their
# transposes, and `R' = (S_1 C)'` for the cotangent of the resting waves
# the endpoint reading sees, `S_1` the scatter of one stage.
struct RationalCoupling{M}
    nstates::Int
    Md::M
    Mx::M
    Ms::M
    Ez::M
    Ed::M
    Ex::M
    Mdt::M
    Mxt::M
    Mst::M
    Ezt::M
    Edt::M
    Ext::M
    Rt::M
    # the output matrices of every block on the host, `(nports, nstates)`,
    # for the resting waves
    Cblk::SparseMatrixCSC{Float64,Int}
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
    Ss = blockdiag(blockscatter, blockscatter)
    # the reflected waves from the incident waves and the states, the
    # states' update, and the output matrices
    ui, uj, uv = Int[], Int[], Float64[]
    zi, zj, zv = Int[], Int[], Float64[]
    ei, ej, ev = Int[], Int[], Float64[]
    fi, fj, fv = Int[], Int[], Float64[]
    ci, cj, cv = Int[], Int[], Float64[]
    prows = cumsum([0; [length(b.signal) for b in p.blocks]])
    for rs in stages
        b = p.blocks[rs.block]
        np, nzb = length(b.signal), size(b.A, 1)
        pr, zb = prows[rs.block], b.zbase
        ucol = (i, q) -> (i - 1)*nports + pr + q
        CZu = kron(Matrix(1.0I, 2, 2), b.C)*rs.Zu
        CZz = kron(Matrix(1.0I, 2, 2), b.C)*rs.Zz
        for i in 1:2, q in 1:np
            for l in 1:2, r in 1:np
                push!(ui, ucol(i, q)); push!(uj, ucol(l, r)); push!(uv, CZu[(i - 1)*np + q, (l - 1)*np + r])
            end
            for k in 1:nzb
                push!(zi, ucol(i, q)); push!(zj, zb + k); push!(zv, CZz[(i - 1)*np + q, k])
            end
        end
        for k in 1:nzb, l in 1:nzb
            push!(ei, zb + k); push!(ej, zb + l); push!(ev, rs.Ez[k, l])
        end
        for k in 1:nzb, i in 1:2, q in 1:np
            push!(fi, zb + k); push!(fj, ucol(i, q)); push!(fv, rs.Eu[k, (i - 1)*np + q])
        end
        for q in 1:np, k in 1:nzb
            push!(ci, pr + q); push!(cj, zb + k); push!(cv, b.C[q, k])
        end
    end
    CZu = sparse(ui, uj, uv, 2nports, 2nports)
    CZz = sparse(zi, zj, zv, 2nports, nz)
    Ezb = sparse(ei, ej, ev, nz, nz)
    Eub = sparse(fi, fj, fv, nz, 2nports)
    Cblk = sparse(ci, cj, cv, nports, nz)
    Md = Ss*CZu*Wd*Gs
    Mx = Ss*CZu*Wx*Gs
    Ms = Ss*CZz
    Ed = Eub*Wd*Gs
    Ex = Eub*Wx*Gs
    Rt = sparse(transpose(blockscatter*Cblk))
    t = A -> sparse(transpose(A))
    d = A -> devicesparse(sparse(A), backend)
    return RationalCoupling(nz, d(Md), d(Mx), d(Ms), d(Ezb), d(Ed), d(Ex),
        d(t(Md)), d(t(Mx)), d(t(Ms)), d(t(Ezb)), d(t(Ed)), d(t(Ex)), d(Rt), Cblk)
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

# the complex entries of the rational blocks at the complex frequency
# `s` on the pattern: on the block's rows, `-s Lscale C (s I - A)^(-1) B R^(-1/2)`
# on its port nodes and `-C (s I - A)^(-1) B R^(1/2)` on its port
# currents, the frequency dependent parts of the hybrid stamp
function rationalvalues(p::TransientProblem, s, Lscale, Jrs, transposed)
    vals = zeros(ComplexF64, nnz(Jrs))
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
    for b in p.blocks
        size(b.A, 1) == 0 && continue
        n = length(b.signal)
        Sr = b.C*((s*I - b.A) \ b.B)
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

# The projection of a step's endpoint onto the algebraic constraints.
# Along an algebraic direction `z` of the problem, where neither
# capacitance nor conductance acts, the equation `z' (L x + J(x)) = z' b(t)`
# constrains the flux alone; the stages satisfy it, and the rule's endpoint,
# the collocation quadratic extrapolated, does not, by an error the plain
# rule never damps and the rate update divides by the step. A linear
# constraint without a source is an invariant of the rule, since the
# endpoint weights sum to zero; a constraint with a junction on its nodes
# or a source on them is not, and those directions are the columns of `Z`.
# What the projection reads: `Z` and its transpose on the backend, `Z' L`,
# the rows of the junction incidence of the junctions on those nodes
# (`RJp`, and its transpose) with their coefficients, `RJ Z` and `Z' L Z`
# on the host for the small Jacobians `Z' (L + J'(x)) Z`, the inverse
# squared norms of the columns for the rate, and `Z'` of the drive
# injection and the constant current for the residual.
#
# The endpoint's other algebraic unknowns are index one: the rate along
# an inertialess direction a resistor reaches, and the port current of a
# scattering block, which the equations along those directions and rows
# determine linearly given the flux and the rates elsewhere, and which
# the rule's endpoint formulas give only to second order. They are read
# from their equations at the endpoint, `Q (G v + L x + J(x) - b) = 0`
# with `Q` the stack of the resistive directions `Zr'` and the block rows
# `Ea'`, for the rate along `Zr` and the current on `Ea`, whose Jacobian
# `M = [Q G Zr, Q L Ea]` is constant. The rows may not determine every
# unknown: a through, whose `I + S` is singular, fixes the sum of its two
# currents and leaves their difference to the nodes it joins, so the
# reading is the minimum norm correction through the pseudoinverse of
# `M`, which moves the combinations the rows determine and nothing else.
# Neither unknown enters the next step's stage equations, so this is a
# reading of the endpoint, fourth order where the rest of the state is,
# which the port waves and the histories of the lines read.
struct GaussProjection{M, MP, VP}
    # the algebraic directions as their supports, the projected junctions,
    # the directions `Z` as columns, their transpose, the constraints'
    # rows `Zt` times the stiffness, the constraints on the host, the
    # constraints' rows as columns `Zl`, the junction rows of the
    # incidence and their transpose, the junction phases along the
    # directions and along the constraints, the projected junctions'
    # coefficients, the stiffness along the constraints and directions,
    # the rate extraction along the directions and its transpose, and the
    # constraints' rows of the injection, the constant current, the
    # lines' forcing and the blocks' resting waves
    directions::Vector{Vector{Int}}
    pj::Vector{Int}
    Z::M
    Zt::M
    ZtL::M
    Zth::SparseMatrixCSC{Float64,Int}
    Zl::M
    RJp::MP
    RJpt::MP
    RJZ::Matrix{Float64}
    RJZl::Matrix{Float64}
    lmoljp::Vector{Float64}
    lmoljpdev::VP
    # the current-phase relations of the projected junctions, in the order
    # of `pj`. The projection evaluates them on the host, from `hphi`, and
    # moves the result where a device array needs it, so one host table
    # serves both.
    relationsp::JunctionRelations{Matrix{Float64},Vector{Bool}}
    ZLZ::Matrix{Float64}
    Zrate::M
    Zratet::M
    Ztinj::SparseMatrixCSC{Float64,Int}
    Ztconstant::Vector{Float64}
    Ztline::SparseMatrixCSC{Float64,Int}
    Ztblock::SparseMatrixCSC{Float64,Int}
    # the reading of the endpoint: the capacitor free islands and the
    # block current rows, their indicators, the rows `Q = [Zr'; Ea']`,
    # the rate system's pseudoinverse on its range, and the rows of the
    # drives, the constant, the lines' forcing and the resting waves
    readrows::Vector{Vector{Int}}
    auxrows::Vector{Int}
    Zr::M
    Ea::M
    Q::M
    Qt::M
    QG::M
    QL::M
    RJQ::Matrix{Float64}
    Minv::Matrix{Float64}
    Qinj::SparseMatrixCSC{Float64,Int}
    Qconst::Vector{Float64}
    Qline::SparseMatrixCSC{Float64,Int}
    Qblock::SparseMatrixCSC{Float64,Int}
end

# The projection of a problem under the scaled matrices, or nothing
# where no direction needs one. A direction whose constraint no junction,
# drive, line or block enters is linear with a constant right hand side,
# which the stages keep exactly, so it is left alone; the others are
# projected. The reading covers every capacitor free island and block
# current row through the rate system's pseudoinverse, which determines
# what the equations determine and leaves the algebraic directions to
# the projection.
function gaussprojection(p::TransientProblem, G::SparseMatrixCSC, L::SparseMatrixCSC, RJ::SparseMatrixCSC,
        lmolj::Vector{Float64}, injection::SparseMatrixCSC, constant::Vector{Float64}, lineinjection::SparseMatrixCSC,
        blockscatter::SparseMatrixCSC, backend)
    n = length(p)
    Zall, Ztall, Rall = sparse(p.directions), sparse(p.constraints), sparse(p.rateextraction)
    droptol!(Zall, 1e-14); droptol!(Ztall, 1e-14); droptol!(Rall, 1e-14)
    touched = Int[]
    for c in 1:size(Zall, 2)
        junctions = nnz(RJ*Zall[:, c]) > 0
        row = Ztall[c:c, :]
        driven = nnz(row*injection) > 0 || nnz(row*lineinjection) > 0 || nnz(row*blockscatter) > 0
        (junctions || driven) && push!(touched, c)
    end
    directions = p.algebraic[touched]
    Zh, Zth, Zrateh = Zall[:, touched], Ztall[touched, :], Rall[touched, :]
    Zlh = sparse(transpose(Zth))
    # the reading: the capacitor free islands and the block current rows
    readrows = [z for z in p.inertialess if all(<=(p.Nnodal), z)]
    auxrows = [b.auxbase + q for b in p.blocks for q in eachindex(b.signal)]
    (isempty(directions) && isempty(readrows) && isempty(auxrows)) && return nothing
    indicator = dirs -> sparse(reduce(vcat, dirs; init = Int[]), reduce(vcat, [fill(c, length(z)) for (c, z) in enumerate(dirs)]; init = Int[]),
        ones(sum(length, dirs; init = 0)), n, length(dirs))
    support = [i for c in 1:size(Zh, 2) for i in findnz(Zh[:, c])[1]]
    pj = sort!(unique!([i for node in vcat(support, reduce(vcat, readrows; init = Int[])) for i in rowvals(RJ)[nzrange(RJ, node)]]))
    RJph = RJ[pj, :]
    RJZ = Matrix(RJph*Zh)
    RJZl = Matrix(RJph*Zlh)
    ZLZ = Matrix(Zth*L*Zh)
    Zrh = indicator(readrows)
    Eah = sparse(auxrows, 1:length(auxrows), ones(length(auxrows)), n, length(auxrows))
    Qh = sparse(transpose(hcat(Zrh, Eah)))
    Minv = ratesystem(G, L, Zrh, Eah).Minv
    d = A -> devicesparse(A, backend)
    return GaussProjection(directions, pj, d(Zh), d(sparse(transpose(Zh))), d(Zth*L), Zth, d(Zlh), d(RJph), d(sparse(transpose(RJph))),
        RJZ, RJZl, lmolj[pj], tobackend(backend, lmolj[pj]),
        isnothing(p.relations) ? emptyrelations(zeros(0)) :
            hostrelations(p.relations, pj),
        ZLZ, d(Zrateh), d(sparse(transpose(Zrateh))),
        Zth*injection, Vector(Zth*constant), Zth*lineinjection, Zth*blockscatter,
        readrows, auxrows, d(Zrh), d(Eah), d(Qh), d(sparse(transpose(Qh))), d(Qh*G), d(Qh*L), Matrix(RJph*Zrh), Minv,
        Qh*injection, Vector(Qh*constant), Qh*lineinjection, Qh*blockscatter)
end

# the Jacobians of the constraints along the projected directions,
# `Z' (L + J'(x)) Z`, one small matrix per column of the projected
# junctions' phases
function projectionmatrices(pr::GaussProjection, hphi::AbstractMatrix)
    return [pr.ZLZ .+ transpose(pr.RJZl)*(Diagonal(pr.lmoljp .*
        derivativeat(pr.relationsp, view(hphi, :, j)))*pr.RJZ) for j in axes(hphi, 2)]
end

# the buffers of the projection over `N` conditions and `m` columns: the
# projected junctions' phases per condition on the backend and the host,
# their products with `m` columns, the coefficients along the directions
# on the host and the backend, and two full columns of work
struct ProjectionWork{A, H, B}
    phip::A
    hphi::H
    phim::A
    zwork::B
    zwork2::B
    g::Matrix{Float64}
    alpha::Matrix{Float64}
    dalpha::B
    work::B
    work2::B
    qwork::B
    qwork2::B
    gq::Matrix{Float64}
    theta::Matrix{Float64}
    thetar::B
    thetaa::B
end

# What a Gauss-Legendre system holds beyond the trapezoidal one: the
# coefficients, the imaginary part of the stage matrix on the Jacobian's
# pattern, the complex matrix the factorization reads, and the projection
# onto the algebraic constraints, or nothing where none needs it.
struct GaussStage{V, M, R, SM, SV}
    coefficients::GaussCoefficients
    imvals::V
    cjacobian::M
    # The projection onto the algebraic constraints, or nothing without
    # any, and the grouped coupling of the rational blocks, or nothing
    # without any, as unions within the backend's sparse matrix and vector
    # types: the stage's type, and with it the system's, the stepper's and
    # every response's, is then the backend's alone, so a circuit with a
    # constraint or a block runs on the code compiled for one without,
    # and a presence check is a branch rather than a specialization.
    projection::Union{Nothing, GaussProjection{SM, SM, SV}}
    # the complex entries of the rational blocks at the stage frequency
    # on the pattern, and the blocks' stage algebra
    rationalvals::R
    rational::Vector{RationalStage}
    coupling::Union{Nothing, RationalCoupling{SM}}
    # The only constructor, and it takes the parameters: `SM` and `SV`
    # appear in the union fields alone, so a circuit with neither a
    # projection nor a coupling passes `nothing` for both and leaves them
    # with nothing to infer from. `gaussstage` reads them off the backend.
    GaussStage{V, M, R, SM, SV}(coefficients, imvals, cjacobian, projection,
            rationalvals, rational, coupling) where {V, M, R, SM, SV} =
        new{V, M, R, SM, SV}(coefficients, imvals, cjacobian, projection,
            rationalvals, rational, coupling)
end

# the stage with its union fields' types taken from the backend
function gaussstage(gc, imvals, cjacobian, projection, rationalvals, stages, coupling, backend)
    SM = typeof(devicesparse(sparse(zeros(1, 1)), backend))
    SV = typeof(tobackend(backend, zeros(1)))
    return GaussStage{typeof(imvals), typeof(cjacobian), typeof(rationalvals), SM, SV}(gc, imvals, cjacobian,
        projection, rationalvals, stages, coupling)
end

function projectionwork(pr::GaussProjection, backend, n::Int, N::Int, m::Int)
    k, npj = length(pr.directions), length(pr.pj)
    kr, ka = length(pr.readrows), length(pr.auxrows)
    kq = kr + ka
    allocate = (dims...) -> KernelAbstractions.zeros(backend, Float64, dims...)
    return ProjectionWork(allocate(npj, N), zeros(npj, N), allocate(npj, m), allocate(k, m), allocate(k, m),
        zeros(k, m), zeros(k, m), allocate(k, m), allocate(n, m), allocate(n, m),
        allocate(kq, m), allocate(kq, m), zeros(kq, m), zeros(kq, m), allocate(kr, m), allocate(ka, m))
end

# The index one unknowns of the endpoint read from their equations: the
# residual `Q (G v + L x + J(x) - b)` from two products on the backend and
# the junction term on the host, the small solve, and the rate along the
# resistive directions and the current on the block rows moved by it.
# `hsin` is the sine of the projected junctions' phases at the endpoint
# (or, for a linearized reading, the cosine times the flux direction's
# phase) on the host, `gq` the drive along `Q` on the host.
function endpointread!(v, x, pr::GaussProjection, pw::ProjectionWork, junction, gq)
    kr = length(pr.readrows)
    isempty(pr.readrows) && isempty(pr.auxrows) && return nothing
    stepmul!(pw.qwork, pr.QG, v)
    stepmul!(pw.qwork2, pr.QL, x)
    pw.qwork .+= pw.qwork2
    copyto!(pw.gq, pw.qwork)
    view(pw.gq, 1:kr, :) .+= transpose(pr.RJQ)*junction
    pw.gq .-= gq
    pw.theta .= -(pr.Minv*pw.gq)
    if kr > 0
        copyto!(pw.thetar, pw.theta[1:kr, :])
        stepmul!(pw.work, pr.Zr, pw.thetar)
        v .+= pw.work
    end
    if !isempty(pr.auxrows)
        copyto!(pw.thetaa, pw.theta[kr + 1:end, :])
        stepmul!(pw.work, pr.Ea, pw.thetaa)
        x .+= pw.work
    end
    return nothing
end

# `Q'` of an injection of targets, on the host, for the reading of a
# direction's current at the endpoint time
function endpointinjection(pr::GaussProjection, injection::SparseMatrixCSC)
    n = size(injection, 1)
    kr, ka = length(pr.readrows), length(pr.auxrows)
    Zrh = sparse(reduce(vcat, pr.readrows; init = Int[]), reduce(vcat, [fill(c, length(z)) for (c, z) in enumerate(pr.readrows)]; init = Int[]),
        ones(sum(length, pr.readrows; init = 0)), n, kr)
    Eah = sparse(pr.auxrows, 1:ka, ones(ka), n, ka)
    return sparse(transpose(hcat(Zrh, Eah)))*injection
end

# The transpose of the reading. The cotangents of the rate along the
# resistive directions and of the current on the block rows,
# `[Zr' vbar; Ea' xbar]`, give through the transposed small solve the
# vector `w = Q' zeta` along the rows, which moves the cotangents of the
# rate and the flux before the reading by `G' w` and `(L + J'(x))' w`,
# with the junction stiffness at the endpoint from the cosines `cosp` of
# the projected junctions' phases on the backend, and the current at the
# endpoint time by `-inj' w`. Leaves `w` in the work buffer and returns
# whether there was anything to do.
function endpointreadtranspose!(vbar, xbar, pr::GaussProjection, pw::ProjectionWork, sys, cosp, nobj, N)
    kr = length(pr.readrows)
    kq = size(pw.gq, 1)
    kq == 0 && return false
    stepmul!(pw.qwork, pr.Q, vbar)
    stepmul!(pw.qwork2, pr.Q, xbar)
    copyto!(pw.gq, pw.qwork)
    copyto!(pw.theta, pw.qwork2)
    view(pw.gq, kr + 1:kq, :) .= view(pw.theta, kr + 1:kq, :)
    pw.theta .= -(transpose(pr.Minv)*pw.gq)
    copyto!(pw.qwork, pw.theta)
    stepmul!(pw.work, pr.Qt, pw.qwork)
    stepmul!(pw.work2, sys.Gt, pw.work)
    vbar .+= pw.work2
    stepmul!(pw.work2, sys.Lt, pw.work)
    xbar .+= pw.work2
    if !isempty(pr.pj)
        stepmul!(pw.phim, pr.RJp, pw.work)
        reshape(pw.phim, :, nobj, N) .*= pr.lmoljpdev .* reshape(cosp, :, 1, N)
        stepmul!(pw.work2, pr.RJpt, pw.phim)
        xbar .+= pw.work2
    end
    return true
end

# the rate along the projected directions replaced by the derivative at
# the endpoint of the cubic through the state, the two stage increments
# and the endpoint: `v += Z R (vrec - v)`, `R` extracting the rate along
# each direction
function projectrate!(v, pr::GaussProjection, pw::ProjectionWork, gc::GaussCoefficients, h, d1, d2, xnew, x)
    w = gc.endrate
    pw.work .= (w[2]/h) .* d1 .+ (w[3]/h) .* d2 .+ (w[4]/h) .* (xnew .- x)
    stepmul!(pw.zwork, pr.Zrate, pw.work)
    stepmul!(pw.zwork2, pr.Zrate, v)
    pw.zwork .-= pw.zwork2
    stepmul!(pw.work, pr.Z, pw.zwork)
    v .+= pw.work
    return v
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
