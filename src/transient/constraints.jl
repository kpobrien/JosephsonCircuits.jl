# The algebraic constraints of a circuit in time: the projection of a
# step's endpoint onto the constraints a junction or a source touches,
# under every rule, the reading of the index one unknowns at a
# Gauss-Legendre endpoint, and the reading of the rate along every
# algebraic direction from the differentiated constraints, which the
# rules apply wherever they report a state and whose tangent reads with
# the constraints perturbed.

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
# on the host for the small Jacobians `Z' (L + J'(x)) Z`, and `Z'` of the
# drive injection and the constant current for the residual. The rate
# along the directions is not the projection's: it is read from the
# differentiated constraint wherever the state is reported
# (`readprojectedrate!`), with the same small Jacobians.
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
struct ConstraintProjection{M, MP, VP}
    # the algebraic directions as their supports, the projected junctions,
    # the directions `Z` as columns and their transpose `Zt`, the
    # constraints' rows times the stiffness `ZcL`, on the host `Zch`, as
    # columns `Zl` and on the backend `Zc`, as `InvariantRate` names them,
    # the junction rows of the incidence and their transpose, the junction
    # phases along the directions and along the constraints, the projected
    # junctions' coefficients, the stiffness along the constraints and
    # directions, and the constraints' rows of the injection, the constant
    # current, the lines' forcing and the blocks' resting waves
    directions::Vector{Vector{Int}}
    pj::Vector{Int}
    Z::M
    Zt::M
    ZcL::M
    Zch::SparseMatrixCSC{Float64,Int}
    Zl::M
    Zc::M
    RJp::MP
    RJpt::MP
    RJpabs::MP
    RJZl::SparseMatrixCSC{Float64,Int}
    lmoljp::Vector{Float64}
    lmoljpdev::VP
    # the current-phase relations of the projected junctions, in the order
    # of `pj`. The projection evaluates them on the host, from `hphi`, and
    # moves the result where a device array needs it, so one host table
    # serves both.
    relationsp::JunctionRelations{Matrix{Float64},Vector{Int}}
    # the pattern of the small Jacobians `Z' (L + J'(x)) Z` with the values
    # of `Z' L Z`, and for each projected junction the positions of the
    # values its stiffness adds to and their coefficients, from `jptr[p]`
    # to `jptr[p + 1] - 1` (see projectionplan)
    jacobian::SparseMatrixCSC{Float64,Int}
    jptr::Vector{Int}
    jidx::Vector{Int}
    jcoef::Vector{Float64}
    Zcinj::SparseMatrixCSC{Float64,Int}
    Zcconstant::Vector{Float64}
    Zcline::SparseMatrixCSC{Float64,Int}
    Zcblock::SparseMatrixCSC{Float64,Int}
    # the entrywise magnitudes of `Z' L`, of `RJZl'` and of the rows of
    # the injection, the lines' forcing and the blocks' scatter, which
    # bound the rounding of the residual, whose terms cancel at the
    # solution; `RJpabs` above is that of `RJp`, which bounds the
    # rounding of the projected junctions' phases
    ZcLabs::M
    RJZltabs::SparseMatrixCSC{Float64,Int}
    Zcinjabs::SparseMatrixCSC{Float64,Int}
    Zclineabs::SparseMatrixCSC{Float64,Int}
    Zcblockabs::SparseMatrixCSC{Float64,Int}
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
    RJQ::SparseMatrixCSC{Float64,Int}
    Minv::RatePseudoinverse
    Qinj::SparseMatrixCSC{Float64,Int}
    Qconst::Vector{Float64}
    Qline::SparseMatrixCSC{Float64,Int}
    Qblock::SparseMatrixCSC{Float64,Int}
    # Where a block is on a direction, its states move the constraint at a
    # rate the differentiated constraint would have to solve for together
    # with the block's port currents, so the rate along the projected
    # directions is then the derivative at the endpoint of the cubic
    # through the state, the stages and the projected endpoint, all on the
    # constraint, third order: `cubic`, and the extraction of the rate
    # along each direction from a rate of the state, `R = (Z' Z)^-1 Z'`,
    # with its transpose.
    cubic::Bool
    Zrate::M
    Zratet::M
end

# the buffers of the projection over `N` conditions and `m` columns: the
# projected junctions' phases per condition on the backend and the host,
# their products with `m` columns, the coefficients along the directions
# on the host and the backend, and two full columns of work
struct ProjectionWork{A, H, B}
    phip::A
    hphi::H
    hphiabs::H
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
    # the conditions whose projection has not reached its floor
    open::Vector{Bool}
    # host work of the residual and of the small solves: the projected
    # junctions' currents and slopes per condition, `(npj, N)`, the
    # floors per condition, the terms along the directions, `(k, m)`, and
    # one condition's Jacobian on the projection's pattern
    current::Matrix{Float64}
    slope::Matrix{Float64}
    floor::Vector{Float64}
    kwork::Matrix{Float64}
    jacobian::SparseMatrixCSC{Float64,Int}
    # the sparse factorizations of the Jacobians and the fill reducing
    # ordering of their pattern, none without directions: one,
    # refactorized for every condition at every solve, or, for the Newton
    # iteration of a stepper's endpoint, one kept per condition (see
    # `projectendpoint!`), nothing until its condition's first
    factors::Vector{Union{Nothing,KLU.KLUFactorization{Float64,Int}}}
    ordering::Union{Nothing,FillOrdering}
    # per condition, for that iteration: the residual before its last
    # correction, whether its next correction refreshes its factorization,
    # and the factorizations it has made; and one condition's slopes
    previous::Vector{Float64}
    refresh::Vector{Bool}
    factorizations::Vector{Int}
    dcol::Vector{Float64}
    # the buffers of the rate system's factorized blocks (see ratework)
    ratework::Vector{NTuple{3,Matrix{Float64}}}
    # the drive of a stepper's endpoint along the constraints' rows and
    # the magnitudes of its terms, `(k, N)`, along the reading's rows,
    # `(kq, N)`, and the magnitudes of the drive values
    gb::Matrix{Float64}
    gbabs::Matrix{Float64}
    gqdrive::Matrix{Float64}
    hostabs::Matrix{Float64}
end

# The algebraic directions, as columns of `Zall` with their constraints
# as rows of `Zcall`, split into those with a junction on their nodes or a
# drive, line or block on them, which the projection handles, and the
# rest, whose linear constraint with a constant right hand side the
# stages keep exactly. The rate along each set is read separately, so a
# direction the stiffness couples to a touched one, through either's
# constraint, is projected with it: the touched set is closed under the
# couplings `Zc L Z`, and a subcircuit nothing couples to the touched
# ones keeps its invariant reading.
function touchedalgebraic(Zall, Zcall, L, RJ, injection, lineinjection, blockscatter, stateful)
    touched, blocked = Int[], Int[]
    # the drives, the lines' forcing and the blocks' scatter along each
    # constraint's row, one column per row
    alongrows = A -> sparse(transpose(Zcall*A))
    drives, lines, onblock = alongrows(injection), alongrows(lineinjection), alongrows(blockscatter)
    # the junctions each direction moves, one column per direction, without
    # those whose two nodes it moves alike
    moved = dropzeros!(RJ*Zall)
    for c in 1:size(Zall, 2)
        junctions = !isempty(nzrange(moved, c))
        block = !isempty(nzrange(onblock, c))
        driven = !isempty(nzrange(drives, c)) || !isempty(nzrange(lines, c)) || block
        (junctions || driven) && push!(touched, c)
        # a block without states, an open, a short or a through, moves
        # no constraint; one with states does
        any(q -> stateful[rowvals(onblock)[q]] && nonzeros(onblock)[q] != 0, nzrange(onblock, c)) && push!(blocked, c)
    end
    if !isempty(touched) && length(touched) < size(Zall, 2)
        # the union of the couplings' pattern with its transpose's: the
        # entries' signs follow the bases of the directions and of the
        # constraints, which are chosen apart, and a coupling would cancel
        # against its transpose where they differ
        coupling = abs.(droptol!(Zcall*L*Zall, 1e-14))
        coupling = coupling + sparse(transpose(coupling))
        reached = falses(size(Zall, 2))
        reached[touched] .= true
        queue = copy(touched)
        while !isempty(queue)
            c = popfirst!(queue)
            for d in rowvals(coupling)[nzrange(coupling, c)]
                reached[d] && continue
                reached[d] = true
                push!(queue, d)
            end
        end
        touched = findall(reached)
    end
    return touched, setdiff(1:size(Zall, 2), touched), blocked
end

# whether each column of the blocks' scatter, a block port, belongs to a
# block with states
function statefulcolumns(p::TransientProblem, blockscatter)
    cols = [size(b.A, 1) > 0 for b in p.blocks for q in eachindex(b.signal)]
    return length(cols) == size(blockscatter, 2) ? cols : trues(size(blockscatter, 2))
end

# The partition of a problem's algebraic directions, computed once for
# the projection and the invariant reading to share: the directions and
# the constraints with their negligible entries dropped, and the
# touched, untouched and blocked ones of `touchedalgebraic`, so that no
# coupling the stiffness keeps crosses the two readings.
function algebraicpartition(p::TransientProblem, L::SparseMatrixCSC, RJ::SparseMatrixCSC, injection::SparseMatrixCSC,
        lineinjection::SparseMatrixCSC, blockscatter::SparseMatrixCSC)
    Zall, Zcall = sparse(p.directions), sparse(p.constraints)
    droptol!(Zall, 1e-14); droptol!(Zcall, 1e-14)
    touched, untouched, blocked = touchedalgebraic(Zall, Zcall, L, RJ, injection, lineinjection, blockscatter,
        statefulcolumns(p, blockscatter))
    return (; Zall, Zcall, touched, untouched, blocked)
end

# The rate along the algebraic directions the projection leaves alone.
# Their constraint `Z' L x` is an invariant of the stages, and its rate
# `Z' L v = 0` of the rule in exact arithmetic; the rule's rate update
# carries the rounding of every stage solve along such a direction
# forward, since nothing damps it there, so wherever the state is
# reported the rate along these directions is read from the
# differentiated constraint, `v - Z (Zc L Z)^-1 Zc L v`. The directions
# `Z` and their transpose, their constraints' rows `Zc`, those rows as
# columns `Zl` and the rows times the stiffness `Zc L` are held on the
# backend, and `Zc L Z`, which has the sparsity of the couplings between
# the constraints, factorized once on the host.
struct InvariantRate{M, F}
    Z::M
    Zt::M
    Zc::M
    Zl::M
    ZcL::M
    factor::F
end

function invariantrate(p::TransientProblem, L::SparseMatrixCSC, partition, backend)
    untouched = partition.untouched
    isempty(untouched) && return nothing
    Zh, Zch = partition.Zall[:, untouched], partition.Zcall[untouched, :]
    d = A -> devicesparse(A, backend)
    return InvariantRate(d(Zh), d(sparse(transpose(Zh))), d(Zch), d(sparse(transpose(Zch))), d(Zch*L), lu(Zch*L*Zh))
end

# The work of the reading over `m` columns: the coefficients along the
# directions on the backend, a second set for a tangent's forcing, on a
# device their copy on the host, the solution on the host, and a column
# of the state per column read; and the factorization's own handle, a
# copy sharing the factors but with its own solve workspace and lock.
# One per chunk of a batch, since a chunk reads while another does; on a
# device the reading is a round trip to the host.
struct InvariantRateWork{A, B, F}
    zwork::A
    zforce::A
    zhost::Matrix{Float64}
    zsol::Matrix{Float64}
    work::B
    factor::F
end

function invariantratework(ir::InvariantRate, backend, n::Integer, m::Integer)
    k = size(ir.Zc, 1)
    zwork = KernelAbstractions.zeros(backend, Float64, k, m)
    zhost = backend isa CPU ? zwork : zeros(k, m)
    return InvariantRateWork(zwork, KernelAbstractions.zeros(backend, Float64, k, m), zhost, zeros(k, m),
        KernelAbstractions.zeros(backend, Float64, n, m), copy(ir.factor))
end

# `dst` set to the columns of `v` with their rate along the invariant
# directions read from the differentiated constraint. A tangent along a
# component reads its final rate through the same map with the derivative
# of the constraint's row added: for the stiffness perturbed by `dL`, the
# read rate `w` of the nominal solve keeps `Z' (L + dL) w = 0`, so the
# tangent of the reading carries `-Z (Z' L Z)^-1 Z' dL w` besides the
# reading of the tangent rate itself; `forcing` is then `dL w` per column.
function readinvariantrate!(dst, v, ir::InvariantRate, w::InvariantRateWork, forcing = nothing)
    stepmul!(w.zwork, ir.ZcL, v)
    if !isnothing(forcing)
        stepmul!(w.zforce, ir.Zc, forcing)
        w.zwork .+= w.zforce
    end
    w.zhost === w.zwork || copyto!(w.zhost, w.zwork)
    ldiv!(w.zsol, w.factor, w.zhost)
    copyto!(w.zwork, w.zsol)
    stepmul!(w.work, ir.Z, w.zwork)
    w.work .= v .- w.work
    copyto!(dst, w.work)
    return dst
end

# The rate along the projected directions read from the differentiated
# constraint at the state `x`: `dst = v + Z alpha` with
# `Z' (L + J'(x)) Z alpha = Z' bdot - Z' (L + J'(x)) v`, the Jacobians of
# the projection per condition, the junction term from the projected
# junctions' phases on the host, and `bdotz` the rate of the drive along
# the constraints' rows, `(directions, columns)` on the host. The columns
# of `v` are the directions of every condition, contiguous per condition,
# `x` one column per condition, or `nothing` where the projected
# junctions' phases are already in the work's `hphi`. A tangent's reading
# adds `Z' forcing`,
# with `forcing` the perturbation of the constraint's operator on the
# solve's read rate, `(state, columns)` on the backend, and `curvature`,
# a term formed on the host along the rows.
function readprojectedrate!(dst, v, x, pr::ConstraintProjection, pw::ProjectionWork, bdotz;
        forcing = nothing, curvature = nothing)
    dst === v || copyto!(dst, v)
    (isempty(pr.directions) || pr.cubic) && return dst
    if !isempty(pr.pj) && !isnothing(x)
        stepmul!(pw.phip, pr.RJp, x)
        copyto!(pw.hphi, pw.phip)
    end
    constrainttangent!(pw, pr, v, forcing)
    pw.g .-= bdotz
    isnothing(curvature) || (pw.g .+= curvature)
    constraintcorrection!(pw, pr)
    dst .+= pw.work
    return dst
end

# The correction along the directions from the violation in the work's
# `g`: the coefficients `alpha = -(Z' (L + J'(x)) Z)^-1 g` per condition,
# with the Jacobians at the phases in `hphi`, or with `kept` on the
# factorizations the conditions keep (see `projectionsolve!`), and the
# state's move `Z alpha` into the work's `work`; none for the conditions
# `open` does not hold, when it is given.
function constraintcorrection!(pw::ProjectionWork, pr::ConstraintProjection, open = nothing; kept::Bool = false)
    projectionsolve!(pw, pr, false, open; kept)
    pw.alpha .*= -1
    copyto!(pw.dalpha, pw.alpha)
    stepmul!(pw.work, pr.Z, pw.dalpha)
    return pw.work
end

# The endpoint `x`, one column per condition, projected onto the
# constraints: Newton on the coefficients along the directions, from the
# residual and its floor of `constraintresidual!` and the correction of
# `constraintcorrection!`, until the residual is at its floor, `gb` being
# the drive along the constraints' rows, `(directions, conditions)` on
# the host, and `gbabs` the magnitudes of its terms; a linear constraint
# is met by one correction. Each condition keeps its factorization of
# the Jacobian across corrections and endpoints, as the stage operator
# keeps its own (see `newtonsolve!`), a correction on it a chord
# (simplified Newton) step, and refreshes it at the current state before
# a correction when the last left the residual above `contraction` times
# the one before, or when at the rate it contracted the floor lies beyond
# the corrections left, the floor being many orders below an endpoint's
# first residual, or when it has none of its own yet. A refresh is the
# Jacobian at that state, so poorly contracting corrections refresh as
# often as Newton does, and it is not carried to the next endpoint,
# whose residual always asks for a correction. The work holds one
# factorization per condition, or one for one condition. A condition is
# corrected until its own residual has reached its floor and no further,
# so that its endpoint does not depend on the conditions beside it, and
# the work's `open` holds the conditions which have not. `moved!` is
# called with each move of the state. Returns whether every residual
# reached its floor within `iterations` corrections.
function projectendpoint!(x, pw::ProjectionWork, pr::ConstraintProjection, gb, gbabs, iterations::Int, moved! = nothing;
        contraction::Float64 = NEWTONCONTRACTION)
    open = pw.open
    fill!(pw.previous, Inf)
    for iteration in 1:iterations + 1
        floor = constraintresidual!(pw, pr, x, gb, gbabs)
        for j in axes(gb, 2)
            residual = maximum(abs, view(pw.g, :, j))
            open[j] = !(residual <= floor[j])
            rate = residual/pw.previous[j]
            slow = 0 < rate < 1 && floor[j] < residual && log(floor[j]/residual)/log(rate) > iterations - iteration + 1
            pw.refresh[j] = pw.factorizations[j] == 0 || rate > contraction || slow
            pw.previous[j] = residual
        end
        any(open) || return true
        iteration > iterations && return false
        constraintcorrection!(pw, pr, open; kept = true)
        x .+= pw.work
        isnothing(moved!) || moved!(pw.work)
        isempty(pr.pj) && return true
    end
    return false
end

# The violation of the constraints along the projected directions at the
# state `x`, one column per condition: `Z' (L x + J(x)) - gb` with `gb`
# the drive along the rows, `(directions, conditions)` on the host, into
# the work's `g`, the projected junctions' phases left in `hphi`. Returns
# the roundoff floor of the violation per condition, the unit roundoff
# times the magnitudes of the terms summed, entry by entry, `gbabs` being
# those of the drive: the terms cancel at the solution, two junctions in
# series on an island balancing their currents, so their sum says
# nothing of the rounding in it. A junction's phase is a difference of
# node fluxes which may be far larger than it, in a series array or in
# the voltage state, so the rounding of the phases, `|RJp| |x|`, enters
# through the junctions' stiffness as well.
function constraintresidual!(pw::ProjectionWork, pr::ConstraintProjection, x, gb, gbabs)
    if !isempty(pr.pj)
        stepmul!(pw.phip, pr.RJp, x)
        copyto!(pw.hphi, pw.phip)
    end
    stepmul!(pw.zwork, pr.ZcL, x)
    copyto!(pw.g, pw.zwork)
    pw.work .= abs.(x)
    stepmul!(pw.zwork2, pr.ZcLabs, pw.work)
    copyto!(pw.alpha, pw.zwork2)
    current, slope = pw.current, pw.slope
    relationinto!(current, pr.relationsp, pw.hphi)
    current .*= pr.lmoljp
    slope .= abs.(current)
    mul!(pw.kwork, pr.RJZltabs, slope)
    pw.alpha .+= pw.kwork .+ gbabs
    if !isempty(pr.pj)
        stepmul!(pw.phim, pr.RJpabs, pw.work)
        copyto!(pw.hphiabs, pw.phim)
        derivativeinto!(slope, pr.relationsp, pw.hphi)
        slope .= pr.lmoljp .* abs.(slope) .* pw.hphiabs
        mul!(pw.kwork, pr.RJZltabs, slope)
        pw.alpha .+= pw.kwork
    end
    for j in axes(gb, 2)
        pw.floor[j] = 8eps(Float64)*maximum(view(pw.alpha, :, j))
    end
    mul!(pw.g, transpose(pr.RJZl), current, 1.0, 1.0)
    pw.g .-= gb
    return pw.floor
end

# The violation linearized at the state whose projected junctions' phases
# the work's `hphi` holds: `Z' (L + J'(x)) dx` per column, with
# `Z' forcing` added for a `forcing` `(state, columns)` on the backend,
# into the work's `g`; the columns are the directions of every condition,
# contiguous per condition.
function constrainttangent!(pw::ProjectionWork, pr::ConstraintProjection, dx, forcing = nothing)
    N = size(pw.hphi, 2)
    m = size(pw.g, 2)
    dirs = m ÷ N
    stepmul!(pw.zwork, pr.ZcL, dx)
    if !isnothing(forcing)
        stepmul!(pw.zwork2, pr.Zc, forcing)
        pw.zwork .+= pw.zwork2
    end
    copyto!(pw.g, pw.zwork)
    if !isempty(pr.pj)
        stepmul!(pw.phim, pr.RJp, dx)
        hphim = Array(pw.phim)
        for j in 1:N
            cols = (j - 1)*dirs + 1:j*dirs
            view(pw.g, :, cols) .+= transpose(pr.RJZl)*(pr.lmoljp .*
                derivativeat(pr.relationsp, view(pw.hphi, :, j)) .* view(hphim, :, cols))
        end
    end
    return pw.g
end

# The projection's tangent transposed, at the phases in the work's
# `hphi`: from the cotangent `xbar` of the projected flux, `(state,
# columns)` on the backend, the multiplier `gamma = M'^-1 Z' xbar` along
# the constraints, left in the work's `alpha`, and its image `Zl gamma` in
# the work's `work`, the cotangent of the flux before the projection,
# `xbar -= (L + J'(x))' Zl gamma`, and the cotangent of the projection's
# forcing, `fbar = -Zl gamma`, for the currents and the components the
# caller contracts it with. `cosp`, the derivative of the projected
# junctions' relations at those phases on the backend, is formed here
# where the caller has not formed it.
function constrainttranspose!(xbar, fbar, pw::ProjectionWork, pr::ConstraintProjection, sys, backend; cosp = nothing)
    N = size(pw.hphi, 2)
    m = size(pw.g, 2)
    dirs = m ÷ N
    stepmul!(pw.zwork, pr.Zt, xbar)
    copyto!(pw.g, pw.zwork)
    projectionsolve!(pw, pr, true)
    copyto!(pw.dalpha, pw.alpha)
    stepmul!(pw.work, pr.Zl, pw.dalpha)
    fbar .= .-pw.work
    stepmul!(pw.work2, sys.Lt, pw.work)
    xbar .-= pw.work2
    if !isempty(pr.pj)
        stepmul!(pw.phim, pr.RJp, pw.work)
        c = isnothing(cosp) ? tobackend(backend, derivativeat(pr.relationsp, pw.hphi)) : cosp
        reshape(pw.phim, :, dirs, N) .*= pr.lmoljpdev .* reshape(c, :, 1, N)
        stepmul!(pw.work2, pr.RJpt, pw.phim)
        xbar .-= pw.work2
    end
    return fbar
end

# `Z'` of the rate of the drive at `t` for every condition, into `bdotz`,
# `(directions, conditions)` on the host: the drives by the difference of
# `ratestencil` in a record starting at `t0`, the lines' forced currents
# by `linerate`, a difference of their histories, and the blocks' resting
# waves held. `hv1` and `hv2` are `(drives, conditions)` host work.
function drivedotz!(bdotz, pr::ConstraintProjection, problems, t, t0, delta, hv1, hv2, linerate = nothing)
    offsets, weights = ratestencil(t, t0, delta)
    fill!(hv1, 0)
    for k in eachindex(offsets)
        iszero(weights[k]) && continue
        drivevalues!(hv2, problems, t + offsets[k]*delta)
        hv1 .+= weights[k] .* hv2
    end
    hv1 ./= 12delta
    mul!(bdotz, pr.Zcinj, hv1)
    isnothing(linerate) || (bdotz .+= pr.Zcline*linerate)
    return bdotz
end

# The five point differences which read the rate of the drives at the step
# `delta` of `ratedelta`, as the offsets of their samples in units of
# `delta` and their weights, whose sum over `12 delta` is the rate: the
# central one, and the one sided one on `[t, t + 4 delta]`, fourth order
# both, the one sided one with seven times the rounding.
const CENTRALRATE = ((-2.0, -1.0, 1.0, 2.0, 0.0), (1.0, -8.0, 8.0, -1.0, 0.0))
const ONESIDEDRATE = ((0.0, 1.0, 2.0, 3.0, 4.0), (-25.0, 48.0, -36.0, 16.0, -3.0))
# the difference a record starting at `t0` reads the rate of its drives at
# `t` with: the central one within the record, and the one sided one at its
# start, since the integration reads no drive before the start, and one
# which begins there, `t <= t0 ? 0.0 : ...`, is smooth on that side alone
ratestencil(t, t0, delta) = t - 2delta < t0 ? ONESIDEDRATE : CENTRALRATE

# the step of the difference which reads the rate of the drives, a
# hundredth of the step: the five point difference's truncation at the
# fastest drive a step resolves and its rounding both sit near the
# roundoff of the rate
ratedelta(sys) = 1e-2*sys.h

# The work of the readings over `N` conditions and `m` columns: the
# invariant reading's and the projection's, each nothing without the
# constraints it serves, as unions within one concrete type so that a
# presence check is a branch rather than a specialization; the rows of
# the drive's rate and two sets of drive values on the host; and two
# columns of the state per column for the transpose.
struct RateReadWork{A, F}
    iw::Union{Nothing, InvariantRateWork{A, A, F}}
    pw::Union{Nothing, ProjectionWork{A, Matrix{Float64}, A}}
    bdotz::Matrix{Float64}
    hv1::Matrix{Float64}
    hv2::Matrix{Float64}
    u::A
    t::A
    # the only constructor, and it takes the parameters: `F` appears in
    # the invariant reading's union alone, so a system without one leaves
    # it with nothing to infer from
    RateReadWork{A, F}(iw, pw, bdotz, hv1, hv2, u, t) where {A, F} = new{A, F}(iw, pw, bdotz, hv1, hv2, u, t)
end

function ratereadwork(sys, backend, n, N, m)
    ir, pr = sys.invariant, sys.projection
    nd = length(sys.problem.drives)
    k = isnothing(pr) ? 0 : length(pr.directions)
    u = KernelAbstractions.zeros(backend, Float64, n, m)
    return RateReadWork{typeof(u), SparseArrays.UMFPACK.UmfpackLU{Float64, Int}}(
        isnothing(ir) ? nothing : invariantratework(ir, backend, n, m),
        isnothing(pr) ? nothing : projectionwork(pr, backend, n, N, m),
        zeros(k, m), zeros(nd, N), zeros(nd, N), u, KernelAbstractions.zeros(backend, Float64, n, m))
end

# The reading of the rate along every algebraic direction at the state
# `x`, into `dst`: the invariant constraints through `readinvariantrate!`,
# the projected ones through `readprojectedrate!` with the rows of the
# drive's rate in the work, whichever the system has. `v`, `x` and `dst`
# are `(state, columns)`.
function readrate!(dst, v, x, sys, rw::RateReadWork; forcing = nothing, curvature = nothing)
    ir, pr = sys.invariant, sys.projection
    iw, pw = rw.iw, rw.pw
    if isnothing(ir) || isnothing(iw)
        dst === v || copyto!(dst, v)
    else
        readinvariantrate!(dst, v, ir, iw, forcing)
    end
    (isnothing(pr) || isnothing(pw)) || readprojectedrate!(dst, dst, x, pr, pw, rw.bdotz; forcing, curvature)
    return dst
end

# The reading of a tangent's final rate as the solve's, `(state, columns)`
# over the directions of every condition. Along a perturbation of the
# components the constraints' rows move with it: the stiffness by
# `dL w`, the junction stiffness by the junctions' own perturbation and
# by its curvature along the tangent flux, all on the read final rate
# `w` of the solve. Along a tangent of the currents the rate of the
# currents at the end of the grid is the drive's rate, the derivative of
# the cubic through the last four grid values; `currents` are the
# tangent's on the host (see `TangentCurrents`), one zero column of them
# for a tangent along the components alone, and `injection` the targets'
# injection, scaled, on the host; `linerate`, where the circuit has lines,
# the rate of the currents the tangent's arriving waves force at the end,
# `(line port, columns)` on the backend. Both perturbations of the
# constraints enter as one forcing, `dL w - dbdot`, so a direction the
# nominal problem leaves undriven reads the rate its tangent current
# gives it.
function readtangentrate!(finalrate, dv, dx, sys, sol, N, ndir, perturbation, currents, injection, backend, linerate = nothing)
    n = size(dv, 1)
    o = outputreading(sys, backend, n, N, ndir*N)
    xh = Array(reshape(sol.finalflux, n, N))
    wh = Array(reshape(sol.finalrate, n, N))
    pr = sys.projection
    pw = o.rn.pw
    if !isnothing(pw)
        xbar, wbar = tobackend(backend, copy(xh)), tobackend(backend, copy(wh))
        stepmul!(pw.phip, pr.RJp, xbar)
        copyto!(o.hphi, pw.phip)
        stepmul!(pw.phip, pr.RJp, wbar)
        copyto!(o.er, pw.phip)
    end
    states = !isnothing(perturbation) && perturbation.states
    readoutputs!(o, sys, size(currents.grid, 2), dv, dx, currents, injection, states ? xh : zeros(0, 0), states ? wh : zeros(0, 0),
        perturbation, backend, linerate)
    copyto!(finalrate, o.dvread)
    return finalrate
end

# The work of a tangent's reading of a port on an algebraic direction at
# every time, one concrete object: the readings' work over the columns
# and over the conditions, the read tangent rate, the read rate of the
# solve at the time, the forcing on the host and the backend with the
# lines' part of it on the backend, the curvature term, and the projected
# junctions' phases and read rates on the host with a column of them on
# the backend.
struct OutputReading{A, F}
    rr::RateReadWork{A, F}
    rn::RateReadWork{A, F}
    dvread::A
    wk::A
    gh::Matrix{Float64}
    gdev::A
    lineforcing::A
    ch::Matrix{Float64}
    hphi::Matrix{Float64}
    er::Matrix{Float64}
    erdev::A
end

function outputreading(sys, backend, n, N, m)
    pr = sys.projection
    npj = isnothing(pr) ? 0 : length(pr.pj)
    k = isnothing(pr) ? 0 : length(pr.directions)
    z = (dims...) -> KernelAbstractions.zeros(backend, Float64, dims...)
    return OutputReading(ratereadwork(sys, backend, n, N, m), ratereadwork(sys, backend, n, N, N),
        z(n, m), z(n, N), zeros(n, m), z(n, m), z(n, m), zeros(k, m), zeros(npj, N), zeros(npj, N), z(npj, N))
end

# the tangent rate at grid point `k` read for the port waves, into the
# work's `dvread`: the reading linearized at the projected junctions'
# phases and read rates the work holds, with the tangent currents' rate
# through the stencil, where the circuit has lines the rate `linerate` of
# the currents the tangent's arriving waves force, `(line port, columns)`
# on the backend, which moves the constraints as the drives' rate does,
# and, where the perturbation reads the states, the states `xh` and their
# read rate `wh` on the host
function readoutputs!(o::OutputReading, sys, k, dv, dx, currents, injection, xh, wh, perturbation, backend, linerate = nothing)
    tangentforcing!(o.gh, o.ch, sys, k, currents, injection, o.hphi, o.er, dx, xh, wh, perturbation, o.rr.pw)
    copyto!(o.gdev, o.gh)
    if !isnothing(linerate)
        stepmul!(o.lineforcing, sys.lineinjection, linerate)
        o.gdev .-= o.lineforcing
    end
    isnothing(o.rr.pw) || copyto!(o.rr.pw.hphi, o.hphi)
    return readrate!(o.dvread, dv, nothing, sys, o.rr; forcing = o.gdev, curvature = o.ch)
end

# The work of an adjoint's transpose of that reading: the readings' work,
# the cotangent of the read rate and of the forcing, the read rate of the
# solve at the time, the projected junctions' phases and read rates, and
# two columns over the junctions for the components' contraction.
struct OutputTranspose{A, F}
    rr::RateReadWork{A, F}
    rn::RateReadWork{A, F}
    rbar::A
    fbar::A
    wk::A
    hphi::Matrix{Float64}
    er::Matrix{Float64}
    erdev::A
    ystate::A
    yphi::A
end

function outputtranspose(sys, backend, n, N, m)
    pr = sys.projection
    npj = isnothing(pr) ? 0 : length(pr.pj)
    nj = length(sys.lmolj)
    z = (dims...) -> KernelAbstractions.zeros(backend, Float64, dims...)
    return OutputTranspose(ratereadwork(sys, backend, n, N, m), ratereadwork(sys, backend, n, N, N),
        z(n, m), z(n, m), z(n, N), zeros(npj, N), zeros(npj, N), z(npj, N), z(nj, N), z(nj, N))
end

# The transpose of the reading of a port on an algebraic direction at
# grid point `k` of `nt`, from the cotangent of the read rate in the
# work's `rbar`: the cotangent of the rate before the reading into
# `vbar`, of the flux into `xbar`, of the currents into the ring at the
# points of the stencil, and of the components into `sens` through the
# entries' contraction, on the read rate the work holds and, where the
# perturbation reads the states, the flux `xk`.
function readoutputstranspose!(o::OutputTranspose, sys, k, nt, vbar, xbar, ring, injectiont, targetwork,
        perturbation, pcwork, sens, xk, backend)
    readtranspose!(vbar, xbar, o.fbar, o.rbar, sys, o.rr, o.hphi, o.er, backend)
    idx, wgt = currentratestencil(k, nt)
    if !isempty(idx)
        stepmul!(targetwork, injectiont, o.fbar)
        for i in eachindex(idx)
            ringadd!(ring, idx[i], targetwork, -wgt[i]/sys.h)
        end
    end
    if !isnothing(perturbation)
        pr = sys.projection
        if perturbation.states
            contractentries!(sens, perturbation.entries.L, pcwork.L, o.wk, o.fbar, -1.0)
            stepmul!(o.yphi, sys.RJ, xk)
            derivativeinto!(o.ystate, sys.relations, o.yphi)
            stepmul!(o.yphi, sys.RJ, o.wk)
            o.ystate .*= o.yphi
        else
            yh = zeros(size(o.ystate))
            (isnothing(pr) || isempty(pr.pj)) || (yh[pr.pj, :] .= derivativeat(pr.relationsp, o.hphi) .* o.er)
            copyto!(o.ystate, yh)
        end
        contractentries!(sens, perturbation.entries.J, pcwork.J, o.ystate, o.fbar, -1.0)
    end
    return nothing
end

# The derivative at grid point `k` of `nt` of the polynomial through the
# grid values the stencil reads the currents by, as the points and their
# weights per step: the cubic through the point and the three before it,
# and through the first four points at the first three, so that the
# transpose reaches no column an adjoint has emitted; on a grid of three
# points the quadratic through them, and on one of two the line. Empty
# on a single point, the one zero column of a tangent along the
# components alone.
function currentratestencil(k::Int, nt::Int)
    nt < 2 && return (Int[], Float64[])
    nt == 2 && return ([1, 2], [-1.0, 1.0])
    if nt == 3
        k == 1 && return ([1, 2, 3], [-3.0, 4.0, -1.0] ./ 2)
        k == 2 && return ([1, 2, 3], [-1.0, 0.0, 1.0] ./ 2)
        return ([1, 2, 3], [1.0, -4.0, 3.0] ./ 2)
    end
    k == 1 && return ([1, 2, 3, 4], [-11.0, 18.0, -9.0, 2.0] ./ 6)
    k == 2 && return ([1, 2, 3, 4], [-2.0, -3.0, 6.0, -1.0] ./ 6)
    k == 3 && return ([1, 2, 3, 4], [1.0, -6.0, 3.0, 2.0] ./ 6)
    return ([k - 3, k - 2, k - 1, k], [-2.0, 9.0, -18.0, 11.0] ./ 6)
end

# The forcing of a tangent's reading at grid point `k`, `(state, columns)`
# on the host into `g`, and the curvature term along the projected
# directions into `c`, `(directions, columns)`, or `nothing`: the rate
# of the tangent's `currents` on the grid (see `TangentCurrents`),
# through the stencil and the targets' `injection`; the
# perturbation of the stiffness and of the junctions on the read rate
# `wh` at the state `xh`, `(state, conditions)`, where the perturbation
# reads the states, and on the projected junctions' read rate `er` at
# their phases `hphi`, `(junctions, conditions)`, otherwise; and the
# curvature of the projected junctions' stiffness along the tangent flux
# `dx`, `(state, columns)` on the backend, on `er`. Columns are the
# directions of every condition, contiguous per condition.
function tangentforcing!(g, c, sys, k, currents, injection, hphi, er, dx, xh, wh, perturbation, pw)
    pr = sys.projection
    n, m = size(g)
    N = size(er, 2)
    fill!(g, 0)
    idx, wgt = currentratestencil(k, size(currents.grid, 2))
    if !isempty(idx)
        work = isempty(currents.targets) ? nothing : zeros(size(injection, 2), directions(currents))
        dIdot = sum(wgt[i] .* gridcurrents!(work, currents, idx[i]) for i in eachindex(idx)) ./ sys.h
        g .-= repeat(injection*dIdot, 1, N)
    end
    if !isnothing(perturbation) && !isnothing(perturbation.forcing)
        f = perturbation.forcing
        if perturbation.states
            RJh = sys.problem.RJ
            hr = hostrelations(sys.relations)
            g .+= reshape(f.hdL*wh .+ f.hdJ*(derivativeat(hr, RJh*xh) .* (RJh*wh)), n, m)
        elseif !isnothing(pr) && !isempty(pr.pj)
            y = zeros(length(sys.lmolj), N)
            y[pr.pj, :] .= derivativeat(pr.relationsp, hphi) .* er
            g .+= reshape(f.hdJ*y, n, m)
        end
    end
    if !isnothing(pr) && !isempty(pr.pj) && !pr.cubic
        dirs = m ÷ N
        stepmul!(pw.phim, pr.RJp, dx)
        hphim = Array(pw.phim)
        ns = negsecondat(pr.relationsp, hphi)
        for j in 1:N, d in 1:dirs
            col = (j - 1)*dirs + d
            c[:, col] .= -transpose(pr.RJZl)*(pr.lmoljp .* view(ns, :, j) .* view(hphim, :, col) .* view(er, :, j))
        end
    end
    return g
end

# The transpose of the reading at a time. `rbar`, `(state, columns)` on
# the backend, is the cotangent of the read rate; the cotangent of the
# rate before the reading is added to `vbar`, the cotangent of the flux
# through the curvature of the projected junctions' stiffness to `xbar`,
# and the cotangent of the forcing is left in `fbar`, for the currents
# and the components the caller contracts it with. `hphi` and `er` are
# the projected junctions' phases and read rates, `(junctions,
# conditions)` on the host; the columns are the objectives of every
# condition, contiguous per condition. The projected reading is applied
# after the invariant one, so its transpose comes first.
function readtranspose!(vbar, xbar, fbar, rbar, sys, rw, hphi, er, backend)
    ir, pr = sys.invariant, sys.projection
    pw, iw, u, t = rw.pw, rw.iw, rw.u, rw.t
    m = size(rbar, 2)
    N = size(er, 2)
    copyto!(u, rbar)
    fill!(fbar, 0)
    if !isnothing(pr) && !isempty(pr.directions) && !pr.cubic
        dirs = m ÷ N
        copyto!(pw.hphi, hphi)
        constrainttranspose!(u, fbar, pw, pr, sys, backend)
        if !isempty(pr.pj)
            # the curvature transposed: `-RJp' (lmolj ns er .* RJZl cbar)`
            # with `cbar = -gamma`
            q = pr.RJZl*pw.alpha
            ns = negsecondat(pr.relationsp, hphi)
            hq = reshape(q, :, dirs, N) .* reshape(pr.lmoljp .* ns .* er, :, 1, N)
            copyto!(pw.phim, reshape(hq, :, m))
            stepmul!(pw.work2, pr.RJpt, pw.phim)
            xbar .+= pw.work2
        end
    end
    if !isnothing(ir)
        stepmul!(iw.zwork, ir.Zt, u)
        iw.zhost === iw.zwork || copyto!(iw.zhost, iw.zwork)
        ldiv!(iw.zsol, transpose(iw.factor), iw.zhost)
        copyto!(iw.zwork, iw.zsol)
        stepmul!(iw.work, ir.Zl, iw.zwork)
        fbar .-= iw.work
        stepmul!(t, sys.Lt, iw.work)
        u .-= t
    end
    vbar .+= u
    return fbar
end

# The projection of a problem under the scaled matrices, or nothing
# where no direction needs one. A direction whose constraint no junction,
# drive, line or block enters is linear with a constant right hand side,
# which the stages keep exactly, so it is left alone; the others are
# projected. The reading covers every capacitor free island and block
# current row through the rate system's pseudoinverse, which determines
# what the equations determine and leaves the algebraic directions to
# the projection.
function constraintprojection(p::TransientProblem, G::SparseMatrixCSC, L::SparseMatrixCSC, RJ::SparseMatrixCSC,
        lmolj::Vector{Float64}, injection::SparseMatrixCSC, constant::Vector{Float64}, lineinjection::SparseMatrixCSC,
        blockscatter::SparseMatrixCSC, partition, backend)
    n = length(p)
    Zall, Zcall, touched, blocked = partition.Zall, partition.Zcall, partition.touched, partition.blocked
    directions = p.algebraic[touched]
    Zh, Zch = Zall[:, touched], Zcall[touched, :]
    Zlh = sparse(transpose(Zch))
    Rh = isempty(touched) ? spzeros(0, n) : directionreading(Zh)
    # the reading: the capacitor free islands and the block current rows
    readrows = [z for z in p.inertialess if all(<=(p.Nnodal), z)]
    auxrows = [b.auxbase + q for b in p.blocks for q in eachindex(b.signal)]
    (isempty(directions) && isempty(readrows) && isempty(auxrows)) && return nothing
    indicator = dirs -> sparse(foldl(append!, dirs; init = Int[]), [c for (c, z) in enumerate(dirs) for _ in z],
        ones(sum(length, dirs; init = 0)), n, length(dirs))
    support = [i for c in 1:size(Zh, 2) for i in findnz(Zh[:, c])[1]]
    pj = sort!(unique!([i for node in vcat(support, foldl(append!, readrows; init = Int[])) for i in rowvals(RJ)[nzrange(RJ, node)]]))
    RJph = RJ[pj, :]
    RJZl = RJph*Zlh
    jacobian, jptr, jidx, jcoef = projectionplan(Zch*L*Zh, RJZl, RJph*Zh)
    Zrh = indicator(readrows)
    Eah = sparse(auxrows, 1:length(auxrows), ones(length(auxrows)), n, length(auxrows))
    Qh = sparse(transpose(hcat(Zrh, Eah)))
    Minv = ratesystem(G, L, Zrh, Eah).pinv
    d = A -> devicesparse(A, backend)
    return ConstraintProjection(directions, pj, d(Zh), d(sparse(transpose(Zh))), d(Zch*L), Zch, d(Zlh), d(Zch), d(RJph), d(sparse(transpose(RJph))),
        d(abs.(RJph)),
        RJZl, lmolj[pj], tobackend(backend, lmolj[pj]),
        isnothing(p.relations) ? emptyrelations(zeros(0)) :
            hostrelations(p.relations, pj),
        jacobian, jptr, jidx, jcoef, Zch*injection, Vector(Zch*constant), Zch*lineinjection, Zch*blockscatter,
        d(abs.(Zch*L)), abs.(sparse(transpose(RJZl))), abs.(Zch*injection), abs.(Zch*lineinjection), abs.(Zch*blockscatter),
        readrows, auxrows, d(Zrh), d(Eah), d(Qh), d(sparse(transpose(Qh))), d(Qh*G), d(Qh*L), RJph*Zrh, Minv,
        Qh*injection, Vector(Qh*constant), Qh*lineinjection, Qh*blockscatter,
        !isempty(blocked), d(Rh), d(sparse(transpose(Rh))))
end

# `(Z' Z)^-1 Z'`, which reads the coefficients along the directions `Z`
# of a vector of the state, sparse: `Z' Z` couples only the directions
# which share a node, so it is solved block by block of those couplings,
# each on the nodes its directions span
function directionreading(Z::SparseMatrixCSC)
    n, k = size(Z)
    Zt = sparse(transpose(Z))
    ZtZ = Zt*Z
    rr, rc, rv = Int[], Int[], Float64[]
    at = zeros(Int, k)
    for (D, _) in blockcomponents(ZtZ)
        nodes = sort!(unique!([rowvals(Z)[q] for d in D for q in nzrange(Z, d)]))
        B = densesub!(at, ZtZ, D, D) \ densesub!(at, Zt, D, nodes)
        for (j, node) in enumerate(nodes), (i, d) in enumerate(D)
            iszero(B[i, j]) || (push!(rr, d); push!(rc, node); push!(rv, B[i, j]))
        end
    end
    return sparse(rr, rc, rv, k, n)
end

# the rate along the projected directions replaced by the derivative at
# the endpoint of the cubic through the state, the two stage increments
# and the endpoint, where a block is on a direction: `v += Z R (vrec - v)`,
# `R` extracting the rate along each direction
function projectrate!(v, pr::ConstraintProjection, pw::ProjectionWork, gc, h, d1, d2, xnew, x)
    w = gc.endrate
    pw.work .= (w[2]/h) .* d1 .+ (w[3]/h) .* d2 .+ (w[4]/h) .* (xnew .- x)
    stepmul!(pw.zwork, pr.Zrate, pw.work)
    stepmul!(pw.zwork2, pr.Zrate, v)
    pw.zwork .-= pw.zwork2
    stepmul!(pw.work, pr.Z, pw.zwork)
    v .+= pw.work
    return v
end

# The small solves of the projection per condition, `M_j^-1 g`, or with
# `transposed` `M_j'^-1 g`, on the columns of condition `j` of the work's
# `g` into its `alpha`, with `M_j = Z' (L + J'(x)) Z` the Jacobian of the
# constraints along the projected directions at the condition's phases
# in `hphi`, refactorized (see `factorprojection!`), or with `kept` on
# the factorization the condition keeps, refactorized where the work's
# `refresh` asks (see `projectendpoint!`); zero for the conditions `open`
# does not hold, when it is given.
function projectionsolve!(pw::ProjectionWork, pr::ConstraintProjection, transposed::Bool, open = nothing; kept::Bool = false)
    N = size(pw.hphi, 2)
    dirs = size(pw.g, 2) ÷ N
    isempty(pw.factors) && return pw.alpha
    for j in 1:N
        a = view(pw.alpha, :, (j - 1)*dirs + 1:j*dirs)
        if !isnothing(open) && !open[j]
            fill!(a, 0)
            continue
        end
        f = length(pw.factors) == 1 ? 1 : j
        if !kept || pw.refresh[j]
            factorprojection!(pw, pr, j, f)
            kept && (pw.factorizations[j] += 1)
        end
        F = pw.factors[f]::KLU.KLUFactorization{Float64,Int}
        a .= view(pw.g, :, (j - 1)*dirs + 1:j*dirs)
        transposed ? ldiv!(transpose(F), a) : ldiv!(F, a)
    end
    return pw.alpha
end

# Condition `j`'s Jacobian at its phases in the work's `hphi`, assembled
# into the work's `jacobian` on the projection's pattern, factorized into
# the work's factorization `f`: refactorized on the pivots of its last,
# or a new one where it has none or `fresh` asks for new pivots.
function factorprojection!(pw::ProjectionWork, pr::ConstraintProjection, j::Int, f::Int; fresh::Bool = false)
    M, d = pw.jacobian, pw.dcol
    values = nonzeros(M)
    copyto!(values, nonzeros(pr.jacobian))
    if !isempty(pr.pj)
        derivativeinto!(d, pr.relationsp, view(pw.hphi, :, j))
        for p in eachindex(d)
            w = pr.lmoljp[p]*d[p]
            for q in pr.jptr[p]:pr.jptr[p + 1] - 1
                values[pr.jidx[q]] += w*pr.jcoef[q]
            end
        end
    end
    F = pw.factors[f]
    pw.factors[f] = (fresh || isnothing(F)) ? kluordered(M, pw.ordering; allowsingular = true) :
        refactorize!(KLUfactorization(), F, M)
    return pw
end

# The factorizations a stepper's projection keeps for its endpoints (see
# `projectendpoint!`) set to its conditions' Jacobians at the state `x`,
# one column per condition: refactorized, or with `fresh` factorized
# anew, where the stepper sets its stage operator to a state, at the
# start and at every checkpoint, so that a replay from a checkpoint
# starts on the factorizations the solve had there.
function refreshprojection!(pw::ProjectionWork, pr::ConstraintProjection, x; fresh::Bool = false)
    isempty(pw.factors) && return pw
    if !isempty(pr.pj)
        stepmul!(pw.phip, pr.RJp, x)
        copyto!(pw.hphi, pw.phip)
    end
    for j in axes(pw.hphi, 2)
        factorprojection!(pw, pr, j, length(pw.factors) == 1 ? 1 : j; fresh)
        pw.factorizations[j] += 1
    end
    return pw
end

# The pattern of the projection's small Jacobians
# `Z' (L + J'(x)) Z = ZLZ + RJZl' Diagonal(lmolj J'(phi)) RJZ`, the
# couplings of the directions through the stiffness and through the
# junctions they share, with the values of `ZLZ`; and for each junction
# `p`, the positions in the pattern's values of the entries its stiffness
# adds to, `jidx[jptr[p]:jptr[p + 1] - 1]`, with the coefficients its
# stiffness multiplies there, so a condition's Jacobian is assembled in
# place, each junction adding to the couplings of the directions on its
# nodes.
function projectionplan(ZLZ::SparseMatrixCSC, RJZl::SparseMatrixCSC, RJZ::SparseMatrixCSC)
    k = size(ZLZ, 1)
    pattern(A) = SparseMatrixCSC(size(A)..., SparseArrays.getcolptr(A), rowvals(A), ones(nnz(A)))
    S = pattern(ZLZ) + sparse(transpose(pattern(RJZl)))*pattern(RJZ)
    M = SparseMatrixCSC(k, k, copy(SparseArrays.getcolptr(S)), copy(rowvals(S)), zeros(nnz(S)))
    position(i, j) = (r = nzrange(M, j); r[searchsortedfirst(view(rowvals(M), r), i)])
    for j in 1:k, q in nzrange(ZLZ, j)
        nonzeros(M)[position(rowvals(ZLZ)[q], j)] += nonzeros(ZLZ)[q]
    end
    RJZlt, RJZt = sparse(transpose(RJZl)), sparse(transpose(RJZ))
    jptr, jidx, jcoef = [1], Int[], Float64[]
    for p in 1:size(RJZ, 1)
        for qa in nzrange(RJZlt, p), qb in nzrange(RJZt, p)
            push!(jidx, position(rowvals(RJZlt)[qa], rowvals(RJZt)[qb]))
            push!(jcoef, nonzeros(RJZlt)[qa]*nonzeros(RJZt)[qb])
        end
        push!(jptr, length(jidx) + 1)
    end
    return M, jptr, jidx, jcoef
end

function projectionwork(pr::ConstraintProjection, backend, n::Int, N::Int, m::Int; kept::Bool = false)
    k, npj = length(pr.directions), length(pr.pj)
    kr, ka = length(pr.readrows), length(pr.auxrows)
    kq = kr + ka
    allocate = (dims...) -> KernelAbstractions.zeros(backend, Float64, dims...)
    # the small Jacobians' factorization, analyzed on their pattern once
    # and first factorized at the junctions' stiffness at zero phase; each
    # condition refactorizes it with its own values. With `kept` each
    # condition keeps a factorization of its own, on the same ordering,
    # first factorized at its own values.
    jacobian = copy(pr.jacobian)
    factors = Union{Nothing,KLU.KLUFactorization{Float64,Int}}[]
    ordering = nothing
    if k > 0
        ordering = fillordering(KLUfactorization(), jacobian)
        if kept
            append!(factors, fill(nothing, N))
        else
            values = nonzeros(jacobian)
            for p in 1:npj, q in pr.jptr[p]:pr.jptr[p + 1] - 1
                values[pr.jidx[q]] += pr.lmoljp[p]*pr.jcoef[q]
            end
            push!(factors, kluordered(jacobian, ordering; allowsingular = true))
        end
    end
    return ProjectionWork(allocate(npj, N), zeros(npj, N), zeros(npj, N), allocate(npj, m), allocate(k, m), allocate(k, m),
        zeros(k, m), zeros(k, m), allocate(k, m), allocate(n, m), allocate(n, m),
        allocate(kq, m), allocate(kq, m), zeros(kq, m), zeros(kq, m), allocate(kr, m), allocate(ka, m), fill(false, N),
        zeros(npj, N), zeros(npj, N), zeros(N), zeros(k, m), jacobian, factors, ordering, fill(Inf, N), fill(false, N),
        zeros(Int, N), zeros(npj), ratework(pr.Minv, m), zeros(k, N), zeros(k, N), zeros(kq, N), zeros(size(pr.Zcinj, 2), N))
end

# The index one unknowns of the endpoint read from their equations: the
# residual `Q (G v + L x + J(x) - b)` from two products on the backend and
# the junction term on the host, the small solve, and the rate along the
# resistive directions and the current on the block rows moved by it.
# `junction` is the current of the projected junctions at the endpoint,
# `Lscale/Lj` times their relation at its phases (or, for a linearized
# reading, times their slope and the flux direction's phase) on the host,
# `gq` the drive along `Q` on the host.
function endpointread!(v, x, pr::ConstraintProjection, pw::ProjectionWork, junction, gq)
    kr = length(pr.readrows)
    isempty(pr.readrows) && isempty(pr.auxrows) && return nothing
    stepmul!(pw.qwork, pr.QG, v)
    stepmul!(pw.qwork2, pr.QL, x)
    pw.qwork .+= pw.qwork2
    copyto!(pw.gq, pw.qwork)
    mul!(view(pw.gq, 1:kr, :), transpose(pr.RJQ), junction, 1.0, 1.0)
    pw.gq .-= gq
    ratesolve!(pw.theta, pr.Minv, pw.gq, pw.ratework)
    pw.theta .*= -1
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

# `Q` of an injection of targets, on the host, for the reading of a
# direction's current at the endpoint time
endpointinjection(pr::ConstraintProjection, injection::SparseMatrixCSC) = hostsparse(pr.Q)*injection

# The transpose of the reading. The cotangents of the rate along the
# resistive directions and of the current on the block rows,
# `[Zr' vbar; Ea' xbar]`, give through the transposed small solve the
# vector `w = Q' zeta` along the rows, which moves the cotangents of the
# rate and the flux before the reading by `G' w` and `(L + J'(x))' w`,
# with the junction stiffness at the endpoint from the cosines `cosp` of
# the projected junctions' phases on the backend, and the current at the
# endpoint time by `-inj' w`. Leaves `w` in the work buffer and returns
# whether there was anything to do.
function endpointreadtranspose!(vbar, xbar, pr::ConstraintProjection, pw::ProjectionWork, sys, cosp, nobj, N)
    kr = length(pr.readrows)
    kq = size(pw.gq, 1)
    kq == 0 && return false
    stepmul!(pw.qwork, pr.Q, vbar)
    stepmul!(pw.qwork2, pr.Q, xbar)
    copyto!(pw.gq, pw.qwork)
    copyto!(pw.theta, pw.qwork2)
    view(pw.gq, kr + 1:kq, :) .= view(pw.theta, kr + 1:kq, :)
    ratesolve!(pw.theta, pr.Minv, pw.gq, pw.ratework; transposed = true)
    pw.theta .*= -1
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
