# The integrator: the trapezoidal rule on the second order system, one
# implicit step at a time, each solved by Newton on the real Jacobian of
# harmonic balance at one mode with the linear term of the step folded in.
# The system is scaled by the problem's `Lscale`, and lives on the backend
# the solve runs on: the products are the package's sparse products, the
# Jacobian is assembled by its plan, and the factorization is its KLU or
# cuDSS.

"""
    AbstractTransientIntegrator

The time stepping rules [`transientsolve`](@ref) chooses between,
[`Trapezoidal`](@ref), [`GaussLegendre`](@ref) and [`BackwardEuler`](@ref).
"""
abstract type AbstractTransientIntegrator end

"""
    Trapezoidal()

The trapezoidal rule on the flux and on its rate (Newmark with the
averaging parameters): second order, and free of numerical damping, so a
lossless LC oscillation keeps its energy. The default of
[`transientsolve`](@ref).
"""
struct Trapezoidal <: AbstractTransientIntegrator end

"""
    BackwardEuler()

Backward Euler on the flux and on its rate: first order and strongly
damping, a deliberately dissipative reference for checking that a result
does not depend on the stepping rule.
"""
struct BackwardEuler <: AbstractTransientIntegrator end

# The coefficients of a rule. A step solves
#
#     [alpha*C + beta*G + L] x + J(x) = b(t) + r(x_n, v_n, ...),
#
# with `C`, `G` and `L` the scaled matrices and `J` the junction current;
# the right hand side `r` and the rate update are the rule's.
stepcoefficients(::Trapezoidal, h) = (4/h^2, 2/h)
stepcoefficients(::BackwardEuler, h) = (1/h^2, 1/h)

# The step's factorization when none is given: KLU on the CPU, and cuDSS
# on a device without the iterative refinement the harmonic balance solver
# asks of it. In the step, Newton's residual check catches what a solve
# leaves; the tangent and the adjoint solve once per step without such a
# check, and the device parity tests are the evidence that this is
# enough on the circuits tried, not a bound on the general case.
transientfactorization(backend::Backend) = backend isa CPU ? KLUfactorization() :
    CUDSSFactorization(ir_n_steps = 0)

# The products of a step. A device product is launched without the full
# device synchronization the package's `mul!` performs after each kernel:
# the work of a step is ordered by the stream it is launched on, and a
# host read (a norm, a copy) waits for it, so the synchronization only
# holds the host back between launches.
stepmul!(y, A::SparseMatrixCSC, x) = mul!(y, A, x)
function stepmul!(y::AbstractVector, A::DeviceValuedSparseMatrix, x::AbstractVector)
    backend = KernelAbstractions.get_backend(A.nzval)
    devicecsrmulkernel!(backend, 64)(y, rowpointer(A), columnindices(A), A.nzval, x;
        ndrange = size(A, 1))
    return y
end
function stepmul!(Y::AbstractMatrix, A::DeviceValuedSparseMatrix, X::AbstractMatrix)
    backend = KernelAbstractions.get_backend(A.nzval)
    devicecsrmulmatrixkernel!(backend, (64, 1))(Y, rowpointer(A), columnindices(A), A.nzval, X;
        ndrange = (size(A, 1), size(X, 2)))
    return Y
end

# a host sparse matrix, or its values on the backend in the row order a
# device product and a device factorization read
devicesparse(A::SparseMatrixCSC, ::CPU) = A
function devicesparse(A::SparseMatrixCSC, backend::Backend)
    At = sparse(transpose(A))
    pattern = DeviceSparsePattern(tobackend(backend, copy(SparseArrays.getcolptr(At))),
        tobackend(backend, copy(rowvals(At))), size(At)...)
    return DeviceValuedSparseMatrix(pattern, tobackend(backend, copy(nonzeros(At))))
end

"""
    TransientSystem

The scaled system of a [`TransientProblem`](@ref) at one step size and
one stepping rule, on one backend: the padded and scaled capacitance,
conductance and augmented inverse inductance matrices, the junction
incidence and coefficients, the drive injection, the port readout, the
real Jacobian plan of harmonic balance at one mode carrying the step's
linear term, the Jacobian it fills and the factorization of it. Built by
[`transientsystem`](@ref); [`transientsolve`](@ref), the tangent and the
adjoint all step on it.
"""
struct TransientSystem{B, M, MJ, V, VC, J, P, F, G, RM, RB}
    problem::TransientProblem
    backend::B
    method::AbstractTransientIntegrator
    h::Float64
    Lscale::Float64
    alpha::Float64
    beta::Float64
    # the scaled matrices, on the backend, and the combinations a step
    # reads: `K = alpha*C + beta*G + L` of the residual and the operator,
    # and `A = alpha*C + beta*G - L` (backward Euler: `C/h^2 + G/h`) and
    # `B = (4/h) C` (backward Euler: `C/h`) of the right hand side, so a
    # residual or a right hand side is one product per matrix rather than
    # one per term, which on a device is one launch rather than five
    C::M
    G::M
    L::M
    K::M
    A::M
    B::M
    # the transposes of the conductance and the stiffness, which the
    # scattering blocks make unsymmetric, and whether they are
    Gt::M
    Lt::M
    symmetric::Bool
    # the junction incidence, its transpose and the coefficients Lscale/Lj
    RJ::MJ
    RJt::MJ
    lmolj::V
    # the drive injection, scaled, and the constant current
    injection::M
    constant::V
    # the rational blocks: the gather of the port rates and the port
    # currents of every block, `(2 nports, n)`, and the scatter of a
    # wave source on their rows, scaled
    blockgather::M
    blockscatter::M
    blockgathert::M
    blockscattert::M
    # the lines: the injection of the currents their arriving waves force,
    # scaled, and the gather of the rates across their ports
    lineinjection::M
    linegather::M
    linescatter::M
    # the port terminals as a difference matrix, the impedances and the
    # termination conductances
    ports::M
    portst::M
    portdrives::M
    portimpedances::V
    portconductances::V
    # the Jacobian plan, its matrix and the factorization
    plan::P
    jacobian::J
    factorization::F
    # the cosine of the junction phases the plan reads
    cosphi::VC
    # The current-phase relation of every junction, on the backend, and a
    # work vector for the Horner loop, which reads the phases while it
    # writes its result. The empty table is every junction sinusoidal,
    # which is every circuit that does not ask for another, and the steps
    # then take the `sin` and `cos` they always did. A table rather than a
    # `nothing` so that the type of a system does not depend on what the
    # circuit holds.
    relations::JunctionRelations{RM, RB}
    relationwork::V
    # the stage data of the Gauss-Legendre rule, nothing for the others; a
    # type parameter, so that the step loop reading the tableau is typed
    gauss::G
end

"""
    TransientReuse()

What a transient solve builds and a later solve, tangent or adjoint of
the same problem at the same step, rule and backend takes over rather
than building again: the scaled [`TransientSystem`](@ref) with its
Jacobian plan, the factorization of the step matrix as it was left, and
the Krylov workspace of the iterative step. Pass one as `reuse` to
[`transientsolve`](@ref), [`transienttangent`](@ref) and
[`transientadjoint`](@ref); a solve whose problem, step, rule or backend
differ from the kept system replaces it. The counterpart of the
harmonic balance solver's reuse between the solves of an [`hbcache`](@ref).
The kept objects are mutable and belong to one solve at a time: do not
share a reuse between concurrent solves.
"""
mutable struct TransientReuse
    system::Any
    factor::Any
    workspace::Any
end
TransientReuse() = TransientReuse(nothing, nothing, nothing)

# the kept system when it is the one asked for, or a new one, kept
function transientsystem(reuse::Union{Nothing,TransientReuse}, p::TransientProblem, h::Real,
        method::AbstractTransientIntegrator, backend::Backend, factorization::AbstractFactorization)
    if !isnothing(reuse) && !isnothing(reuse.system)
        s = reuse.system
        if s.problem === p && s.h == Float64(h) && s.method === method &&
                s.backend == backend && s.factorization == factorization
            return s
        end
        reuse.factor = nothing
        reuse.workspace = nothing
    end
    s = transientsystem(p, h, method, backend, factorization)
    isnothing(reuse) || (reuse.system = s)
    return s
end

"""
    transientsystem(problem, h, method, backend, factorization)

The [`TransientSystem`](@ref) of `problem` at the step `h` under `method`
on `backend`, factorizing with `factorization`.
"""
function transientsystem(p::TransientProblem, h::Real, method::AbstractTransientIntegrator,
        backend::Backend, factorization::AbstractFactorization)
    nm = p.matrices
    Ntot = length(p)
    Nnodal, Naux = p.Nnodal, p.Naux
    # the scale of the equations is the problem's, so that a state means
    # the same at every step
    Lscale = p.Lscale
    C, G, L, lineE = transientlinearmatrices(nm, p.coupledbranches, p.graph.Rbn, p.gaugeindices, p.blocks, p.lines,
        Lscale, Nnodal, Naux)
    symmetric = isempty(p.blocks)
    for l in p.lines
        l.delay >= h || throw(ArgumentError(
            lazy"the transmission line at $(l.path) has a delay of $(l.delay) s, shorter than the step $(h) s; a step reads accepted history only, so reduce dt below the delay."))
    end
    Rbnm = hcat(nm.Rbnm, spzeros(eltype(nm.Rbnm), size(nm.Rbnm, 1), Naux))
    gc = method isa GaussLegendre ? gausscoefficients() : nothing
    if method isa GaussLegendre
        # the real part of the complex stage matrix is the plan's linear
        # term; the rule's endpoint reads no combined matrix
        alpha, beta = real((gc.mu/h)^2), real(gc.mu/h)
        K = spaddkeepzeros(spaddkeepzeros(scaledsparse(alpha, C), scaledsparse(beta, G)), L)
        A, B = K, K
    else
        alpha, beta = stepcoefficients(method, h)
        K = spaddkeepzeros(spaddkeepzeros(scaledsparse(alpha, C), scaledsparse(beta, G)), L)
        A, B = if method isa Trapezoidal
            alpha .* C .+ beta .* G .- L, (4/h) .* C
        else
            alpha .* C .+ beta .* G, (1/h) .* C
        end
    end
    # the junction rows of the incidence matrix and the coefficients
    Ljb = nm.Ljb
    RJ = SparseMatrixCSC{Float64,Int}(Rbnm[Ljb.nzind, :])
    RJt = sparse(transpose(RJ))
    lmolj = Float64[Lscale/Ljb.nzval[i] for i in eachindex(Ljb.nzval)]
    # the drives and the constant current, scaled: a current I enters the
    # node equations as Lscale*I/phi0
    injection = (Lscale/phi0) .* p.injection
    constant = (Lscale/phi0) .* p.constantcurrent
    lineinjection = (Lscale/phi0) .* lineE
    linegather = sparse(transpose(lineE))
    blockgather, blockscatter = transientblockmaps(p, Ntot, Lscale)
    # the port voltages as differences of node rates
    np = length(p.portimpedances)
    prow, pcol, pval = Int[], Int[], Float64[]
    for k in 1:np
        p.portpositive[k] > 0 && (push!(prow, k); push!(pcol, p.portpositive[k]); push!(pval, 1.0))
        p.portnegative[k] > 0 && (push!(prow, k); push!(pcol, p.portnegative[k]); push!(pval, -1.0))
    end
    ports = sparse(prow, pcol, pval, np, Ntot)
    portst = sparse(transpose(ports))
    # the drive current into each port: the drives on the port, summed
    driven = [k for (k, d) in enumerate(p.drives) if d.portindex > 0]
    portdrives = sparse([p.drives[k].portindex for k in driven], driven, ones(length(driven)), np, length(p.drives))
    # the real Jacobian of harmonic balance at one real mode: its pattern
    # is that of the step matrix and the junction pairs, its linear term
    # the step matrix, and the junction term the cosines the plan reads
    devicej = !(backend isa CPU)
    ami = fill(1, 1, 1)
    amc = zeros(Int, 1, 1)
    layout = ModeLayout([true], Ntot)
    Z = spzeros(Float64, Ntot, Ntot)
    Nbranches = p.graph.Nbranches
    Jrs, _ = realjacobianstructure(ami, amc, Ljb, Rbnm, 1, Nbranches, K, Z, Z,
        layout, layout, Float64; transposed = devicej, backend)
    junctions = junctionstructure(Float64, ami, amc, Ljb, Lscale, Rbnm, 1,
        Nbranches, 1, backend)
    wzero = Diagonal(zeros(Ntot))
    plan = planstructurerealjacobian(Jrs, Float64, junctions, K, Z, Z, wzero,
        wzero, layout, layout, backend; transposed = devicej)
    jacobian = devicej ? DeviceValuedSparseMatrix(Jrs, tobackend(backend, zeros(Float64, nnz(Jrs)))) : Jrs
    cosphi = tobackend(backend, ones(ComplexF64, length(Ljb.nzval)))
    d = A -> devicesparse(A, backend)
    v = x -> tobackend(backend, x)
    gauss = if method isa GaussLegendre
        imvals = gaussimaginary(gc, h, C, G, Jrs, devicej)
        cjacobian = devicej ? DeviceValuedSparseMatrix(Jrs, tobackend(backend, zeros(ComplexF64, nnz(Jrs)))) :
            SparseMatrixCSC(size(Jrs)..., copy(SparseArrays.getcolptr(Jrs)), copy(rowvals(Jrs)), zeros(ComplexF64, nnz(Jrs)))
        rationalvals = rationalvalues(p, gc.mu/h, Lscale, Jrs, devicej)
        stages = rationalstages(p, gc, h)
        coupling = isempty(stages) ? nothing : rationalcoupling(p, stages, gc, h, Lscale, blockgather, blockscatter, backend)
        gaussstage(gc, v(imvals), cjacobian, gaussprojection(p, G, L, RJ, lmolj, injection, constant, lineinjection, blockscatter, backend),
            v(rationalvals), stages, coupling, backend)
    else
        nothing
    end
    return TransientSystem(p, backend, method, Float64(h), Lscale, alpha, beta,
        d(C), d(G), d(L), d(K), d(A), d(B), d(sparse(transpose(G))), d(sparse(transpose(L))), symmetric,
        d(RJ), d(RJt), v(lmolj), d(injection), v(constant), d(blockgather), d(blockscatter),
        d(sparse(transpose(blockgather))), d(sparse(transpose(blockscatter))), d(lineinjection), d(linegather), d(lineE),
        d(ports), d(portst), d(portdrives), v(p.portimpedances), v(p.portconductances), plan, jacobian,
        factorization, cosphi, torelations(p.relations, v(zeros(0))),
        v(zeros(length(Ljb.nzval))), gauss)
end

# a sparse matrix scaled with every stored entry kept, zero or not: the
# pattern of the step operator must hold every entry a block's rows can
# take, and a broadcast drops the ones that vanish at a particular value,
# as an ideal through's conductance entries do
scaledsparse(a, A::SparseMatrixCSC) = SparseMatrixCSC(A.m, A.n, copy(SparseArrays.getcolptr(A)), copy(rowvals(A)), a .* nonzeros(A))

# The maps of the rational blocks: the gather of the rate across each
# port and of each port current, as `(2 nports, n)` with the rates first,
# and the scatter of a source `2 eta` on the block's port current rows,
# scaled as the rows are, `(n, nports)`
function transientblockmaps(p::TransientProblem, Ntot::Int, Lscale::Float64)
    gi, gj, gv = Int[], Int[], Float64[]
    si, sj, sv = Int[], Int[], Float64[]
    offset = 0
    for b in p.blocks
        n = length(b.signal)
        for q in 1:n
            b.signal[q] > 0 && (push!(gi, offset + q); push!(gj, b.signal[q]); push!(gv, 1.0))
            b.ref[q] > 0 && (push!(gi, offset + q); push!(gj, b.ref[q]); push!(gv, -1.0))
            push!(gi, offset + n + q); push!(gj, b.auxbase + q); push!(gv, 1.0)
            push!(si, b.auxbase + q); push!(sj, offset ÷ 2 + q); push!(sv, 2Lscale/phi0)
        end
        offset += 2n
    end
    nports = sum(b -> length(b.signal), p.blocks; init = 0)
    return sparse(gi, gj, gv, 2nports, Ntot), sparse(si, sj, sv, Ntot, nports)
end

# The scaled linear matrices of a problem: the augmentation of the
# coupled inductors' constitutive rows and the gauge rows folded into
# the inverse inductance, as hbnlsolve does, the padded and scaled
# capacitance, conductance and stiffness, the scattering blocks'
# constitutive rows on the rates and the port currents with the
# currents' Kirchhoff couplings, and the lines' conductances at their
# impedance; also the lines' incidence. The same matrices classify the
# directions of a problem at its construction and step it.
function transientlinearmatrices(nm::CircuitMatrices, coupledbranches, Rbn, gaugeindices, blocks, lines, Lscale, Nnodal, Naux)
    Ntot = Nnodal + Naux
    AmnaL = calcAmnaind(coupledbranches, nm.Lb, nm.Mb, Rbn, 1, Nnodal, Ntot, Lscale)
    Amna = real.(spaddkeepzeros(calcAmna(gaugeindices, Ntot), AmnaL))
    C = mnapad(Lscale .* real.(nm.Cnm), Naux)
    G = mnapad(Lscale .* real.(nm.Gnm), Naux)
    L = spaddkeepzeros(mnapad(Lscale .* real.(nm.invLnm), Naux), Amna)
    if !isempty(blocks)
        blockG, blockL = transientblockstamps(blocks, Ntot, Lscale)
        G = spaddkeepzeros(G, blockG)
        L = spaddkeepzeros(L, blockL)
    end
    if !isempty(lines)
        lineG, lineE = transientlinestamps(lines, Ntot, Lscale)
        G = spaddkeepzeros(G, lineG)
    else
        lineE = spzeros(Ntot, 0)
    end
    return C, G, L, lineE
end

# The stamps of the lines: the conductance `Lscale/Z` across each port,
# and the incidence of a port's forced current into its signal node and
# out of its reference node, whose transpose reads the rate across the
# port. Returns the conductance contribution and the incidence.
function transientlinestamps(lines::Vector{TransientLine}, Ntot::Int, Lscale::Float64)
    gi, gj, gv = Int[], Int[], Float64[]
    ei, ej, ev = Int[], Int[], Float64[]
    for (l, line) in enumerate(lines), q in 1:2
        g = Lscale/line.Z
        s, r = line.signal[q], line.ref[q]
        col = 2(l - 1) + q
        s > 0 && (push!(gi, s); push!(gj, s); push!(gv, g); push!(ei, s); push!(ej, col); push!(ev, 1.0))
        r > 0 && (push!(gi, r); push!(gj, r); push!(gv, g); push!(ei, r); push!(ej, col); push!(ev, -1.0))
        (s > 0 && r > 0) && (push!(gi, s); push!(gj, r); push!(gv, -g); push!(gi, r); push!(gj, s); push!(gv, -g))
    end
    return sparse(gi, gj, gv, Ntot, Ntot), sparse(ei, ej, ev, Ntot, 2length(lines))
end

# The stamps of the scattering blocks in the scaled equations. A block's
# constitutive row for port `p` is `Lscale (I - S) R^(-1/2) E' xdot -
# (I + S) R^(1/2) u = 0`, the hybrid equation of the linearized solver
# with the voltage `phi0 xdot` and the current `phi0 u / Lscale` in the
# units of the state, so its conductance entries go to `G` on the nodes
# of the ports and its current entries to `L` on the port unknowns; the
# current of a port enters its signal node and leaves its reference node
# in the node equations, as the branch currents of the coupled inductors
# do. Returns the conductance and the stiffness contributions.
function transientblockstamps(blocks::Vector{TransientBlock}, Ntot::Int, Lscale::Float64)
    gi, gj, gv = Int[], Int[], Float64[]
    li, lj, lv = Int[], Int[], Float64[]
    for b in blocks
        n = length(b.signal)
        Bb = Lscale .* (I - b.S) .* transpose(1 ./ sqrt.(b.R))
        Cb = (I + b.S) .* transpose(sqrt.(b.R))
        for q in 1:n
            aux = b.auxbase + q
            b.signal[q] > 0 && (push!(li, b.signal[q]); push!(lj, aux); push!(lv, 1.0))
            b.ref[q] > 0 && (push!(li, b.ref[q]); push!(lj, aux); push!(lv, -1.0))
            for r in 1:n
                push!(li, aux); push!(lj, b.auxbase + r); push!(lv, -Cb[q, r])
                b.signal[r] > 0 && (push!(gi, aux); push!(gj, b.signal[r]); push!(gv, Bb[q, r]))
                b.ref[r] > 0 && (push!(gi, aux); push!(gj, b.ref[r]); push!(gv, -Bb[q, r]))
            end
        end
    end
    return sparse(gi, gj, gv, Ntot, Ntot), sparse(li, lj, lv, Ntot, Ntot)
end

# the junction phases of a state and the junction current they drive,
# `RJ' * (Lscale/Lj .* sin.(RJ*x))`; the phases are kept for the Jacobian
function junctionphases!(phi, sys::TransientSystem, x)
    stepmul!(phi, sys.RJ, x)
    return phi
end
function junctioncurrent!(y, sys::TransientSystem, phi, work)
    r = sys.relations
    if allsinusoidal(r)
        work .= sys.lmolj .* sin.(phi)
    else
        relationinto!(sys.relationwork, r, phi)
        work .= sys.lmolj .* sys.relationwork
    end
    stepmul!(y, sys.RJt, work)
    return y
end

# the drift `G*v + L*x + J(x)` at a state, the part of the equations that
# is not the inertia, with the junction current given
function transientdrift!(y, sys::TransientSystem, x, v, junction, work)
    stepmul!(y, sys.G, v)
    stepmul!(work, sys.L, x)
    y .+= work .+ junction
    return y
end

# the scaled node current of the drives at a time, and the drive values
function drivecurrent!(b, sys::TransientSystem, values, hostvalues, t)
    p = sys.problem
    @inbounds for (k, d) in enumerate(p.drives)
        hostvalues[k] = d.current(t)
    end
    all(isfinite, hostvalues) || throw(ArgumentError(lazy"a source returned a nonfinite current at t = $(t) s."))
    values === hostvalues || copyto!(values, hostvalues)
    stepmul!(b, sys.injection, values)
    b .+= sys.constant
    return b
end

# the Jacobian of the step at the phases, assembled by the plan from the
# cosines and factorized, or refactorized on the analysis of the first
function stepjacobian!(sys::TransientSystem, phi, factor)
    derivativeinto!(sys.cosphi, sys.relations, phi)
    assemblerealjacobian!(nonzeros(sys.jacobian), sys.plan, sys.cosphi)
    if isnothing(factor)
        return factorize(sys.factorization, sys.jacobian)
    else
        return refactorize!(sys.factorization, factor, sys.jacobian)
    end
end

"""
    TransientSolution

Result of [`transientsolve`](@ref). `times` is in seconds; `voltage`
holds the port voltages in Volts and `incident` and `outgoing` the real
instantaneous power waves in sqrt(W), one row per port in compiled port
order and one column per saved time. `initialflux`, `initialrate`,
`finalflux` and `finalrate` always hold the first and the last scaled
state. With `record = :phases` or `:states`, `phases` holds the junction
phases the tangent and the adjoint read, at every saved time under the
trapezoidal rule and at the two stages of each step under
[`GaussLegendre`](@ref), as `(junction, stage, time)`, and `endphases`
the phases at the projected endpoints of the junctions on the algebraic
directions the rule projects, or `nothing` without any, and `linewaves`
the wave leaving each port of each transmission line at every step, as
`(port, time)` in sqrt(W), or `nothing` without lines, or with
checkpoints when the history before each of them, kept as their
`waves`, is smaller than the record, which is kept instead once the
longest delay exceeds the checkpoint interval; on the solve's backend
as every record is; with
`record = :states`, `flux` and `rate` hold the scaled node fluxes and
their rates at every saved time as well; with `record = :checkpoints`,
`checkpoints` holds the state and the stage predictor every
`checkpointevery` steps, from which the responses replay the steps
between, so that a record of any length costs the checkpoints and one
window of phases. `stats` counts the
steps, the Newton corrections, the numeric factorizations, the retries
of a rejected correction after a fresh factorization at its base point,
and the Krylov iterations. The arrays live on the backend of the solve.
"""
struct TransientSolution{P, M, V}
    problem::P
    method::AbstractTransientIntegrator
    dt::Float64
    times::Vector{Float64}
    voltage::M
    incident::M
    outgoing::M
    phases::Any
    endphases::Any
    linewaves::Any
    # the flux and rate records untyped, as the phases are, so that the
    # record level does not make a new solution type and the responses
    # compile once
    flux::Any
    rate::Any
    checkpoints::Any
    initialflux::V
    initialrate::V
    finalflux::V
    finalrate::V
    blockstates::Any
    # the initial waves of the lines and states of the blocks, the rest
    # of the state the solve started from, which the noise checks
    initialwaves::Any
    initialstates::Any
    stats::NamedTuple
end

# the levels of the record: the port waves, the junction phases the
# responses read, or the whole state history
function recordlevel(record::Symbol, saveevery)
    record in (:ports, :phases, :states, :checkpoints) || throw(ArgumentError("record must be :ports, :phases, :states or :checkpoints."))
    (record == :ports || saveevery == 1) || throw(ArgumentError("a record of the phases, the states or the checkpoints needs saveevery = 1, since the derivatives read every step."))
    return record in (:phases, :states), record == :states, record == :checkpoints
end

# the port voltages and waves of a state, from the node rates and the
# drive currents into the ports
function portwaves!(voltage, incident, outgoing, sys::TransientSystem, v, drivevalues, work)
    stepmul!(work, sys.ports, v)
    voltage .= phi0 .* work
    stepmul!(work, sys.portdrives, drivevalues)
    current = work .- sys.portconductances .* voltage
    z = sys.portimpedances
    incident .= (voltage .+ z .* current) ./ (2 .* sqrt.(z))
    outgoing .= (voltage .- z .* current) ./ (2 .* sqrt.(z))
    return nothing
end

# the checks of the arguments every solve makes, and the uniform grid
function transientgrid(tspan, dt, maxsteps, saveevery, maxiters, rtol, atol, linearsolver, reuse)
    length(tspan) == 2 || throw(ArgumentError("tspan must have two endpoints."))
    t0, tf = Float64(tspan[1]), Float64(tspan[2])
    (isfinite(t0) && isfinite(tf) && tf > t0) || throw(ArgumentError("tspan must be finite and increasing."))
    (isfinite(dt) && dt > 0 && t0 + dt > t0) || throw(ArgumentError("dt must be finite, positive, and resolvable at tspan[1]."))
    (saveevery > 0 && maxiters > 0 && maxsteps > 0) || throw(ArgumentError("saveevery, maxiters and maxsteps must be positive."))
    (isfinite(rtol) && isfinite(atol) && rtol >= 0 && atol > 0) || throw(ArgumentError("rtol must be nonnegative and atol positive."))
    count = (tf - t0)/Float64(dt)
    (isfinite(count) && count <= maxsteps) || throw(ArgumentError(lazy"the transient needs more than maxsteps = $(maxsteps) steps."))
    nearest = round(count)
    abs(count - nearest) <= 8eps(count) && (count = nearest)
    nsteps = max(1, ceil(Int, count))
    h = (tf - t0)/nsteps
    isnothing(linearsolver) || linearsolver isa AbstractHBLinearSolver || throw(ArgumentError(
        "linearsolver must be nothing (a direct solve) or one of the package's Krylov solvers, such as GMRES()."))
    isnothing(reuse) || reuse isa TransientReuse || throw(ArgumentError("reuse must be a TransientReuse or nothing."))
    return t0, tf, nsteps, h
end

"""
    transientsolve(problems::AbstractVector{TransientProblem}, tspan; dt,
        method = GaussLegendre(), backend = CPU(), factorization = nothing,
        reuse = nothing, initialstate = transientstate of each, saveevery = 1,
        record = :ports, checkpointevery = 0, rtol = 1e-9, atol = 1e-10,
        maxiters = 15, maxsteps = 10^7)

The same circuit under every drive condition of `problems`, built by
[`transientproblem`](@ref)`(problem; sources)` so that they share one
compiled circuit and drive the same targets, stepped as one system under
[`GaussLegendre`](@ref): the states are matrices with the conditions as
columns, every product takes all conditions in one call, the junction
stiffness and the factorization are one per condition, KLU on the CPU
and the uniform cuDSS batch on a device, and the Newton engine accepts
and refreshes per condition. On a device this is where the throughput
lies, since the launches of a step serve every condition. `initialstate`
is one pair for all conditions or a vector of pairs. Returns a
[`TransientBatchSolution`](@ref), whose `solution[j]` is the ordinary
solution of condition `j`.

On the host the conditions are also split across the threads of the
session and their chunks stepped at once, since the conditions of a batch
are independent of one another. The whole step parallelizes that way, not
only the assembly, the factorization and the solve, and the chunks fill
the arrays of the batch in place, so the split costs no copy and changes
no result: any layout of chunks gives the same bits. Start Julia with
`-t` to use it. A device keeps the one chunk its uniform batch already
is.
"""
function transientsolve(problems::AbstractVector{TransientProblem}, tspan; dt::Real,
        method::AbstractTransientIntegrator = GaussLegendre(), backend::Backend = CPU(),
        factorization = nothing, reuse = nothing, initialstate = nothing,
        saveevery::Integer = 1, record::Symbol = :ports, checkpointevery::Integer = 0, rtol::Real = 1e-9,
        atol::Real = 1e-10, maxiters::Integer = 15, maxsteps::Integer = 10^7)
    p = batchcompatible(problems)
    method isa GaussLegendre || throw(ArgumentError("a batch of conditions steps under GaussLegendre()."))
    t0, tf, nsteps, h = transientgrid(tspan, dt, maxsteps, saveevery, maxiters, rtol, atol, nothing, reuse)
    states = if isnothing(initialstate)
        [transientstate(q) for q in problems]
    elseif initialstate isa AbstractVector && !isempty(initialstate) && first(initialstate) isa Union{Tuple,NamedTuple}
        length(initialstate) == length(problems) || throw(DimensionMismatch("give one initial state per condition, or one for all."))
        collect(initialstate)
    else
        fill(initialstate, length(problems))
    end
    fact = isnothing(factorization) ? transientfactorization(backend) : factorization
    sys = transientsystem(reuse, p, h, method, backend, fact)
    return gaussbatchintegrate(sys, collect(problems), t0, tf, nsteps, states, saveevery, record, checkpointevery,
        Float64(rtol), Float64(atol), maxiters, reuse)
end

"""
    transientsolve(problem, tspan; dt, method = Trapezoidal(), backend = CPU(),
        factorization = nothing, linearsolver = nothing,
        initialstate = transientstate(problem), saveevery = 1,
        record = :ports, checkpointevery = 0, rtol = 1e-9, atol = 1e-10,
        maxiters = 15, maxsteps = 10^7)

Integrate the circuit in physical time on a uniform grid of step at most
`dt`, shortened slightly to land on `tspan[2]`. `method` is the stepping
rule, [`Trapezoidal`](@ref) by default, [`GaussLegendre`](@ref) for
fourth order at the same step, or [`BackwardEuler`](@ref); `backend` is
where the solve runs, and `factorization` the sparse factorization, KLU
on the CPU and cuDSS on a CUDA device by default. The
Jacobian's pattern is fixed and its symbolic analysis done once; a linear
circuit factorizes once. `rtol` and `atol` control the Newton residual,
not the temporal error, and a Newton step that fails throws.

With `linearsolver = GMRES()` (or another of the package's Krylov
solvers) each Newton correction is solved iteratively and matrix free,
with the last factorization of the step matrix as the preconditioner: the
factorization is then refreshed only when the iteration count says it
has drifted, rather than at every step the junction phases move. The
Krylov workspace and the preconditioner persist for the whole solve.

By default only the port outputs and the first and last states are
kept, which is what the port responses need; `saveevery` decimates the
saved samples without changing the grid and keeps both endpoints.
`record = :phases` with `saveevery = 1` records the junction phases at
every step, the least [`transienttangent`](@ref),
[`transientadjoint`](@ref) and the noise need, `record = :states` the
whole flux and rate history as well, and `record = :checkpoints`, under
[`GaussLegendre`](@ref), only the state every `checkpointevery` steps
(the square root of the step count by default), from which the
responses replay the steps between at the cost of one more solve, so
the memory of a record of any length is bounded.
`initialstate` is the pair of [`transientstate`](@ref); the solver checks
that it satisfies the algebraic equations of the circuit at the start.
"""
function transientsolve(p::TransientProblem, tspan; dt::Real,
        method::AbstractTransientIntegrator = Trapezoidal(), backend::Backend = CPU(),
        factorization = nothing, linearsolver = nothing, reuse = nothing,
        initialstate = transientstate(p),
        saveevery::Integer = 1, record::Symbol = :ports, checkpointevery::Integer = 0, rtol::Real = 1e-9,
        atol::Real = 1e-10, maxiters::Integer = 15, maxsteps::Integer = 10^7)
    t0, tf, nsteps, h = transientgrid(tspan, dt, maxsteps, saveevery, maxiters, rtol, atol, linearsolver, reuse)
    (isempty(p.blocks) && isempty(p.lines)) || method isa GaussLegendre || throw(ArgumentError(
        "a circuit with scattering blocks or transmission lines steps under GaussLegendre()."))
    fact = isnothing(factorization) ? transientfactorization(backend) : factorization
    sys = transientsystem(reuse, p, h, method, backend, fact)
    if method isa GaussLegendre
        isnothing(linearsolver) || throw(ArgumentError("the Gauss-Legendre rule solves its stages on its complex factorization; it takes no linearsolver."))
        return unbatch(gaussbatchintegrate(sys, [p], t0, tf, nsteps, [initialstate], saveevery, record, checkpointevery,
            Float64(rtol), Float64(atol), maxiters, reuse))
    end
    record == :checkpoints && throw(ArgumentError("checkpoints are a record of the Gauss-Legendre rule."))
    return transientintegrate(sys, t0, tf, nsteps, initialstate, saveevery, record,
        Float64(rtol), Float64(atol), maxiters, linearsolver, reuse)
end

# The scaled node current of the drives and the constant sources of a
# problem at a time, on the host, for the checks made once at the start;
# the problem is a condition's, which the system holds only in its
# injection and its scale
function hostdrivecurrent(sys::TransientSystem, t, p::TransientProblem = sys.problem, linevalues = zeros(2length(p.lines)),
        blockwaves = nothing)
    values = [d.current(t) for d in p.drives]
    all(isfinite, values) || throw(ArgumentError(lazy"a source returned a nonfinite current at t = $(t) s."))
    b = hostsparse(sys.injection)*values .+ (sys.Lscale/phi0) .* p.constantcurrent .+ hostsparse(sys.lineinjection)*linevalues
    isnothing(blockwaves) || (b .+= hostsparse(sys.blockscatter)*blockwaves)
    return b
end

# The consistency of a state with the algebraic equations of a condition.
# Along each inertialess direction `z` of the problem the equations are
# `z' (G v + L x + J(x)) = z' b(t)`, which nothing integrates, and along
# each algebraic direction, where the conductance vanishes as well, the
# rate is read by the same equation differentiated,
# `z' ((L + J'(x)) v) = z' b'(t)`, so a rate violating it would ring
# rather than decay. Returns the largest violation of both and the scale
# to compare it with, on the host.
function transientconsistency(sys::TransientSystem, x, v, t, p::TransientProblem = sys.problem, linevalues = zeros(2length(p.lines)),
        blockwaves = nothing)
    isempty(p.inertialess) && return 0.0, 1.0
    G, L, RJ = hostsparse(sys.G), hostsparse(sys.L), hostsparse(sys.RJ)
    lmolj = Array(sys.lmolj)
    xh, vh = Array(x), Array(v)
    b = hostdrivecurrent(sys, t, p, linevalues, blockwaves)
    phi = RJ*xh
    hr = hostrelations(sys.relations)
    gv, lx, j = G*vh, L*xh, transpose(RJ)*(lmolj .* relationat(hr, phi))
    r = gv .+ lx .+ j .- b
    # the violations are relative to the terms balanced
    scale = max(norm(b, Inf), norm(gv, Inf), norm(lx, Inf), norm(j, Inf), 1.0)
    worst = 0.0
    # the rate of the drives, by a central difference far below the step
    delta = 1e-3*sys.h
    bdot = (hostdrivecurrent(sys, t + delta, p, linevalues, blockwaves) .- hostdrivecurrent(sys, t - delta, p, linevalues, blockwaves)) ./ (2delta)
    lv, jv = L*vh, transpose(RJ)*(lmolj .* derivativeat(hr, phi) .* (RJ*vh))
    ratescale = max(norm(bdot, Inf), norm(lv, Inf), norm(jv, Inf), 1.0)
    lv .+= jv .- bdot
    for z in p.inertialess
        worst = max(worst, abs(sum(view(r, z)))/scale)
    end
    # the constraints differentiated: the currents cancel in them
    size(p.constraints, 1) == 0 || (worst = max(worst, norm(p.constraints*lv, Inf)/ratescale))
    return worst, 1.0
end

function transientintegrate(sys::TransientSystem, t0, tf, nsteps, initialstate,
        saveevery, record, rtol, atol, maxiters, linearsolver = nothing, reuse = nothing)
    savephases, savestates, _ = recordlevel(record, saveevery)
    p = sys.problem
    backend = sys.backend
    n, np, nd = length(p), length(p.portimpedances), length(p.drives)
    h, alpha, beta = sys.h, sys.alpha, sys.beta
    trapezoidal = sys.method isa Trapezoidal
    x0, v0 = initialstate
    (length(x0) == n && length(v0) == n) || throw(DimensionMismatch(lazy"the state has $(n) entries; use transientstate."))
    (all(isfinite, x0) && all(isfinite, v0)) || throw(ArgumentError("the initial state must be finite."))
    allocate = () -> KernelAbstractions.zeros(backend, Float64, n)
    x, v = tobackend(backend, Float64.(collect(x0))), tobackend(backend, Float64.(collect(v0)))
    xnew, residual, rhs, work = allocate(), allocate(), allocate(), allocate()
    junction, junctionnew, correction, trial = allocate(), allocate(), allocate(), allocate()
    # a trial of the line search keeps its own residual, phases and
    # junction current, so that the base point's stay what a refreshed
    # Jacobian and correction are built from
    trialresidual, trialjunction = allocate(), allocate()
    b, bprev = allocate(), allocate()
    nj = length(sys.lmolj)
    phi, phinew, trialphi, jwork = [KernelAbstractions.zeros(backend, Float64, nj) for _ in 1:4]
    hostvalues = zeros(nd)
    values = tobackend(backend, zeros(nd))
    prevvalues = zeros(nd)
    portwork = KernelAbstractions.zeros(backend, Float64, np)
    # the drive and the junction current at the start, and the check of
    # the algebraic equations along the inertialess directions
    drivecurrent!(b, sys, values, hostvalues, t0)
    copyto!(prevvalues, hostvalues)
    junctionphases!(phi, sys, x)
    junctioncurrent!(junction, sys, phi, jwork)
    violation, _ = transientconsistency(sys, x, v, t0)
    violation <= atol + rtol || throw(ArgumentError(
        "the initial state violates the algebraic equations of the circuit along a direction without capacitance to ground (a node no capacitor touches, a capacitive island, a coupled inductor or gauge row); supply a consistent transientstate, or start the drive from an equilibrium."))
    # The Jacobian and its factorization at the start. A factorization is
    # kept across steps and corrections while Newton converges in one
    # correction with it, which the predictor from the last increment makes
    # the rule for a junction whose phase moves little per step; when a
    # second correction is needed, or the first fails, the Jacobian is
    # reassembled at the current iterate and factorized, and the next step
    # then assembles at its predictor before its first correction rather
    # than trying the stale factorization again. A linear circuit therefore
    # factorizes once, a weakly driven junction nearly so, and a strongly
    # driven one once per step with one correction.
    # a kept factorization is refreshed rather than analyzed again
    factor = stepjacobian!(sys, phi, isnothing(reuse) ? nothing : reuse.factor)
    factorizations, corrections, retries = 1, 0, 0
    nonlinear = nj > 0
    stalefailed = false
    # The iterative step: the step matrix applied matrix free at the
    # current iterate, the last factorization applied as the preconditioner,
    # and the workspace kept for the whole solve. The factorization is
    # refreshed when a solve needs more than a few iterations, or fails,
    # which is the evidence that the phases have moved away from it.
    iterative = !isnothing(linearsolver)
    products = allocate()
    operator = FunctionOperator((y, v) -> begin
        stepmul!(y, sys.K, v)
        junctionproduct!(products, sys, phinew, v, jwork); y .+= products
        return y
    end, n)
    workspace = if !iterative
        nothing
    elseif !isnothing(reuse) && reuse.workspace isa GMRESWorkspace && size(reuse.workspace.V, 1) == n
        reuse.workspace
    else
        GMRESWorkspace(residual, 20)
    end
    preconditioner = (z, r) -> myldiv!(z, factor, r)
    krylovlimit = 4
    krylov = 0
    # what the Newton engine calls: the residual at the base point and at
    # a trial, each keeping its own phases and junction current, the
    # refresh of the factorization at the base point's phases, and the
    # solve with the kept factorization, direct or preconditioned Krylov
    baseresidual! = (norms, r, y) -> (norms[1] = stepresidual!(r, sys, y, phinew, junctionnew, jwork, work, rhs); nothing)
    trialresidual! = (norms, r, y) -> (norms[1] = stepresidual!(r, sys, y, trialphi, trialjunction, jwork, work, rhs); nothing)
    # the refresh refactorizes in place: KLU and cuDSS both return the
    # factorization they were given, so the closure need not rebind it
    refresh! = () -> (stepjacobian!(sys, phinew, factor); nothing)
    accept! = mask -> (copyto!(phinew, trialphi); copyto!(junctionnew, trialjunction); nothing)
    newtonwork = NewtonWork(backend, 1)
    tolerance = [0.0]
    solve! = if iterative
        (c, r) -> begin
            out = hblinearsolve!(linearsolver, c, operator, r, workspace, preconditioner;
                rtol = 1e-4, atol = 0.0, maxrestarts = 5)
            (nonlinear && (!out.converged || out.iterations > krylovlimit), out.iterations)
        end
    else
        (c, r) -> (myldiv!(c, factor, r); (false, 0))
    end
    # the saved outputs
    nsaved = cld(nsteps, saveevery) + 1
    times = Vector{Float64}(undef, nsaved)
    voltage, incident, outgoing = [KernelAbstractions.allocate(backend, Float64, np, nsaved) for _ in 1:3]
    phases = savephases ? KernelAbstractions.allocate(backend, Float64, nj, nsaved) : nothing
    flux = savestates ? KernelAbstractions.allocate(backend, Float64, n, nsaved) : nothing
    rate = savestates ? KernelAbstractions.allocate(backend, Float64, n, nsaved) : nothing
    times[1] = t0
    portwaves!(view(voltage, :, 1), view(incident, :, 1), view(outgoing, :, 1), sys, v, values, portwork)
    savephases && copyto!(view(phases, :, 1), phi)
    if savestates
        copyto!(view(flux, :, 1), x)
        copyto!(view(rate, :, 1), v)
    end
    initialflux, initialrate = copy(x), copy(v)
    saved = 1
    lastincrement = allocate()
    for step in 1:nsteps
        t = step == nsteps ? tf : t0 + step*h
        copyto!(bprev, b)
        drivecurrent!(b, sys, values, hostvalues, t)
        # the right hand side of the step from the previous state. For the
        # trapezoidal rule on x' = v and C v' = b - G v - L x - J(x), with
        # the rate update v_{n+1} = (2/h)(x_{n+1} - x_n) - v_n substituted,
        #   (b_{n+1} + b_n) + [alpha*C + beta*G - L] x_n + (4/h) C v_n - J(x_n);
        # for backward Euler
        #   b_{n+1} + C (x_n/h^2 + v_n/h) + G x_n/h.
        stepmul!(rhs, sys.A, x)
        stepmul!(work, sys.B, v)
        if trapezoidal
            rhs .+= work .+ b .+ bprev .- junction
        else
            rhs .+= work .+ b
        end
        scale = max(norm(rhs, Inf), 1.0)
        tolerance[1] = atol + rtol*scale
        # Newton on x_{n+1}, from the previous increment
        xnew .= x .+ lastincrement
        converged, fresh, ncorr, nfact, nretry, nkrylov, stale = newtonsolve!(xnew, correction, trial,
            residual, trialresidual, baseresidual!, trialresidual!, refresh!, solve!,
            tolerance, maxiters, stalefailed, nonlinear, iterative, newtonwork; accept!)
        corrections += ncorr
        factorizations += nfact
        retries += nretry
        krylov += nkrylov
        converged || error(lazy"the Newton solve of step $(step) at t = $(t) s did not converge; reduce dt or check the initial state and the circuit.")
        stalefailed = iterative ? stale : fresh
        # the rate of the new state and the state itself
        if trapezoidal
            v .= (2/h) .* (xnew .- x) .- v
        else
            v .= (xnew .- x) ./ h
        end
        lastincrement .= xnew .- x
        copyto!(x, xnew)
        copyto!(phi, phinew)
        copyto!(junction, junctionnew)
        if step % saveevery == 0 || step == nsteps
            saved += 1
            times[saved] = t
            portwaves!(view(voltage, :, saved), view(incident, :, saved), view(outgoing, :, saved), sys, v, values, portwork)
            savephases && copyto!(view(phases, :, saved), phi)
            if savestates
                copyto!(view(flux, :, saved), x)
                copyto!(view(rate, :, saved), v)
            end
        end
    end
    KernelAbstractions.synchronize(backend)
    if !isnothing(reuse)
        reuse.factor = factor
        reuse.workspace = workspace
    end
    return TransientSolution(p, sys.method, h, times, voltage, incident, outgoing,
        phases, nothing, nothing, flux, rate, nothing, initialflux, initialrate, copy(x), copy(v), nothing, nothing, nothing,
        (; steps = nsteps, newtoncorrections = corrections, factorizations,
            retries, kryloviterations = krylov, rtol, atol, maxiters))
end

# the residual of the step at a trial state, `K*x + J(x) - rhs`, with the
# phases and the junction current of the trial kept for the Jacobian and
# the next step
function stepresidual!(residual, sys::TransientSystem, x, phi, junction, jwork, work, rhs)
    junctionphases!(phi, sys, x)
    junctioncurrent!(junction, sys, phi, jwork)
    stepmul!(residual, sys.K, x)
    residual .+= junction .- rhs
    return norm(residual, Inf)
end

# The Newton solve of one implicit step, shared by the stepping rules,
# on a batch of conditions held as the columns of the state. `x` is the
# base point and is left at the converged state; `baseresidual!(norms, r, y)`
# and `trialresidual!(norms, r, y)` evaluate the residual at `y` into `r`
# and its norm per column into `norms`, each keeping its own cache of what
# the Jacobian reads, so that a rejected trial never overwrites the base
# point's; `refresh!()` assembles and factorizes the Jacobians at the base
# point; `solve!(c, r)` solves for the correction and returns whether the
# solve found the factorization stale, and its iteration count;
# `accept!(mask)` makes the trial's cache the base point's on the columns
# of `mask`, whose residual is then adopted rather than evaluated again.
# The policy: a kept factorization is tried for the first correction
# unless the last step found it stale; a second correction, or a rejected
# first one, refreshes it at the base point and rebuilds the correction
# from the base point's residual; the line search halves a rejected
# column's correction. With `simplified` the factorization is an
# approximation of the Jacobian even when fresh, as the Gauss-Legendre
# stage matrix is, so a second correction is the rule rather than a sign
# of staleness; the refresh is then decided by the contraction of the
# residual between corrections, a poor one on any column meaning the
# frozen operator has drifted. Returns whether every column converged,
# whether the factorization was refreshed, and the counts. `work` holds
# the host vectors of the per column bookkeeping and the mask on the
# backend.
struct NewtonWork{V, M}
    norms::Vector{Float64}
    trialnorms::Vector{Float64}
    previous::Vector{Float64}
    steps::Vector{Float64}
    stepsdev::V
    mask::M
    accepted::Vector{Bool}
    newly::Vector{Bool}
end
function NewtonWork(backend, ncolumns)
    return NewtonWork(zeros(ncolumns), zeros(ncolumns), fill(Inf, ncolumns), ones(ncolumns),
        tobackend(backend, ones(ncolumns)), tobackend(backend, fill(false, ncolumns)), fill(false, ncolumns), fill(false, ncolumns))
end
# the columns of `dst` where `mask` holds are replaced by `src`'s, in one
# broadcast, on a vector or on the columns of a matrix or of an array whose
# columns are the second dimension
maskcolumns!(dst::AbstractVector, src, mask) = (dst .= ifelse.(mask, src, dst); dst)
maskcolumns!(dst::AbstractMatrix, src, mask) = (dst .= ifelse.(transpose(mask), src, dst); dst)
maskcolumns!(dst::AbstractArray{<:Any,3}, src, mask) = (dst .= ifelse.(reshape(mask, 1, :, 1), src, dst); dst)
scalecolumns(steps, a::AbstractVector) = steps .* a
scalecolumns(steps, a::AbstractMatrix) = transpose(steps) .* a
scalecolumns(steps, a::AbstractArray{<:Any,3}) = reshape(steps, 1, :, 1) .* a

function newtonsolve!(x, correction, trial, residual, trialresidual, baseresidual!,
        trialresidual!, refresh!, solve!, tol::AbstractVector, maxiters, stalefailed, nonlinear, iterative,
        work::NewtonWork; simplified::Bool = false, accept! = nothing)
    norms, trialnorms, previous, steps = work.norms, work.trialnorms, work.previous, work.steps
    converged = false
    fresh = false
    corrections, factorizations, retries, krylov = 0, 0, 0, 0
    laststale = false
    fill!(previous, Inf)
    adopted = false
    for iteration in 0:maxiters
        adopted || baseresidual!(norms, residual, x)
        adopted = false
        if all(j -> isfinite(norms[j]) && norms[j] <= tol[j], eachindex(norms))
            converged = true
            break
        end
        iteration == maxiters && break
        refreshnow = if simplified
            stalefailed || any(j -> norms[j] > tol[j] && norms[j] > 0.25*previous[j], eachindex(norms))
        else
            stalefailed || (!iterative && iteration >= 1)
        end
        if nonlinear && !fresh && refreshnow
            refresh!()
            factorizations += 1
            fresh = true
        end
        copyto!(previous, norms)
        stale, iterations = solve!(correction, residual)
        krylov += iterations
        laststale = stale
        corrections += 1
        # the line search, per column: an accepted column is taken into the
        # base point through the mask, a rejected one halves its step
        fill!(steps, 1.0)
        accepted, newly = work.accepted, work.newly
        fill!(accepted, false)
        for _ in 0:12
            # an accepted column has a zero step, so its trial is its base
            copyto!(work.stepsdev, steps)
            trial .= x .- scalecolumns(work.stepsdev, correction)
            trialresidual!(trialnorms, trialresidual, trial)
            fill!(newly, false)
            for j in eachindex(norms)
                accepted[j] && continue
                # a column already converged is accepted as it is
                if norms[j] <= tol[j]
                    accepted[j] = true
                elseif isfinite(trialnorms[j]) && (trialnorms[j] <= tol[j] || trialnorms[j] < (1 - 1e-4*steps[j])*norms[j])
                    accepted[j] = true
                    newly[j] = true
                end
                accepted[j] && (steps[j] = 0.0)
            end
            copyto!(work.mask, newly)
            maskcolumns!(x, trial, work.mask)
            if !isnothing(accept!)
                maskcolumns!(residual, trialresidual, work.mask)
                accept!(work.mask)
                for j in eachindex(norms)
                    newly[j] && (norms[j] = trialnorms[j])
                end
                adopted = true
            end
            all(accepted) && break
            if !fresh && nonlinear && !iterative
                # the retry: the Jacobian and the correction at the base
                # point, whose residual and cache the trial left alone
                refresh!()
                factorizations += 1
                retries += 1
                fresh = true
                solve!(correction, residual)
                corrections += 1
                continue
            end
            for j in eachindex(norms)
                accepted[j] || (steps[j] *= 0.5)
            end
            # a rejected column keeps its base point: its step is zero for
            # the accepted ones now, so their trial columns are untouched
        end
        all(accepted) || break
    end
    return converged, fresh, corrections, factorizations, retries, krylov, laststale
end

"""
    transientdemodulate(solution, port, frequency; quantity = :outgoing,
        window = t -> 1.0)

The complex peak amplitude of the saved trace of `port` at `frequency` in
Hz: the integral of `2*window(t)*trace(t)*exp(-2pi*im*frequency*t)` over
the integral of `window`, by trapezoidal quadrature on the saved samples.
`quantity` is `:voltage`, `:incident` or `:outgoing`. A smooth window
suppresses leakage from a strong pump; the saved rate must resolve the
carrier. On a device solution the port trace is downloaded once.
"""
function transientdemodulate(sol::TransientSolution, port::Integer, frequency::Real;
        quantity::Symbol = :outgoing, window = t -> 1.0)
    isfinite(frequency) || throw(ArgumentError("the frequency must be finite."))
    quantity in (:voltage, :incident, :outgoing) || throw(ArgumentError(lazy"quantity must be :voltage, :incident or :outgoing, not $(quantity)."))
    p = findfirst(port_ -> port_.number == port, sol.problem.circuit.ports)
    isnothing(p) && throw(ArgumentError(lazy"there is no port $(port)."))
    signal = Array(view(getproperty(sol, quantity), p, :))
    total, weight = 0.0im, 0.0
    tprev = sol.times[1]
    wprev = Float64(window(tprev))
    yprev = wprev*signal[1]*cispi(-2frequency*tprev)
    for k in 2:length(sol.times)
        t = sol.times[k]
        w = Float64(window(t))
        y = w*signal[k]*cispi(-2frequency*t)
        dt = t - tprev
        total += dt*(y + yprev)
        weight += dt*(w + wprev)/2
        tprev, wprev, yprev = t, w, y
    end
    (isfinite(weight) && weight > 0) || throw(ArgumentError("the window must have a positive finite integral."))
    return total/weight
end
