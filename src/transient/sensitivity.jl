# The tangent and the exact discrete adjoint of the steps taken: the
# derivative of a recorded trajectory on its own grid, linearized about the
# full loaded state, with respect to currents injected at the ports and to
# the initial state. The step matrix is symmetric, since the circuit
# matrices, the augmentation and the junction term are, so the adjoint
# solves use the factorization of the step itself.

# the outputs of a state as linear maps: the port voltage is `phi0*P*v`,
# a wave is `(V + s*Z*(I - G*V))/(2 sqrt(Z))` with `s = +1` incident and
# `-1` outgoing, so a wave is `cv .* V + cd .* I` per port
function outputcoefficients(sys::TransientSystem, quantity::Symbol)
    quantity in (:voltage, :incident, :outgoing) || throw(ArgumentError(
        lazy"quantity must be :voltage, :incident or :outgoing, not $(quantity)."))
    z, g = Array(sys.portimpedances), Array(sys.portconductances)
    np = length(z)
    cv, cd = ones(np), zeros(np)
    if quantity != :voltage
        s = quantity == :incident ? 1.0 : -1.0
        cv .= (1 .- s .* z .* g) ./ (2 .* sqrt.(z))
        cd .= s .* sqrt.(z) ./ 2
    end
    return tobackend(sys.backend, cv), tobackend(sys.backend, cd)
end

# the junction Jacobian applied to a direction, or to the columns of
# directions at once, `RJ' * (lmolj .* cos(phi) .* (RJ*d))`
function junctionproduct!(y, sys::TransientSystem, phi, d, work)
    stepmul!(work, sys.RJ, d)
    r = sys.relations
    if allsinusoidal(r)
        work .*= sys.lmolj .* cos.(phi)
    else
        derivativeinto!(sys.relationwork, r, phi)
        work .*= sys.lmolj .* sys.relationwork
    end
    stepmul!(y, sys.RJt, work)
    return y
end

function recordedsolution(sol::TransientSolution)
    sol.method isa WRspice && throw(ArgumentError(
        "the derivatives and the noise replay the package's own stepping rules; solve with Trapezoidal() or GaussLegendre() rather than WRspice()."))
    isnothing(sol.phases) && isnothing(sol.checkpoints) && throw(ArgumentError(
        "the derivatives read the junction phases at every step: solve with record = :phases, :states or :checkpoints and saveevery = 1."))
    isnothing(sol.phases) && !(sol.method isa GaussLegendre) && throw(ArgumentError("checkpoints are a record of the Gauss-Legendre rule."))
    length(sol.times) == sol.stats.steps + 1 || throw(ArgumentError("the complete trajectory is needed."))
    return nothing
end

"""
    transientinjection(problem, targets)

The unscaled injection of a unit current at each of `targets` into the node
equations, one sparse column per target: a port number injects into the
port's positive terminal, as a port source does, and a component name
injects the component's own current, out of its first terminal and into
its second, as a named source does, whatever the component is, which is
how a bath of a resistor is placed. The default targets of
[`transienttangent`](@ref) and [`transientadjoint`](@ref) are the ports
in compiled order.
"""
function transientinjection(p::TransientProblem, targets)
    psc = p.circuit
    rows, cols, vals = Int[], Int[], Float64[]
    for (k, target) in enumerate(targets)
        if target isa Integer
            q = findfirst(port -> port.number == target, psc.ports)
            isnothing(q) && throw(ArgumentError(lazy"there is no port $(target)."))
            n1, n2 = p.portpositive[q], p.portnegative[q]
        elseif target isa Union{AbstractString,Symbol}
            c = get(psc.componentnamedict, String(target), 0)
            c > 0 || throw(ArgumentError(lazy"there is no component $(target)."))
            n1, n2 = psc.nodeindices[2, c] - 1, psc.nodeindices[1, c] - 1
        else
            throw(ArgumentError(lazy"a target is a port number or a component name, not $(target)."))
        end
        n1 > 0 && (push!(rows, n1); push!(cols, k); push!(vals, 1.0))
        n2 > 0 && (push!(rows, n2); push!(cols, k); push!(vals, -1.0))
    end
    return sparse(rows, cols, vals, length(p), length(targets))
end

# the ports in compiled order, as targets
porttargets(p::TransientProblem) = [port.number for port in p.circuit.ports]

# the problem a plan or a bath is built against, from whatever names it
transientproblemof(p::TransientProblem) = p
transientproblemof(s::TransientSolution) = s.problem
transientproblemof(b::TransientBatchSolution) = first(b.problems)

# the compiled index of a port from its number, which is what the user
# facing functions take; the traces of a solution are in compiled order
function portindex(p::TransientProblem, number)
    number isa Integer || throw(ArgumentError(lazy"a port is named by its number, not $(number)."))
    q = findfirst(port -> port.number == number, p.circuit.ports)
    isnothing(q) && throw(ArgumentError(lazy"there is no port $(number)."))
    return q
end
portrows(p::TransientProblem, ports) = Int[portindex(p, number) for number in ports]

# which ports the targets are, zero for a component, for the direct
# feedthrough of a target's current into a port wave
targetports(p::TransientProblem, targets) = [t isa Integer ? findfirst(port -> port.number == t, p.circuit.ports) : 0 for t in targets]

"""
    ComponentPerturbation

The derivative of the scaled equations of a [`TransientProblem`](@ref) with
respect to a relative perturbation `p -> r*p` at `r = 1` of the value of
each of a set of named components, which the tangent and the adjoint of a
recorded transient carry as directions and objectives of their own (see
[`transientsensitivity`](@ref)). The equations are affine in `C`, `1/R`,
`1/L` and `1/Lj`, so the derivative is the component's own contribution to
each matrix, with a sign, as [`SensitivityStamp`](@ref) has it for the
linearized solve, built from the same classification
([`componentstamp`](@ref)) so that the two solvers support the same
components and reject the same ones. A component which is a port's own
termination moves the port's reference impedance and conductance with
it, as the linearized solve has it, which enters the port waves directly
(see [`directcoefficients`](@ref)).

# Fields
- `names`: the component names, in the order of the directions.
- `entries`, `hostentries`: the derivatives as their stored entries, one
    [`PerturbationEntries`](@ref) per kind on the backend and on the
    host, which the adjoint contracts its multipliers against; a
    component touches a few entries, so the contraction and its work are
    the size of the entries, not of the state times the components.
- `forcing`: for a tangent, the derivatives of the scaled capacitance,
    conductance, stiffness and junction current of every component
    stacked, `dC`, `dG`, `dL`, `dJ` on the backend with rows
    `(c - 1)n + 1` to `cn` holding component `c`'s, so one product gives
    every component's forcing at once, and `hdG`, `hdL`, `hdJ` on the
    host for the endpoint of a Gauss-Legendre step; the tangent's
    directions are dense right hand sides in any case. `nothing` for an
    adjoint, which never forms them.
- `states`: whether any component reads the state; the junctions alone
    read only the phases the record always holds.
- `ports`: the compiled port whose termination each component is, the
    row of its trace, or zero.
"""
struct ComponentPerturbation{E, H, F}
    names::Vector{String}
    n::Int
    entries::E
    hostentries::H
    forcing::F
    states::Bool
    ports::Vector{Int}
end

"""
    PerturbationEntries

The stored entries of one kind of derivative of the components, for the
contraction of an adjoint: the forcing of component `c` is `-dM_c A` for
the derivative `dM_c` and the state `A` its kind reads, so its
contraction against the multipliers `M` is the sum over the entries
`(i, j, v)` of `dM_c` of `-v A[j] M[i]`, taken as the gather of the rows
`j` of `A` and `i` of `M` by two selection matrices, their product, and
the sum into the components by a matrix carrying the values, every one
the size of the entries.
"""
struct PerturbationEntries{M}
    gatherstate::M
    gathermultiplier::M
    sum::M
    count::Int
end

# the stored entries of one kind of derivative: the row `i`, the column
# `j` and the value `v` of each, and the component `c` it belongs to, of
# a derivative with `k` columns
struct PerturbationTriplets
    i::Vector{Int}
    j::Vector{Int}
    v::Vector{Float64}
    c::Vector{Int}
    k::Int
end
PerturbationTriplets(k::Int) = PerturbationTriplets(Int[], Int[], Float64[], Int[], k)
function Base.push!(t::PerturbationTriplets, M::SparseMatrixCSC, c::Int)
    i, j, v = findnz(M)
    for q in eachindex(v)
        iszero(v[q]) && continue
        push!(t.i, i[q]); push!(t.j, j[q]); push!(t.v, v[q]); push!(t.c, c)
    end
    return t
end

# the entries of one kind of derivative, on the backend
function perturbationentries(t::PerturbationTriplets, n::Int, nc::Int, backend)
    m = length(t.v)
    d = A -> devicesparse(A, backend)
    return PerturbationEntries(d(sparse(1:m, t.j, ones(m), m, t.k)), d(sparse(1:m, t.i, ones(m), m, n)),
        d(sparse(t.c, 1:m, -t.v, nc, m)), m)
end

# the stacked derivative of one kind, `(ncomponents n, k)`, from its
# entries, for the forcing of a tangent
stackedderivative(t::PerturbationTriplets, n::Int, nc::Int) = sparse((t.c .- 1) .* n .+ t.i, t.j, t.v, nc*n, t.k)

# the component's index in the flat table from its name, however the
# table is keyed
function componentindex(psc::CompiledCircuit, name)
    key = String(name)
    idx = get(psc.componentnamedict, key, 0)
    iszero(idx) && (idx = get(psc.componentnamedict, Symbol(key), 0))
    iszero(idx) && throw(ArgumentError(lazy"The component $(name) is not in this circuit."))
    return idx
end

# the junction incidence of a problem and the junction coefficients
# `Lscale/Lj`, on the host, as the system builds them
function junctionincidence(p::TransientProblem)
    nm = p.matrices
    Ljb = nm.Ljb
    Rbnm = hcat(nm.Rbnm, spzeros(eltype(nm.Rbnm), size(nm.Rbnm, 1), p.Naux))
    RJ = SparseMatrixCSC{Float64,Int}(Rbnm[Ljb.nzind, :])
    return RJ, Float64[p.Lscale/Ljb.nzval[j] for j in eachindex(Ljb.nzval)]
end

# The perturbation of the named components of a problem, on the backend:
# the stored entries of every component's derivative, and, with
# `forcing`, the stacked derivatives a tangent's forcing is one product
# of; an adjoint asks for the entries alone, so nothing it holds grows
# with the state times the components. The entries are read here, once,
# from each component's stamp padded to the state, which no step builds.
function componentperturbation(p::TransientProblem, names, backend; forcing::Bool)
    psc, cg, nm = p.circuit, p.graph, p.matrices
    n, Nnodal, Lscale = length(p), p.Nnodal, p.Lscale
    Ljb = nm.Ljb
    nj = length(Ljb.nzval)
    RJ, lmolj = junctionincidence(p)
    RJt = sparse(transpose(RJ))
    lookups = componentlookups(p.coupledbranches, Ljb)
    isempty(names) && throw(ArgumentError("name at least one component."))
    nc = length(names)
    tC, tG, tL, tJ = PerturbationTriplets(n), PerturbationTriplets(n), PerturbationTriplets(n), PerturbationTriplets(nj)
    ports = zeros(Int, nc)
    for (c, name) in enumerate(names)
        idx = componentindex(psc, name)
        kind, info = componentstamp(idx, psc, cg, nm, lookups, 1, psc.Nnodes)
        pad = M -> mnapad(SparseMatrixCSC{Float64,Int}(real.(M)), n - Nnodal)
        # the capacitance enters as `r C`, the others inversely, so their
        # derivatives carry the minus of `d(1/(r p))/dr`
        if kind == :C
            push!(tC, Lscale .* pad(info), c)
        elseif kind == :G
            push!(tG, -Lscale .* pad(info), c)
        elseif kind == :invL
            push!(tL, -Lscale .* pad(info), c)
        else
            push!(tJ, sparse(findnz(RJt[:, info])[1], fill(info, nnz(RJt[:, info])), -lmolj[info] .* findnz(RJt[:, info])[2], n, nj), c)
        end
        # a port's own termination: the environments are listed by port
        # number, the traces by compiled port
        pos = findfirst(==(idx), nm.portenvironmentindices)
        isnothing(pos) || (ports[c] = portindex(p, nm.portnumbers[pos]))
    end
    states = !isempty(tC.v) || !isempty(tG.v) || !isempty(tL.v)
    entries = b -> (C = perturbationentries(tC, n, nc, b), G = perturbationentries(tG, n, nc, b),
        L = perturbationentries(tL, n, nc, b), J = perturbationentries(tJ, n, nc, b))
    stacked = if forcing
        d = A -> devicesparse(A, backend)
        dG, dL, dJ = stackedderivative(tG, n, nc), stackedderivative(tL, n, nc), stackedderivative(tJ, n, nc)
        (dC = d(stackedderivative(tC, n, nc)), dG = d(dG), dL = d(dL), dJ = d(dJ), hdG = dG, hdL = dL, hdJ = dJ)
    else
        nothing
    end
    return ComponentPerturbation(String.(collect(names)), n, entries(backend), entries(CPU()), stacked, states, ports)
end

# the quantities a step's residual reads, over `N` conditions: the
# acceleration, the rate, the flux and the relation values, their
# derivative's forcing when a tangent asks for one, and the state and the
# stage increments of a Gauss-Legendre step
function perturbationwork(cp::ComponentPerturbation, sys::TransientSystem, N::Int; forcing::Bool)
    allocate = (dims...) -> KernelAbstractions.zeros(sys.backend, Float64, dims...)
    nc, n = length(cp.names), cp.n
    nj = length(sys.lmolj)
    return (F = forcing ? allocate(nc*n, N) : nothing, work = forcing ? allocate(nc*n, N) : nothing,
        a = allocate(n, N), w = allocate(n, N), X = allocate(n, N), f = allocate(nj, N), fprev = allocate(nj, N),
        x = allocate(n, N), v = allocate(n, N), delta = allocate(n, N, 2))
end

# The forcing of the linearized equations by the perturbations, the
# negative of the derivative of the residual, `-(dC a + dG w + dL X + dJ f)`
# for the quantities the residual reads, each `(n, N)` over the
# conditions, into `F`, `(ncomponents n, N)`, whose reshaping to
# `(n, ncomponents N)` is the batch layout of the directions.
function perturbationforcing!(F, cp::ComponentPerturbation, pw, work)
    fill!(F, 0)
    st = cp.forcing
    if cp.states
        stepmul!(work, st.dC, pw.a); F .-= work
        stepmul!(work, st.dG, pw.w); F .-= work
        stepmul!(work, st.dL, pw.X); F .-= work
    end
    stepmul!(work, st.dJ, pw.f); F .-= work
    return F
end

# the work of the contraction of every kind of entries over `N`
# conditions and `nobj` objectives: the gathered state and multipliers,
# their product and the sum into the components
function contractionwork(cp::ComponentPerturbation, entries, N::Int, nobj::Int, backend)
    allocate = (dims...) -> KernelAbstractions.zeros(backend, Float64, dims...)
    nc = length(cp.names)
    one = e -> (GA = allocate(e.count, N), GM = allocate(e.count, nobj*N), prod = allocate(e.count, nobj), out = allocate(nc, nobj))
    return (C = one(entries.C), G = one(entries.G), L = one(entries.L), J = one(entries.J))
end

# `S[:, :, j] += weight * (the forcing of condition j)' * (the multipliers
# of condition j)` for one kind of entries, on the state `A`, `(rows, N)`,
# and the multipliers `M`, `(n, nobj N)`, the objectives of a condition
# contiguous in the columns
function contractentries!(S, e::PerturbationEntries, cw, A, M, weight)
    e.count == 0 && return S
    stepmul!(cw.GA, e.gatherstate, A)
    stepmul!(cw.GM, e.gathermultiplier, M)
    N, nobj = size(A, 2), size(cw.prod, 2)
    for j in 1:N
        cw.prod .= view(cw.GA, :, j) .* view(cw.GM, :, (j - 1)*nobj + 1:j*nobj)
        stepmul!(cw.out, e.sum, cw.prod)
        view(S, :, :, j) .+= weight .* cw.out
    end
    return S
end

# the contraction of a step's quantities against its multipliers, every
# kind of entries
function contractperturbation!(S, cp::ComponentPerturbation, entries, cw, pw, M, weight)
    if cp.states
        contractentries!(S, entries.C, cw.C, pw.a, M, weight)
        contractentries!(S, entries.G, cw.G, pw.w, M, weight)
        contractentries!(S, entries.L, cw.L, pw.X, M, weight)
    end
    contractentries!(S, entries.J, cw.J, pw.f, M, weight)
    return S
end

# the state at the start of the Gauss-Legendre step from time `k` to
# `k + 1` and the step's stage increments, from the window's record or
# replay
function stagestates!(pw, window, k::Int)
    window.state!(pw.x, pw.v, k)
    window.stages!(pw.delta, k)
    return nothing
end

# The quantities stage `i` of a Gauss-Legendre step reads: the stage
# residual's acceleration, rate and flux from the increments and the
# state at the start of the step, and the relation values at the stage
# phases `phi`.
function stagequantities!(pw, cp::ComponentPerturbation, sys::TransientSystem, gc::GaussCoefficients, phi, i::Int)
    if cp.states
        h = sys.h
        d1, d2 = stage(pw.delta, 1), stage(pw.delta, 2)
        pw.a .= (entry(gc.ainv2, i, 1) .* d1 .+ entry(gc.ainv2, i, 2) .* d2) ./ h^2 .- (gc.ainvone[i]/h) .* pw.v
        pw.w .= (entry(gc.ainv, i, 1) .* d1 .+ entry(gc.ainv, i, 2) .* d2) ./ h
        pw.X .= pw.x .+ stage(pw.delta, i)
    end
    relationinto!(pw.f, sys.relations, phi)
    return pw
end

# The quantities the step from time `k - 1` to `k` of the trapezoidal or
# backward Euler rule reads, from the recorded states and phases: the
# residual's terms in `x_{k-1}`, `x_k` and `v_{k-1}` differentiated, which
# for the trapezoidal rule is the rule's own second difference on the
# capacitance and both endpoints on the rest, and for backward Euler the
# endpoint.
function stepquantities!(pw, cp::ComponentPerturbation, sys::TransientSystem, sol, k::Int)
    h, alpha, beta = sys.h, sys.alpha, sys.beta
    trapezoidal = sys.method isa Trapezoidal
    if cp.states
        xk, xp = view(sol.flux, :, k), view(sol.flux, :, k - 1)
        vp = view(sol.rate, :, k - 1)
        if trapezoidal
            pw.a .= alpha .* (xk .- xp) .- (4/h) .* vp
            pw.w .= beta .* (xk .- xp)
            pw.X .= xk .+ xp
        else
            pw.a .= (xk .- xp) ./ h^2 .- vp ./ h
            pw.w .= (xk .- xp) ./ h
            pw.X .= xk
        end
    end
    relationinto!(pw.f, sys.relations, view(sol.phases, :, k))
    if trapezoidal
        relationinto!(pw.fprev, sys.relations, view(sol.phases, :, k - 1))
        pw.f .+= pw.fprev
    end
    return pw
end

# The work of the endpoint of a Gauss-Legendre step, whose projection and
# reading are host computations: the state on the backend and on the
# host, the junction phases and relation values, the host junction
# incidence and relation table, the reading's rows `Q`, and for a tangent
# the forcing `(n, ncomponents N)` with its stacked buffer, for an adjoint
# the contraction work.
function endpointwork(cp::ComponentPerturbation, sys::TransientSystem, p::TransientProblem, N::Int, nobj::Int; forcing::Bool)
    n, nj, nc = cp.n, length(sys.lmolj), length(cp.names)
    pr = sys.gauss.projection
    RJ, _ = junctionincidence(p)
    relations = isnothing(p.relations) ? emptyrelations(zeros(0)) : hostrelations(p.relations)
    Q = endpointinjection(pr, sparse(1.0I, n, n))
    return (F = forcing ? zeros(n, nc*N) : nothing, work = forcing ? zeros(nc*n, N) : nothing,
        cw = forcing ? nothing : contractionwork(cp, cp.hostentries, N, nobj, CPU()),
        x = KernelAbstractions.zeros(sys.backend, Float64, n, N), v = KernelAbstractions.zeros(sys.backend, Float64, n, N),
        xh = zeros(n, N), vh = zeros(n, N), phi = zeros(nj, N), f = zeros(nj, N), RJ, relations, Q)
end

# The quantities the endpoint of a Gauss-Legendre step at time `k` reads,
# on the host: the state from the window when the perturbation reads it,
# and the relation values of every junction at the endpoint; the
# junctions alone read the projected junctions' recorded phases `hphi`,
# the rows `pj`, which are the only ones the projection and the reading
# see.
function endpointquantities!(pe, cp::ComponentPerturbation, window, hphi, pj, k::Int)
    if cp.states
        window.state!(pe.x, pe.v, k)
        copyto!(pe.xh, pe.x)
        copyto!(pe.vh, pe.v)
        mul!(pe.phi, pe.RJ, pe.xh)
    else
        fill!(pe.phi, 0)
        pe.phi[pj, :] .= hphi
    end
    relationinto!(pe.f, pe.relations, pe.phi)
    return pe
end

# the forcing at the endpoint for a tangent, `-(dG v + dL x + dJ f)`, on
# the host
function endpointforcing!(pe, cp::ComponentPerturbation)
    n = cp.n
    fill!(pe.work, 0)
    if cp.states
        mul!(pe.work, cp.forcing.hdG, pe.vh)
        mul!(pe.work, cp.forcing.hdL, pe.xh, 1.0, 1.0)
    end
    mul!(pe.work, cp.forcing.hdJ, pe.f, 1.0, 1.0)
    pe.F .= .-reshape(pe.work, n, :)
    return pe.F
end

# the contraction at the endpoint for an adjoint, against the host
# multipliers `M`, `(n, nobj N)`
function endpointcontract!(S, pe, cp::ComponentPerturbation, M, weight)
    e, cw = cp.hostentries, pe.cw
    if cp.states
        contractentries!(S, e.G, cw.G, pe.vh, M, weight)
        contractentries!(S, e.L, cw.L, pe.xh, M, weight)
    end
    contractentries!(S, e.J, cw.J, pe.f, M, weight)
    return S
end

"""
    directcoefficients(sys, quantity)

The derivative of the port wave coefficients of `outputcoefficients`
with respect to a relative perturbation of a port's own termination,
which moves the port's reference impedance and its termination
conductance together, `z -> r z` and `g -> g/r`, as the linearized solve
moves them: `-cv/2` and `cd/2` for a wave, since `z g` is unchanged, and
nothing for the voltage.
"""
function directcoefficients(sys::TransientSystem, quantity::Symbol)
    cv, cd = outputcoefficients(sys, quantity)
    quantity == :voltage && return zero(cv), zero(cd)
    return -cv ./ 2, cd ./ 2
end

# the current the drives of every condition push into each port at time
# `t`, on the host, `(nports, N)`
function portdrivecurrents!(I, problems, t)
    fill!(I, 0)
    for (j, q) in enumerate(problems), d in q.drives
        d.portindex > 0 && (I[d.portindex, j] += d.current(t))
    end
    return I
end

# The direct term of the port waves of a tangent along the components
# which are port terminations: at time `k`, for each such component `c`
# of port `q`, the derivative of the wave coefficients times the recorded
# port voltage and the drive current of each condition, into the
# `(nports, ndir N)` outputs of each quantity through the host matrices
# `direct`, one per quantity.
struct DirectTerm{V, B}
    ports::Vector{Int}
    dcv::Dict{Symbol, Vector{Float64}}
    dcd::Dict{Symbol, Vector{Float64}}
    I::Matrix{Float64}
    V::Matrix{Float64}
    direct::Matrix{Float64}
    dev::V
    backend::B
end
function directterm(cp::ComponentPerturbation, sys::TransientSystem, N::Int, ndir::Int)
    any(>(0), cp.ports) || return nothing
    np = length(sys.portimpedances)
    dcv = Dict{Symbol, Vector{Float64}}(); dcd = Dict{Symbol, Vector{Float64}}()
    for q in (:voltage, :incident, :outgoing)
        a, b = directcoefficients(sys, q)
        dcv[q], dcd[q] = Array(a), Array(b)
    end
    return DirectTerm(cp.ports, dcv, dcd, zeros(np, N), zeros(np, N), zeros(np, ndir*N),
        KernelAbstractions.zeros(sys.backend, Float64, np, ndir*N), sys.backend)
end
# the term of time `k` added to the outputs `outs` of the three quantities
function adddirectterm!(outs, dt::DirectTerm, sol, problems, k::Int)
    copyto!(dt.V, view(sol.voltage, :, k, :))
    portdrivecurrents!(dt.I, problems, sol.times[k])
    N, ndir = size(dt.I, 2), length(dt.ports)
    for (s, q) in enumerate((:voltage, :incident, :outgoing))
        q == :voltage && continue
        fill!(dt.direct, 0)
        for j in 1:N, (c, port) in enumerate(dt.ports)
            port > 0 || continue
            dt.direct[port, (j - 1)*ndir + c] = dt.dcv[q][port]*dt.V[port, j] + dt.dcd[q][port]*dt.I[port, j]
        end
        copyto!(dt.dev, dt.direct)
        outs[s] .+= dt.dev
    end
    return nothing
end

# The direct term of an adjoint's objective along the port termination
# components: the weights of the objective on the port waves times the
# derivative of the wave coefficients times the recorded port voltage and
# the drive current, summed over the record, into `S`, `(nc, nobj, N)` on
# the host.
function adddirectsensitivity!(S, cp::ComponentPerturbation, sys::TransientSystem, sol, problems, weights, quantity::Symbol)
    (quantity != :voltage && any(>(0), cp.ports)) || return S
    dcv, dcd = Array.(directcoefficients(sys, quantity))
    np, nt = size(weights, 1), size(weights, 2)
    N, nobj = length(problems), (ndims(weights) == 3 ? size(weights, 3) : 1)
    w = reshape(Float64.(collect(weights)), np, nt, nobj)
    V = reshape(Array(sol.voltage), np, nt, N)
    I = zeros(np, N)
    for k in 1:nt
        portdrivecurrents!(I, problems, sol.times[k])
        for j in 1:N, (c, port) in enumerate(cp.ports)
            port > 0 || continue
            d = dcv[port]*V[port, k, j] + dcd[port]*I[port, j]
            for o in 1:nobj
                S[c, o, j] += w[port, k, o]*d
            end
        end
    end
    return S
end

# `X = F \ B` on every column at once where the factorization solves
# matrices, KLU here and cuDSS in its extension, and column by column on
# one that solves vectors
function matrixsolve!(X, factor, B)
    for j in axes(B, 2)
        myldiv!(view(X, :, j), factor, view(B, :, j))
    end
    return X
end
matrixsolve!(X, factor::KLU.KLUFactorization, B) = myldiv!(X, factor, B)
matrixsolve!(X, factor::Transpose{<:Any,<:KLU.KLUFactorization}, B) = myldiv!(X, factor, B)

"""
    transienttangent(solution, currents; targets = the ports,
        initialstate = nothing, factorization = nothing, reuse = nothing,
        outputsink = nothing)

The tangent of a recorded transient along a perturbation: `currents[q, k]`
is an additional Norton current in Amperes at target `q` (a port number,
or a component name, see [`transientinjection`](@ref)) at recorded time
`k`, and `initialstate` an optional pair of perturbations of the scaled
flux and rate at the start. A third dimension of `currents` is a set of
directions propagated together, each step's factorization serving them
all. Under [`GaussLegendre`](@ref) a current on the grid is read at the
stage times through a cubic Lagrange stencil; a current given as
`currents[q, s, k, direction]`, with `s = 1` its value at recorded time
`k` and `s = 2, 3` its values at the two stage times of the step from
`k` to `k + 1`, is read as it is, so a pulse keeps its support and a
tone its exact phase at the stages; the trapezoidal rule reads the grid
values of either form. Returns `(voltage, incident, outgoing,
finalflux, finalrate)` in the units of the solve, on its backend, with
the directions as the trailing dimension; with an `outputsink`, a
function `outputsink(k, voltage, incident, outgoing)` receiving the
three port by direction matrices of recorded time `k` on the backend,
valid until the next call, the histories are not stored and those three
are `nothing`, so a measurement of a long record needs no memory per
time. The linearization is about the full recorded state, so the loaded
junction phases enter every response; both the trajectory and its grid
are held fixed.
"""
function transienttangent(sol::TransientSolution, currents::Union{Nothing,AbstractArray{<:Real}};
        targets = porttargets(sol.problem), initialstate = nothing,
        factorization = nothing, reuse = nothing, outputsink = nothing, perturbation = nothing)
    recordedsolution(sol)
    p = sol.problem
    backend = KernelAbstractions.get_backend(sol.finalflux)
    fact = isnothing(factorization) ? transientfactorization(backend) : factorization
    sys = transientsystem(reuse, p, sol.dt, sol.method, backend, fact)
    n, np, nt = length(p), length(p.portimpedances), length(sol.times)
    nq = length(targets)
    ndir = tangentdirections(currents, perturbation, nq, nt)
    # invoked dynamically on the untyped kept system (see transientsolve)
    sys.method isa GaussLegendre && return Base.invokelatest(gaussbatchtangent, sol, currents, targets, initialstate, sys, outputsink, perturbation)
    h = sys.h
    trapezoidal = sys.method isa Trapezoidal
    # the components' forcing of each step, from the recorded states, and
    # the direct term of the port waves along a port's own termination
    pwork = isnothing(perturbation) ? nothing : perturbationwork(perturbation, sys, 1; forcing = true)
    dterm = isnothing(perturbation) ? nothing : directterm(perturbation, sys, 1, ndir)
    injection = devicesparse((sys.Lscale/phi0) .* transientinjection(p, targets), backend)
    ports = tobackend(backend, [q > 0 ? q : 0 for q in targetports(p, targets)])
    # the trapezoidal step reads the grid values of a staged current; a
    # tangent along the components alone carries one zero column of them
    dI = isnothing(currents) ? KernelAbstractions.zeros(backend, Float64, nq, 1, ndir) :
        tobackend(backend, reshape(Float64.(collect(ndims(currents) == 4 ? selectdim(currents, 2, 1) : currents)), nq, nt, ndir))
    kcol = isnothing(currents) ? (k -> 1) : identity
    allocate = (dims...) -> KernelAbstractions.zeros(backend, Float64, dims...)
    dx, dv, dxnew, rhs, work = allocate(n, ndir), allocate(n, ndir), allocate(n, ndir), allocate(n, ndir), allocate(n, ndir)
    if !isnothing(initialstate)
        x0, v0 = initialstate
        (size(x0, 1) == n && size(v0, 1) == n) || throw(DimensionMismatch(lazy"the state has $(n) entries."))
        copyto!(dx, reshape(Float64.(collect(x0)), n, :))
        copyto!(dv, reshape(Float64.(collect(v0)), n, :))
    end
    nj = length(sys.lmolj)
    phi, jwork = allocate(nj), allocate(nj, ndir)
    voltage, incident, outgoing = isnothing(outputsink) ? [allocate(np, nt, ndir) for _ in 1:3] : (nothing, nothing, nothing)
    coefficients = Dict(q => outputcoefficients(sys, q) for q in (:voltage, :incident, :outgoing))
    portwork = allocate(np, ndir)
    # the current of the targets that are ports, into their port waves
    portmap = devicesparse(sparse([q for q in targetports(p, targets) if q > 0],
        [k for (k, q) in enumerate(targetports(p, targets)) if q > 0], ones(count(>(0), targetports(p, targets))), np, nq), backend)
    directwork = allocate(np, ndir)
    outwork = isnothing(outputsink) ? nothing : [allocate(np, ndir) for _ in 1:3]
    function outputs!(k)
        stepmul!(portwork, sys.ports, dv)
        portwork .*= phi0
        stepmul!(directwork, portmap, view(dI, :, kcol(k), :))
        outs = isnothing(outputsink) ? (view(voltage, :, k, :), view(incident, :, k, :), view(outgoing, :, k, :)) : outwork
        for (s, q) in enumerate((:voltage, :incident, :outgoing))
            cv, cd = coefficients[q]
            outs[s] .= cv .* portwork .+ cd .* directwork
        end
        isnothing(dterm) || adddirectterm!(outs, dterm, sol, [p], k)
        isnothing(outputsink) || outputsink(k, outwork[1], outwork[2], outwork[3])
        nothing
    end
    outputs!(1)
    factor = isnothing(reuse) ? nothing : reuse.factor
    for k in 2:nt
        # the right hand side of the tangent step, from the previous
        # tangent and the current perturbations: for the trapezoidal rule
        #   [alpha*C + beta*G - L - J'(x_n)] dx_n + (4/h) C dv_n + db_{n+1} + db_n,
        # for backward Euler
        #   C (dx_n/h^2 + dv_n/h) + G dx_n/h + db_{n+1}
        stepmul!(rhs, sys.A, dx)
        stepmul!(work, sys.B, dv); rhs .+= work
        if trapezoidal
            junctionproduct!(work, sys, view(sol.phases, :, k - 1), dx, jwork); rhs .-= work
            stepmul!(work, injection, view(dI, :, kcol(k - 1), :)); rhs .+= work
        end
        stepmul!(work, injection, view(dI, :, kcol(k), :)); rhs .+= work
        if !isnothing(perturbation)
            stepquantities!(pwork, perturbation, sys, sol, k)
            perturbationforcing!(pwork.F, perturbation, pwork, pwork.work)
            rhs .+= reshape(pwork.F, n, ndir)
        end
        # the step matrix at the recorded phases, and the solves
        copyto!(phi, view(sol.phases, :, k))
        factor = stepjacobian!(sys, phi, factor)
        matrixsolve!(dxnew, factor, rhs)
        if trapezoidal
            dv .= (2/h) .* (dxnew .- dx) .- dv
        else
            dv .= (dxnew .- dx) ./ h
        end
        copyto!(dx, dxnew)
        outputs!(k)
    end
    KernelAbstractions.synchronize(backend)
    isnothing(reuse) || (reuse.factor = factor)
    squeeze = a -> isnothing(a) ? nothing : (!isnothing(currents) && ndims(currents) == 2) ? reshape(a, size(a)[1:end-1]...) : a
    return (; voltage = squeeze(voltage), incident = squeeze(incident), outgoing = squeeze(outgoing),
        finalflux = squeeze(copy(dx)), finalrate = squeeze(copy(dv)))
end

# the directions of a tangent: those of its currents, or, along the
# components alone with no currents, one per component
function tangentdirections(currents, perturbation, nq, nt)
    if isnothing(currents)
        isnothing(perturbation) && throw(ArgumentError("a tangent needs currents, or components to differentiate along."))
        return length(perturbation.names)
    end
    ndir = currentshape(currents, nq, nt)
    isnothing(perturbation) || ndir == length(perturbation.names) || throw(DimensionMismatch(
        "the perturbation gives one direction per component."))
    return ndir
end

# The columns of an adjoint's currents at the recorded times are final a
# few steps after a step first touches them, the trapezoidal step touching
# its own two endpoints and the Gauss-Legendre stencil up to four grid
# points: a ring of five columns holds the pending ones, and a column is
# handed to the sink, with the direct feedthrough of a port's current into
# its own wave added, once no remaining step touches it, in decreasing
# time. A sink that stores every column is the default; the noise sinks
# each column into its bath contraction and stores none.
struct CurrentRing{A, M, S, F}
    slots::A
    column::M
    sink::S
    feedthrough!::F
end
function CurrentRing(backend, nq, nobj, sink, feedthrough!)
    slots = KernelAbstractions.zeros(backend, Float64, nq, nobj, 5)
    column = KernelAbstractions.zeros(backend, Float64, nq, nobj)
    return CurrentRing(slots, column, sink, feedthrough!)
end
ringadd!(ring::CurrentRing, j, values, weight) = (view(ring.slots, :, :, mod1(j, 5)) .+= weight .* values; nothing)
function ringemit!(ring::CurrentRing, j)
    slot = view(ring.slots, :, :, mod1(j, 5))
    ring.column .= slot
    ring.feedthrough!(ring.column, j)
    ring.sink(j, ring.column)
    fill!(slot, 0)
    return nothing
end

# the pieces of an adjoint every rule shares: the transposed injection and
# port map, the output coefficients, the feedthrough, and the ring with the
# storing sink when none is given
function adjointsetup(sol::TransientSolution, weights, quantity::Symbol, targets, sys::TransientSystem, sink)
    p = sol.problem
    backend = sys.backend
    n, np, nt = length(p), length(p.portimpedances), length(sol.times)
    nq, nobj = length(targets), (ndims(weights) == 3 ? size(weights, 3) : 1)
    w = tobackend(backend, reshape(Float64.(collect(weights)), np, nt, nobj))
    cv, cd = outputcoefficients(sys, quantity)
    injectiont = devicesparse(sparse(transpose((sys.Lscale/phi0) .* transientinjection(p, targets))), backend)
    tp = targetports(p, targets)
    portmapt = devicesparse(sparse([k for (k, q) in enumerate(tp) if q > 0], [q for q in tp if q > 0],
        ones(count(>(0), tp)), nq, np), backend)
    allocate = (dims...) -> KernelAbstractions.zeros(backend, Float64, dims...)
    directwork, targetwork = allocate(np, nobj), allocate(nq, nobj)
    feedthrough! = (column, k) -> begin
        directwork .= cd .* view(w, :, k, :)
        stepmul!(targetwork, portmapt, directwork)
        column .+= targetwork
        nothing
    end
    currents = isnothing(sink) ? allocate(nq, nt, nobj) : nothing
    store = isnothing(sink) ? (k, values) -> (copyto!(view(currents, :, k, :), values); nothing) : sink
    ring = CurrentRing(backend, nq, nobj, store, feedthrough!)
    return w, cv, injectiont, ring, currents
end

# the shape of a tangent's currents, on the grid or at the stages, and
# the number of directions
function currentshape(currents, nq, nt)
    ongrid = ndims(currents) in (2, 3) && size(currents, 1) == nq && size(currents, 2) == nt
    staged = ndims(currents) == 4 && size(currents, 1) == nq && size(currents, 2) == 3 && size(currents, 3) == nt
    ongrid || staged || throw(DimensionMismatch(
        lazy"currents must have one row per target ($(nq)), one column per recorded time ($(nt)) and optionally a third dimension of directions, or be (targets, 3, times, directions) at the grid and the stages."))
    all(isfinite, currents) || throw(ArgumentError("the current perturbations must be finite."))
    return ndims(currents) == 2 ? 1 : size(currents, ndims(currents))
end

"""
    transientadjoint(solution, weights; quantity = :outgoing, targets = the ports,
        factorization = nothing, reuse = nothing, sink = nothing, stagesink = nothing)

The exact discrete adjoint of `sum(weights .* getproperty(solution,
quantity))` on the recorded grid, `weights` having one row per port and
one column per recorded time and `quantity` being `:voltage`, `:incident`
or `:outgoing`; a third dimension of `weights` is a set of objectives
propagated together on each step's factorization. Returns `(currents,
initialflux, initialrate)`: the derivatives with respect to a Norton
current at every target (a port number or a component name, see
[`transientinjection`](@ref)) and recorded time, in Amperes, including
the direct feedthrough of a port's current into its wave, and with
respect to the scaled initial flux and rate, with the objectives as the
trailing dimension. For a consistent perturbation the contraction equals
the weighted output of [`transienttangent`](@ref); a complex demodulation
is two real objectives, and a time integral carries its quadrature
weights in `weights`. With a `sink`, a function `sink(k, values)`, the
currents are not stored: each column, a targets by objectives matrix on
the backend valid until the next call, is handed to the sink once it is
final, in decreasing recorded time, and `currents` is `nothing`; a
contraction over a long record then needs no memory per time. With a
`stagesink` as well, a function `stagesink(k, i, values)`, the
multipliers of the two stages of each Gauss-Legendre step from `k` to
`k + 1` go to it as they are, at their stage times, and the grid columns
hold only what the grid reads; that pair is the transpose of the staged
form of the tangent's currents, and it is what the noise contracts.
"""
function transientadjoint(sol::TransientSolution, weights::AbstractArray{<:Real};
        quantity::Symbol = :outgoing, targets = porttargets(sol.problem),
        factorization = nothing, reuse = nothing, sink = nothing, stagesink = nothing,
        components = String[])
    recordedsolution(sol)
    p = sol.problem
    backend = KernelAbstractions.get_backend(sol.finalflux)
    fact = isnothing(factorization) ? transientfactorization(backend) : factorization
    sys = transientsystem(reuse, p, sol.dt, sol.method, backend, fact)
    n, np, nt = length(p), length(p.portimpedances), length(sol.times)
    nq, nobj = length(targets), (ndims(weights) == 3 ? size(weights, 3) : 1)
    ndims(weights) in (2, 3) && size(weights, 1) == np && size(weights, 2) == nt || throw(DimensionMismatch(
        lazy"weights must have one row per port ($(np)), one column per recorded time ($(nt)) and optionally a third dimension of objectives."))
    all(isfinite, weights) || throw(ArgumentError("the weights must be finite."))
    perturbation = isempty(components) ? nothing : componentperturbation(p, components, backend; forcing = false)
    isnothing(perturbation) || recordedstates(sol, perturbation)
    # invoked dynamically on the untyped kept system (see transientsolve)
    sys.method isa GaussLegendre && return Base.invokelatest(gaussbatchadjoint, sol, weights, quantity, targets, sys, sink, stagesink, perturbation)
    h = sys.h
    trapezoidal = sys.method isa Trapezoidal
    w, cv, injectiont, ring, currents = adjointsetup(sol, weights, quantity, targets, sys, sink)
    allocate = (dims...) -> KernelAbstractions.zeros(backend, Float64, dims...)
    # the components' forcing of each step, contracted against the step's
    # multipliers on the entries of the components, and the direct term of
    # the objective along a port's own termination
    pwork = isnothing(perturbation) ? nothing : perturbationwork(perturbation, sys, 1; forcing = false)
    pcwork = isnothing(perturbation) ? nothing : contractionwork(perturbation, perturbation.entries, 1, nobj, backend)
    sensitivity = isnothing(perturbation) ? nothing : allocate(length(perturbation.names), nobj, 1)
    isnothing(perturbation) || (sensitivity .+= tobackend(backend,
        adddirectsensitivity!(zeros(length(perturbation.names), nobj, 1), perturbation, sys, sol, [p], weights, quantity)))
    xbar, vbar, lambda, work = allocate(n, nobj), allocate(n, nobj), allocate(n, nobj), allocate(n, nobj)
    nj = length(sys.lmolj)
    phi, jwork = allocate(nj), allocate(nj, nobj)
    portwork = allocate(np, nobj)
    targetwork = allocate(nq, nobj)
    # the adjoint of an output at time k: the rate receives phi0*P'*(cv.*w)
    function output!(k)
        portwork .= cv .* view(w, :, k, :)
        stepmul!(work, sys.portst, portwork)
        vbar .+= phi0 .* work
    end
    output!(nt)
    factor = isnothing(reuse) ? nothing : reuse.factor
    for k in nt:-1:2
        # the rate update feeds the flux adjoint, then the step's solve is
        # transposed on the symmetric step matrix at the recorded phases
        xbar .+= (trapezoidal ? 2/h : 1/h) .* vbar
        copyto!(phi, view(sol.phases, :, k))
        factor = stepjacobian!(sys, phi, factor)
        matrixsolve!(lambda, factor, xbar)
        # the components' forcing of the step against its multipliers
        if !isnothing(perturbation)
            stepquantities!(pwork, perturbation, sys, sol, k)
            contractperturbation!(sensitivity, perturbation, perturbation.entries, pcwork, pwork, lambda, 1.0)
        end
        # the current perturbations the step read; the column at k is
        # final once this step has added to it
        stepmul!(targetwork, injectiont, lambda)
        ringadd!(ring, k, targetwork, 1.0)
        trapezoidal && ringadd!(ring, k - 1, targetwork, 1.0)
        ringemit!(ring, k)
        # the previous state's adjoints; A and B are symmetric
        stepmul!(xbar, sys.A, lambda)
        if trapezoidal
            junctionproduct!(work, sys, view(sol.phases, :, k - 1), lambda, jwork); xbar .-= work
            xbar .-= (2/h) .* vbar
            stepmul!(work, sys.B, lambda)
            vbar .= work .- vbar
        else
            xbar .-= (1/h) .* vbar
            stepmul!(vbar, sys.B, lambda)
        end
        output!(k - 1)
    end
    ringemit!(ring, 1)
    KernelAbstractions.synchronize(backend)
    isnothing(reuse) || (reuse.factor = factor)
    squeeze = a -> isnothing(a) ? nothing : ndims(weights) == 2 ? reshape(a, size(a)[1:end-1]...) : a
    return (; currents = squeeze(currents), initialflux = squeeze(xbar), initialrate = squeeze(vbar), initialwaves = nothing, initialstates = nothing,
        sensitivity = isnothing(sensitivity) ? nothing : squeeze(reshape(sensitivity, size(sensitivity, 1), nobj)))
end

# what a perturbation reads of the record: the junction phases every
# record holds, and for a component other than a junction the states, or
# the checkpoints the Gauss-Legendre rule replays them from
function recordedstates(sol, cp::ComponentPerturbation)
    cp.states || return nothing
    (isnothing(sol.flux) && isnothing(sol.checkpoints)) && throw(ArgumentError(
        "the sensitivity to a capacitor, inductor or resistor reads the recorded states: solve with record = :states, or with record = :checkpoints under GaussLegendre(); a junction alone reads the phases."))
    return nothing
end

"""
    transientsensitivity(solution, names; factorization = nothing,
        reuse = nothing, outputsink = nothing)

The derivative of the port responses of a recorded transient, or of every
condition of a [`TransientBatchSolution`](@ref) on one pass, with respect
to a relative perturbation `r` of the value of each component `names`
names (`p -> r*p` at `r = 1`), as `Ssensitivity` of [`hblinsolve`](@ref)
is of the scattering parameters. Returns `(voltage, incident, outgoing,
finalflux, finalrate)` as [`transienttangent`](@ref) does, with the
components as the trailing dimension, before the conditions of a batch.
The components supported are those of the linearized solve, `C`, `L`, `R`
and `Lj` with numeric values and no mutual coupling. The equations of a
step are differentiated as they were taken, so the derivative is exact
for the recorded steps; a capacitor, inductor or resistor reads the
recorded states, `record = :states`, or `record = :checkpoints` under
[`GaussLegendre`](@ref), which replays them, while a junction reads the
phases every record holds. The adjoint counterpart is the `components`
keyword of [`transientadjoint`](@ref), whose `sensitivity` is the
derivative of its objective with respect to the same perturbations.
"""
function transientsensitivity(sol::Union{TransientSolution,TransientBatchSolution}, names;
        factorization = nothing, reuse = nothing, outputsink = nothing)
    recordedsolution(sol)
    p = transientproblemof(sol)
    backend = KernelAbstractions.get_backend(sol.finalflux)
    cp = componentperturbation(p, names, backend; forcing = true)
    recordedstates(sol, cp)
    return transienttangent(sol, nothing; factorization, reuse, outputsink, perturbation = cp)
end
