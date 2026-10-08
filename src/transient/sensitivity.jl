# The tangent and the exact discrete adjoint of the steps taken: the
# derivative of a recorded trajectory on its own grid, linearized about the
# full loaded state, with respect to currents injected at the ports and to
# the initial state; here the entries of the responses and the steps of
# the trapezoidal and backward Euler rules, and the perturbation of the
# component values, which the Gauss-Legendre responses of batch.jl read as
# well. Under those two rules, which step no scattering block, the step
# matrix is symmetric, since the circuit matrices, the augmentation and
# the junction term are, so the adjoint solves use the factorization of
# the step itself.

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

# the system a response of a recorded solution or batch `s` of the problem
# `p` steps on, the reuse's kept one where it serves
function responsesystem(s, p::TransientProblem, factorization, reuse)
    recordedsolution(s)
    backend = KernelAbstractions.get_backend(s.finalflux)
    return transientsystem(reuse, p, s.dt, s.method, backend, steppingfactorization(factorization, backend))
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
port's first (positive) terminal, as a port source does, and a component
name drives the current through the component from its first terminal to
its second, as a named `CurrentSource` does, whatever the component is:
it draws the current from the node at its first terminal and delivers it
to the node at its second. The default targets of
[`transienttangent`](@ref) and [`transientadjoint`](@ref) are the ports
in the order of their numbers.
"""
function transientinjection(p::TransientProblem, targets)
    psc = p.circuit
    rows, cols, vals = Int[], Int[], Float64[]
    for (k, target) in enumerate(targets)
        if target isa Integer
            q = portindex(p, target)
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

# the ports in the order of their numbers, as targets
porttargets(p::TransientProblem) = [port.number for port in p.ports]

# The targets as the edges of a multigraph whose vertices are the floating
# subnetworks of `p` and the grounded rest of the circuit, 1 the grounded
# rest and 1 + r subnetwork r: a unit current of target `k` enters vertex
# `ends[1, k]` and leaves `ends[2, k]`, as the rows of its injection in the
# problem's units (`targetinjection`) give them, a terminal at ground
# holding none. A target whose terminals an element joins, a port with its
# termination or any component but a current source, has equal ends; one
# whose ends differ drives a net current into a subnetwork, which has no
# path back but through the gauge row, and a response to it is set by
# which node of the subnetwork that row sits on (see `bindsources`).
function targetedges(p::TransientProblem, injection::SparseMatrixCSC)
    vertex = ones(Int, size(injection, 1))
    for (r, island) in enumerate(p.floatingcomponents), i in island
        vertex[i - 1] = 1 + r
    end
    ends = ones(Int, 2, size(injection, 2))
    rows, vals = rowvals(injection), nonzeros(injection)
    for k in axes(injection, 2), t in nzrange(injection, k)
        iszero(vals[t]) || (ends[vals[t] > 0 ? 1 : 2, k] = vertex[rows[t]])
    end
    return ends
end

# The currents of a tangent in the form the steps read them
# (`TangentCurrents`), judged as the solve judges its drives
# (`checkbalance`): their net into each floating subnetwork at every
# recorded time and direction, and under Gauss-Legendre at the stages of
# every step, refused beyond the rounding of its sum. The trapezoidal rule
# reads the grid alone, and the stages a staged current gives at the last
# time begin no step. Along a net current the tangent would differentiate
# a solve that refuses it.
function checktangentbalance(p::TransientProblem, injection::SparseMatrixCSC, currents::TangentCurrents,
        method::AbstractTransientIntegrator)
    (currents.constant || isempty(p.floatingcomponents)) && return nothing
    ends = targetedges(p, injection)
    all(k -> ends[1, k] == ends[2, k], axes(ends, 2)) && return nothing
    terms = 2size(ends, 2)
    net, scale = zeros(1 + length(p.floatingcomponents)), zeros(1 + length(p.floatingcomponents))
    grid, stages = currents.grid, currents.stages
    if !isempty(currents.targets)
        # each direction's waveform at its target, judged in a column of
        # every target as the values of every target are
        column = zeros(size(ends, 2))
        for d in axes(grid, 3), k in axes(grid, 2)
            r = unbalancedtarget!(net, scale, ends, terms, column, currents.targets[d], grid[1, k, d])
            r > 0 && unbalancedtangent(p, currents, r, net[1 + r], 1, k, d)
        end
        if method isa GaussLegendre
            for d in axes(stages, 4), k in 1:size(stages, 3) - 1, i in 1:2
                r = unbalancedtarget!(net, scale, ends, terms, column, currents.targets[d], stages[1, i, k, d])
                r > 0 && unbalancedtangent(p, currents, r, net[1 + r], i + 1, k, d)
            end
        end
        return nothing
    end
    for d in axes(grid, 3), k in axes(grid, 2)
        r = unbalancedsubnetwork!(net, scale, ends, terms, grid, (k, d))
        r > 0 && unbalancedtangent(p, currents, r, net[1 + r], 1, k, d)
    end
    if method isa GaussLegendre
        for d in axes(stages, 4), k in 1:size(stages, 3) - 1, i in 1:2
            r = unbalancedsubnetwork!(net, scale, ends, terms, stages, (i, k, d))
            r > 0 && unbalancedtangent(p, currents, r, net[1 + r], i + 1, k, d)
        end
    end
    return nothing
end

# the judgment of one direction's `value` at its target `q` alone, in the
# work `column` of every target, which it leaves zero
function unbalancedtarget!(net::Vector{Float64}, scale::Vector{Float64}, ends::Matrix{Int}, terms::Int,
        column::Vector{Float64}, q::Int, value::Float64)
    column[q] = value
    r = unbalancedsubnetwork!(net, scale, ends, terms, column, ())
    column[q] = 0.0
    return r
end

# The net currents of the targets' currents `c[:, I...]` into the vertices
# of their graph (`targetedges`) and the sums of their terms' magnitudes,
# a target within one vertex adding to neither: the first floating
# subnetwork whose net is beyond the rounding of its `terms`, or zero.
function unbalancedsubnetwork!(net::Vector{Float64}, scale::Vector{Float64}, ends::Matrix{Int}, terms::Int,
        c::Array{Float64}, I::NTuple{N, Int}) where {N}
    fill!(net, 0.0)
    fill!(scale, 0.0)
    for q in axes(ends, 2)
        u, v = ends[1, q], ends[2, q]
        u == v && continue
        x = c[q, I...]
        net[u] += x
        net[v] -= x
        scale[u] += abs(x)
        scale[v] += abs(x)
    end
    for r in 2:length(net)
        abs(net[r]) > terms*eps(Float64)*scale[r] && return r - 1
    end
    return 0
end

# the refusal of a tangent's net current `net` into floating subnetwork
# `r` at stage `s` (1 the grid value), recorded time `k` and direction
# `d`, named by its index in the array the caller gave
@noinline function unbalancedtangent(p::TransientProblem, currents::TangentCurrents, r::Int, net::Float64,
        s::Int, k::Int, d::Int)
    nodes = join(p.circuit.nodenames[p.floatingcomponents[r]], ", ")
    index = currents.staged ? "$(s), $(k), $(d)" : currents.single ? "$(k)" : "$(k), $(d)"
    throw(ArgumentError(lazy"the currents of the tangent drive a net current of $(net) A into the nodes ($(nodes)) at currents[:, $(index)], which no element connects to ground, so the current has no path back; perturb the sources so that their currents into the nodes cancel, as the solve requires of them."))
end

# The targets of an adjoint judged by the floating subnetworks they feed:
# the derivative along a target that drives a net current into one is
# set by the gauge row, and only a combination of the targets whose net
# currents cancel in every subnetwork means anything, as the sum along two
# sources of one waveform into a subnetwork and out of it. A target is a
# unit current between its two terminals, so the targets are the edges of
# a multigraph on the floating subnetworks and the grounded rest of the
# circuit (`targetedges`), and the balanced combinations are the
# circulations on it. A target takes part in one exactly when its edge
# lies on a cycle. A bridge, whose removal separates its ends, alone
# carries the net current into the side it cuts off, which no other target
# returns, and is refused; a target whose terminals lie in one vertex
# balances alone.
function checkadjointtargets(p::TransientProblem, injection::SparseMatrixCSC, targets)
    islands = p.floatingcomponents
    isempty(islands) && return nothing
    beyond = multigraphbridges(targetedges(p, injection), 1 + length(islands))
    k = findfirst(>(0), beyond)
    isnothing(k) && return nothing
    nodes = join(p.circuit.nodenames[islands[beyond[k] - 1]], ", ")
    target = collect(targets)[k]
    throw(ArgumentError(lazy"the target $(target) drives a net current into the nodes ($(nodes)), which no element connects to ground, and no combination of the targets returns it, so its derivative means nothing: request it with the targets whose currents return its own and combine their derivatives so that the currents cancel, or connect the nodes to ground (a resistor or a capacitor will do)."))
end

# The bridges of the multigraph on the vertices `1:nv` whose edge `k`
# joins `ends[1, k]` and `ends[2, k]`: for each edge, the vertex its
# removal cuts off from the root of a depth-first search, or zero for an
# edge on a cycle or from a vertex to itself. Tarjan's search, without
# recursion: the edge to a vertex of the search tree is a bridge exactly
# when no edge from the vertex's subtree reaches above the vertex. The
# search steps back over the edge it came by alone, so two edges between
# the same vertices are each the other's way back. Each component is
# searched from its lowest vertex, so a bridge cuts off the side without
# vertex 1.
function multigraphbridges(ends::Matrix{Int}, nv::Int)
    # the edges at each vertex `u`, entries `start[u]:start[u + 1] - 1` of
    # `neighbor` and `edge`
    start = zeros(Int, nv + 1)
    start[1] = 1
    for k in axes(ends, 2)
        u, v = ends[1, k], ends[2, k]
        u == v && continue
        start[u + 1] += 1
        start[v + 1] += 1
    end
    cumsum!(start, start)
    neighbor, edge = zeros(Int, start[end] - 1), zeros(Int, start[end] - 1)
    cursor = start[1:nv]
    for k in axes(ends, 2)
        u, v = ends[1, k], ends[2, k]
        u == v && continue
        neighbor[cursor[u]], edge[cursor[u]] = v, k
        neighbor[cursor[v]], edge[cursor[v]] = u, k
        cursor[u] += 1
        cursor[v] += 1
    end
    # the search's path, each vertex with the edge it was reached by and
    # its next entry, and each vertex's discovery time and the earliest
    # its subtree reaches
    discovery, low = zeros(Int, nv), zeros(Int, nv)
    beyond = zeros(Int, size(ends, 2))
    path = Tuple{Int, Int, Int}[]
    clock = 0
    for root in 1:nv
        discovery[root] > 0 && continue
        clock += 1
        discovery[root] = low[root] = clock
        push!(path, (root, 0, start[root]))
        while !isempty(path)
            u, by, i = path[end]
            if i < start[u + 1]
                path[end] = (u, by, i + 1)
                w, k = neighbor[i], edge[i]
                k == by && continue
                if discovery[w] == 0
                    clock += 1
                    discovery[w] = low[w] = clock
                    push!(path, (w, k, start[w]))
                else
                    low[u] = min(low[u], discovery[w])
                end
            else
                pop!(path)
                isempty(path) && continue
                parent = first(path[end])
                low[parent] = min(low[parent], low[u])
                low[u] > discovery[parent] && (beyond[by] = u)
            end
        end
    end
    return beyond
end

# the problem a plan or a bath is built against, from whatever names it
transientproblemof(p::TransientProblem) = p
transientproblemof(s::TransientSolution) = s.problem
transientproblemof(b::TransientBatchSolution) = first(b.problems)

# the row of a port's trace from its number, which is what the user
# facing functions take; the traces of a solution are in the order of
# the ports' numbers
portindex(p::TransientProblem, number) = portindex(p.ports, number)
portrows(p::TransientProblem, ports) = Int[portindex(p, number) for number in ports]

# which ports the targets are, zero for a component, for the direct
# feedthrough of a target's current into a port wave
targetports(p::TransientProblem, targets) = [t isa Integer ? portindex(p, t) : 0 for t in targets]

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
- `ports`: the row of the port whose termination each component is, the
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

# The perturbation of the named components of a problem, on the backend:
# the stored entries of every component's derivative, and, with
# `forcing`, the stacked derivatives a tangent's forcing is one product
# of; an adjoint asks for the entries alone, so nothing it holds grows
# with the state times the components. The entries are read here, once,
# from each component's stamp, which no step builds.
function componentperturbation(p::TransientProblem, names, backend; forcing::Bool)
    psc, nm = p.circuit, p.matrices
    n, Lscale = length(p), p.Lscale
    Ljb = nm.Ljb
    nj = length(Ljb.nzval)
    RJ, lmolj = p.RJ, p.lmolj
    RJt = sparse(transpose(RJ))
    lookups = componentlookups(p.coupledbranches, Ljb)
    isempty(names) && throw(ArgumentError("name at least one component."))
    nc = length(names)
    tC, tG, tL, tJ = PerturbationTriplets(n), PerturbationTriplets(n), PerturbationTriplets(n), PerturbationTriplets(nj)
    ports = zeros(Int, nc)
    for (c, name) in enumerate(names)
        idx = componentindex(psc, name)
        s = componentstamp(idx, psc, nm, lookups, 1)
        # the capacitance enters as `r C`, the others inversely, so their
        # derivatives carry the minus of `d(1/(r p))/dr`; the entries are
        # on the node rows, which lead the state, a junction's those of
        # its column of the incidence, read from that column alone
        if s.kind == :Lj
            info = s.junction
            for q in nzrange(RJt, info)
                v = -lmolj[info]*nonzeros(RJt)[q]
                iszero(v) && continue
                push!(tJ.i, rowvals(RJt)[q]); push!(tJ.j, info); push!(tJ.v, v); push!(tJ.c, c)
            end
        else
            t, sc = s.kind == :C ? (tC, Lscale) : s.kind == :G ? (tG, -Lscale) : (tL, -Lscale)
            for q in eachindex(s.vals)
                v = sc*real(s.vals[q])
                iszero(v) && continue
                push!(t.i, s.rows[q]); push!(t.j, s.cols[q]); push!(t.v, v); push!(t.c, c)
            end
        end
        # a port's own termination: the environments are listed by port
        # number, as the traces are
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

# whether two perturbations built on one system are the same: a response
# builds its perturbation anew from the names it is given, so a kept
# workspace is matched on the names, in order, which resolve to the same
# components and the same entries on the system the workspace holds
sameperturbation(a, b) = isnothing(a) && isnothing(b)
sameperturbation(a::ComponentPerturbation, b::ComponentPerturbation) = a.names == b.names

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

# `X = F \ B` on every column at once, through the package's solve
# (`trysolve!`), which falls back to `\` for a method without an in place
# one (QR); cuDSS's in its extension
matrixsolve!(X, factor, B) = trysolve!(X, factor, B)

# The readings of a trapezoidal or backward Euler record the tangent
# and the adjoint share: the window of its states, which the
# perturbation's endpoint reads, and `readingat!(o, k)`, the projected
# junctions' phases and read rates at time `k` from the record into the
# reading work `o`, and with `withstates` the state where the
# perturbation reads it, into `xk`, with its rate as the solve read it.
function recordreadings(sol::TransientSolution, sys::TransientSystem, withstates::Bool, xk)
    statewindow = ResponseWindow(1:length(sol.times), nothing, nothing, nothing, nothing,
        (x, v, k) -> (copyto!(x, view(sol.flux, :, k)); copyto!(v, view(sol.rate, :, k)); nothing), nothing)
    readingat! = (o, k) -> begin
        prj = sys.projection
        if !isnothing(prj) && !isempty(prj.pj)
            copyto!(o.hphi, Array(view(sol.phases, prj.pj, k)))
            copyto!(o.er, Array(view(sol.endrates, :, k)))
        end
        if withstates
            copyto!(xk, view(sol.flux, :, k))
            copyto!(o.wk, view(sol.rate, :, k))
        end
        nothing
    end
    return statewindow, readingat!
end

"""
    transienttangent(solution, currents; targets = the ports,
        initialstate = nothing, factorization = nothing, reuse = nothing,
        outputsink = nothing, statesink = nothing)

The tangent of a recorded transient along a perturbation: `currents[q, k]`
is an additional Norton current in Amperes at target `q` (a port number,
or a component name, see [`transientinjection`](@ref)) at recorded time
`k`, and `initialstate` an optional tuple of perturbations at the start:
of the scaled flux and rate, `(flux, rate)`, each `(state, direction)`,
and under [`GaussLegendre`](@ref) optionally of the waves the lines
carried before it, `(line port, prehistory column, direction)`, and of
the rational blocks' states, `(state, direction)`, as its third and
fourth members, each with a trailing dimension of the conditions of a
batch when it differs between them. A third dimension of `currents` is a set of
directions propagated together, each step's factorization serving them
all. The currents into a subnetwork no element connects to ground, which
current sources alone feed, must cancel at every recorded time and at
every stage the rule reads, as the solve requires of its sources: a
tangent along a net current into one is refused. Under
[`GaussLegendre`](@ref) a current on the grid is read at the
stage times through a cubic Lagrange stencil, through the line between
the step's grid values on a record shorter than four points; a current
given as `currents[q, s, k, direction]`, with `s = 1` its value at
recorded time `k` and `s = 2, 3` its values at the two stage times of
the step from `k` to `k + 1`, is read as it is, so a pulse keeps its
support and a tone its exact phase at the stages; the trapezoidal rule
reads the grid values of either form. The rate of a current on the
grid, which the reading of the rate along an algebraic direction
carries, is that of the cubic through four grid values, of the
quadratic through three or of the line through two. Returns `(voltage,
incident, outgoing, finalflux, finalrate, finalwaves, finalstates)` in
the units of the solve, on its backend, with the directions as the
trailing dimension, under [`GaussLegendre`](@ref) `finalwaves` the waves
leaving each line port over the delay window before the end, `(line
port, column, direction)` as the third member of `initialstate` takes
them, so that a tangent of the next record continues this one, and
`finalstates` the rational blocks' states at the end, each `nothing`
without lines or blocks or under another rule; with an `outputsink`, a
function `outputsink(k, voltage, incident, outgoing)` receiving the
three port by direction matrices of recorded time `k` on the backend,
valid until the next call, the histories are not stored and those three
are `nothing`, so a measurement of a long record needs no memory per
time. A `statesink(k, flux, rate, states, waves)` receives the
tangent's state at each recorded time `k` in the same way, its scaled
fluxes and rates, its blocks' states as `finalstates` holds them and
under [`GaussLegendre`](@ref) the waves leaving its line ports at the
time as `finalwaves` holds them, state or port by direction, each
`nothing` without blocks or lines.
The linearization is about the full recorded state, so the loaded
junction phases enter every response; both the trajectory and its grid
are held fixed.
"""
Base.@nospecializeinfer function transienttangent(sol::TransientSolution,
        @nospecialize(currents::Union{Nothing,AbstractArray{<:Real}});
        targets = porttargets(sol.problem), initialstate = nothing,
        factorization = nothing, reuse = nothing, outputsink = nothing, statesink = nothing, perturbation = nothing)
    # compiled once whatever the arguments are: the responses run on a
    # batch, of which a solution is one condition, and its entry brings
    # them to the forms the steps take
    @nospecialize targets initialstate factorization reuse outputsink statesink perturbation
    return map(dropcondition, transienttangent(batchof(sol), currents;
        targets, initialstate, factorization, reuse, outputsink, statesink, perturbation))
end

# the tangent of a trapezoidal or backward Euler solve, on the forms the
# batch's entry made of its arguments
function steptangent(sol::TransientSolution, currents::TangentCurrents, injh::SparseMatrixCSC{Float64, Int},
        tp::Vector{Int}, initial::TangentInitial, sys::TransientSystem, factorization::AbstractFactorization, kept,
        @nospecialize(outputsink), @nospecialize(statesink), perturbation, @nospecialize(reuse))
    p = sol.problem
    backend = sys.backend
    n, np, nt = length(p), length(p.portimpedances), length(sol.times)
    nq, ndir = size(injh, 2), directions(currents)
    h = sys.h
    trapezoidal = sys.method isa Trapezoidal
    # the components' forcing of each step, from the recorded states, and
    # the direct term of the port waves along a port's own termination
    pwork = isnothing(perturbation) ? nothing : perturbationwork(perturbation, sys, 1; forcing = true)
    dterm = isnothing(perturbation) ? nothing : directterm(perturbation, sys, 1, ndir)
    injection = devicesparse(injh, backend)
    # the trapezoidal step reads the grid values of a staged current; a
    # tangent along the components alone carries one zero column of them;
    # where each direction drives one target, the currents of every target
    # at a time are its waveform there through the selector of the targets
    dI = tobackend(backend, currents.grid)
    kcol = k -> gridcolumn(currents, k)
    allocate = (dims...) -> KernelAbstractions.zeros(backend, Float64, dims...)
    selector = isempty(currents.targets) ? nothing : tobackend(backend, targetselector(currents.targets, nq))
    targetwork = allocate(isnothing(selector) ? 0 : nq, ndir)
    gridcurrent = k -> isnothing(selector) ? view(dI, :, kcol(k), :) :
        settargetcurrents!(targetwork, selector, view(dI, :, kcol(k), :))
    dx, dv, dxnew, rhs, work = allocate(n, ndir), allocate(n, ndir), allocate(n, ndir), allocate(n, ndir), allocate(n, ndir)
    if initial.given
        copyto!(dx, conditionslice(initial.flux, 1))
        copyto!(dv, conditionslice(initial.rate, 1))
    end
    nj = length(sys.lmolj)
    phi, jwork = allocate(nj), allocate(nj, ndir)
    # the port outputs of every time, stored, or of none where a sink
    # receives them, empty then, as in the Gauss-Legendre tangent
    voltage, incident, outgoing = [allocate(np, isnothing(outputsink) ? nt : 0, ndir) for _ in 1:3]
    coefficients = Dict(q => outputcoefficients(sys, q) for q in (:voltage, :incident, :outgoing))
    portwork = allocate(np, ndir)
    # the current of the targets that are ports, into their port waves
    portmap = devicesparse(targetportmap(tp, np), backend)
    directwork = allocate(np, ndir)
    outwork = [allocate(np, ndir) for _ in 1:3]
    # where a port or the state sink reads a rate along an algebraic
    # direction, the port waves and the sink's rate come from the tangent
    # rate read as the solve's is, the reading linearized at the recorded
    # phases and read rates of the projected junctions, with the tangent
    # currents' rate and the components' perturbation of the constraints;
    # the work is there whether or not anything reads, so that the
    # closure of the outputs is one type
    reading = outputreading(sys, backend, n, 1, ndir)
    reads = sys.portsread || !isnothing(statesink)
    withstates = !isnothing(perturbation) && perturbation.states
    xk = allocate(n, 1)
    # the projection of the endpoint linearized at the recorded endpoint,
    # with the constraints perturbed there by the components and by the
    # tangent currents: `dx -= Z (Z' (L + J'(x)) Z)^-1 Z' ((L + J'(x)) dx + f)`
    pr = sys.projection
    projecting = !isnothing(pr) && !isempty(pr.directions)
    cw = projecting ? projectionwork(pr, backend, n, 1, ndir) : nothing
    pend = (projecting && !isnothing(perturbation)) ? endpointwork(perturbation, sys, p, 1, ndir; forcing = true) : nothing
    fdev, fh = projecting ? (allocate(n, ndir), allocate(n, ndir)) : (nothing, nothing)
    statewindow, readingat! = recordreadings(sol, sys, withstates, xk)
    function outputs!(k)
        if reads
            readingat!(reading, k)
            readoutputs!(reading, sys, k, dv, dx, currents, injh, withstates ? Array(xk) : zeros(0, 0),
                withstates ? Array(reading.wk) : zeros(0, 0), perturbation, backend)
        end
        stepmul!(portwork, sys.ports, reads ? reading.dvread : dv)
        portwork .*= phi0
        stepmul!(directwork, portmap, gridcurrent(k))
        outs = isnothing(outputsink) ? (view(voltage, :, k, :), view(incident, :, k, :), view(outgoing, :, k, :)) : outwork
        for (s, q) in enumerate((:voltage, :incident, :outgoing))
            cv, cd = coefficients[q]
            outs[s] .= cv .* portwork .+ cd .* directwork
        end
        isnothing(dterm) || adddirectterm!(outs, dterm, sol, [p], k)
        isnothing(outputsink) || outputsink(k, outwork[1], outwork[2], outwork[3])
        isnothing(statesink) || statesink(k, dx, reading.dvread, nothing, nothing)
        nothing
    end
    outputs!(1)
    factor = kept
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
            stepmul!(work, injection, gridcurrent(k - 1)); rhs .+= work
        end
        stepmul!(work, injection, gridcurrent(k)); rhs .+= work
        projecting && (fdev .= .-work)
        if !isnothing(perturbation)
            stepquantities!(pwork, perturbation, sys, sol, k)
            perturbationforcing!(pwork.F, perturbation, pwork, pwork.work)
            rhs .+= reshape(pwork.F, n, ndir)
        end
        # the step matrix at the recorded phases, and the solves; without
        # a junction it does not change, and is factorized once
        if nj > 0 || isnothing(factor)
            copyto!(phi, view(sol.phases, :, k))
            factor = stepjacobian!(sys, factorization, phi, factor)
        end
        matrixsolve!(dxnew, factor, rhs)
        if projecting
            if !isnothing(pend)
                endpointquantities!(pend, perturbation, statewindow, Array(view(sol.phases, pr.pj, k)), pr.pj, k)
                endpointforcing!(pend, perturbation)
                copyto!(fh, pend.F)
                fdev .-= fh
            end
            copyto!(cw.hphi, Array(view(sol.phases, pr.pj, k)))
            constrainttangent!(cw, pr, dxnew, fdev)
            constraintcorrection!(cw, pr)
            dxnew .+= cw.work
        end
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
    # the tangent's rate along the algebraic directions, as the
    # Gauss-Legendre tangent reads it
    finalrate = copy(dv)
    (isnothing(sys.invariant) && isnothing(sys.projection)) ||
        readtangentrate!(finalrate, dv, dx, sys, sol, 1, ndir, perturbation, currents, injh, backend)
    squeeze = a -> isnothing(a) ? nothing : currents.single ? reshape(a, size(a)[1:end-1]...) : a
    stored = isnothing(outputsink)
    return (; voltage = stored ? squeeze(voltage) : nothing, incident = stored ? squeeze(incident) : nothing,
        outgoing = stored ? squeeze(outgoing) : nothing, finalflux = squeeze(copy(dx)), finalrate = squeeze(finalrate),
        finalwaves = nothing, finalstates = nothing)
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

# the pieces of the adjoint of a trapezoidal or backward Euler record, of
# one condition: the weights on the backend, the output coefficients, the
# transposed injection, and the ring with the feedthrough, its sink the
# storing one where none is given (see `adjointmaps`)
function adjointsetup(wh::Array{Float64, 3}, quantity::Symbol, injh::SparseMatrixCSC{Float64, Int}, tp::Vector{Int},
        sys::TransientSystem, @nospecialize(sink))
    backend = sys.backend
    np, nt, nobj = size(wh)
    nq = size(injh, 2)
    w = tobackend(backend, wh)
    cv, cd = outputcoefficients(sys, quantity)
    injectiont, portmapt = adjointmaps(injh, tp, np, backend)
    allocate = (dims...) -> KernelAbstractions.zeros(backend, Float64, dims...)
    currents = isnothing(sink) ? allocate(nq, nt, nobj) : nothing
    ring = CurrentRing(backend, nq, nobj, isnothing(sink) ? storingsink(currents) : sink,
        currentfeedthrough(w, cd, portmapt, allocate(np, nobj), allocate(nq, nobj), 1))
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
# currents built in the steps' form, as the noise builds its baths'
# quadratures: one column per recorded time, and each direction's target,
# where it has one, among the targets
function currentshape(c::TangentCurrents, nq, nt)
    rows = isempty(c.targets) ? nq : 1
    (size(c.grid, 1) == rows && size(c.grid, 2) == nt && (isempty(c.targets) || length(c.targets) == directions(c)) &&
        all(q -> 1 <= q <= nq, c.targets)) || throw(DimensionMismatch(
        lazy"the currents need one column per recorded time ($(nt)) and their targets among the $(nq) given."))
    (all(isfinite, c.grid) && all(isfinite, c.stages)) || throw(ArgumentError("the current perturbations must be finite."))
    return directions(c)
end

"""
    transientadjoint(solution, weights; quantity = :outgoing, targets = the ports,
        factorization = nothing, reuse = nothing, sink = nothing, stagesink = nothing,
        components = String[])

The exact discrete adjoint of `sum(weights .* getproperty(solution,
quantity))` on the recorded grid, `weights` having one row per port and
one column per recorded time and `quantity` being `:voltage`, `:incident`
or `:outgoing`; a third dimension of `weights` is a set of objectives
propagated together on each step's factorization. Returns `(currents,
initialflux, initialrate, initialwaves, initialstates, sensitivity)`:
the derivatives with respect to a Norton current at every target (a
port number or a component name, see [`transientinjection`](@ref)) and
recorded time, in Amperes, including the direct feedthrough of a port's
current into its wave; with respect to the scaled initial flux and
rate; under [`GaussLegendre`](@ref) with respect to the lines'
prehistory and the rational blocks' initial states, the cotangents of
the third and fourth members of the tangent's `initialstate`, or
`nothing` without them; and with respect to the relative values of the
named `components`, `(component, objective)`, as
[`transientsensitivity`](@ref) perturbs them, or `nothing` without any;
all with the objectives as the trailing dimension. For a consistent
perturbation the contraction equals
the weighted output of [`transienttangent`](@ref); a complex demodulation
is two real objectives, and a time integral carries its quadrature
weights in `weights`. The derivative along a target that drives a net
current into a subnetwork no element connects to ground, a current
source between it and the rest of the circuit, depends on which node of
the subnetwork the solver takes as its flux reference: only a combination
of the targets whose currents cancel there means anything, as the sum
along two sources of one waveform into it and out of it, and a target
that no such combination includes is refused. With a `sink`, a function
`sink(k, values)`, the currents are not stored: each column, a targets
by objectives matrix on the backend valid until the next call, is handed
to the sink once it is final, in decreasing recorded time, and
`currents` is `nothing`; a contraction over a long record then needs no
memory per time. With a
`stagesink` as well, a function `stagesink(k, i, values)`, the
multipliers of the two stages of each Gauss-Legendre step from `k` to
`k + 1` go to it as they are, at their stage times, in a buffer which,
like the sink's columns, is valid until the next call, and the grid columns
hold only what the grid reads; that pair is the transpose of the staged
form of the tangent's currents, and it is what the noise contracts.
"""
Base.@nospecializeinfer function transientadjoint(sol::TransientSolution, @nospecialize(weights::AbstractArray{<:Real});
        quantity::Symbol = :outgoing, targets = porttargets(sol.problem),
        factorization = nothing, reuse = nothing, sink = nothing, stagesink = nothing,
        components = String[])
    # compiled once whatever the arguments are, on the batch of one
    # condition the solution is (see transienttangent)
    @nospecialize targets factorization reuse sink stagesink components
    return map(dropcondition, transientadjoint(batchof(sol), weights;
        quantity, targets, factorization, reuse, sink, stagesink, components))
end

# the adjoint of a trapezoidal or backward Euler solve, on the forms the
# batch's entry made of its arguments
function stepadjoint(sol::TransientSolution, wh::Array{Float64, 3}, single::Bool, quantity::Symbol,
        injh::SparseMatrixCSC{Float64, Int}, tp::Vector{Int}, sys::TransientSystem, factorization::AbstractFactorization, kept,
        @nospecialize(sink), perturbation, @nospecialize(reuse))
    p = sol.problem
    backend = sys.backend
    n, np, nt = length(p), length(p.portimpedances), length(sol.times)
    nq, nobj = size(injh, 2), size(wh, 3)
    h = sys.h
    trapezoidal = sys.method isa Trapezoidal
    w, cv, injectiont, ring, currents = adjointsetup(wh, quantity, injh, tp, sys, sink)
    allocate = (dims...) -> KernelAbstractions.zeros(backend, Float64, dims...)
    # the components' forcing of each step, contracted against the step's
    # multipliers on the entries of the components, and the direct term of
    # the objective along a port's own termination
    pwork = isnothing(perturbation) ? nothing : perturbationwork(perturbation, sys, 1; forcing = false)
    pcwork = isnothing(perturbation) ? nothing : contractionwork(perturbation, perturbation.entries, 1, nobj, backend)
    sensitivity = isnothing(perturbation) ? nothing : allocate(length(perturbation.names), nobj, 1)
    isnothing(perturbation) || (sensitivity .+= tobackend(backend,
        adddirectsensitivity!(zeros(length(perturbation.names), nobj, 1), perturbation, sys, sol, [p], wh, quantity)))
    xbar, vbar, lambda, work = allocate(n, nobj), allocate(n, nobj), allocate(n, nobj), allocate(n, nobj)
    nj = length(sys.lmolj)
    phi, jwork = allocate(nj), allocate(nj, nobj)
    portwork = allocate(np, nobj)
    targetwork = allocate(nq, nobj)
    # the adjoint of an output at time `k`, as the Gauss-Legendre
    # adjoint's (see `adjointoutput!`), its work there whether or not a
    # port reads, so that the closure is one type
    transposing = outputtranspose(sys, backend, n, 1, nobj)
    withstates = !isnothing(perturbation) && perturbation.states
    xk = allocate(n, 1)
    # the projection of the endpoint transposed (see `constrainttranspose!`):
    # the cotangent of the flux before it, and through the cotangent of
    # the projection's forcing the current at the endpoint time and the
    # components, contracted on the host
    pr = sys.projection
    projecting = !isnothing(pr) && !isempty(pr.directions)
    cw = projecting ? projectionwork(pr, backend, n, 1, nobj) : nothing
    fbar = projecting ? allocate(n, nobj) : nothing
    pend = (projecting && !isnothing(perturbation)) ? endpointwork(perturbation, sys, p, 1, nobj; forcing = false) : nothing
    sensh = isnothing(pend) ? nothing : zeros(length(perturbation.names), nobj, 1)
    statewindow, readingat! = recordreadings(sol, sys, withstates, xk)
    function output!(k)
        portwork .= cv .* view(w, :, k, :)
        stepmul!(work, sys.portst, portwork)
        if !sys.portsread
            vbar .+= phi0 .* work
            return nothing
        end
        o = transposing
        o.rbar .= phi0 .* work
        readingat!(o, k)
        readoutputstranspose!(o, sys, k, nt, vbar, xbar, ring, injectiont, targetwork, perturbation, pcwork, sensitivity,
            xk, backend)
        nothing
    end
    output!(nt)
    factor = kept
    for k in nt:-1:2
        # the rate update feeds the flux adjoint, then the step's solve is
        # transposed on the symmetric step matrix at the recorded phases
        xbar .+= (trapezoidal ? 2/h : 1/h) .* vbar
        if projecting
            copyto!(cw.hphi, Array(view(sol.phases, pr.pj, k)))
            constrainttranspose!(xbar, fbar, cw, pr, sys, backend)
            stepmul!(targetwork, injectiont, fbar)
            ringadd!(ring, k, targetwork, -1.0)
            if !isnothing(pend)
                endpointquantities!(pend, perturbation, statewindow, Array(view(sol.phases, pr.pj, k)), pr.pj, k)
                endpointcontract!(sensh, pend, perturbation, Array(fbar), -1.0)
            end
        end
        if nj > 0 || isnothing(factor)
            copyto!(phi, view(sol.phases, :, k))
            factor = stepjacobian!(sys, factorization, phi, factor)
        end
        matrixsolve!(lambda, factor, xbar)
        # the components' forcing of the step against its multipliers
        if !isnothing(perturbation)
            stepquantities!(pwork, perturbation, sys, sol, k)
            contractperturbation!(sensitivity, perturbation, perturbation.entries, pcwork, pwork, lambda, 1.0)
        end
        # the current perturbations the step read; the column at k is
        # final once this step has added to it, except for the first four,
        # which the reading of a rate along an algebraic direction at the
        # first three times still reaches
        stepmul!(targetwork, injectiont, lambda)
        ringadd!(ring, k, targetwork, 1.0)
        trapezoidal && ringadd!(ring, k - 1, targetwork, 1.0)
        k > 4 && ringemit!(ring, k)
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
    for j in min(4, nt):-1:1
        ringemit!(ring, j)
    end
    KernelAbstractions.synchronize(backend)
    isnothing(reuse) || (reuse.factor = factor)
    isnothing(sensh) || (sensitivity .+= tobackend(backend, sensh))
    squeeze = a -> isnothing(a) ? nothing : single ? reshape(a, size(a)[1:end-1]...) : a
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
finalflux, finalrate, finalwaves, finalstates)` as [`transienttangent`](@ref) does, with the
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
Base.@nospecializeinfer function transientsensitivity(sol::Union{TransientSolution,TransientBatchSolution}, @nospecialize(names);
        factorization = nothing, reuse = nothing, outputsink = nothing)
    # compiled once whatever the arguments are (see transienttangent)
    @nospecialize factorization reuse outputsink
    recordedsolution(sol)
    p = transientproblemof(sol)
    backend = KernelAbstractions.get_backend(sol.finalflux)
    cp = componentperturbation(p, names, backend; forcing = true)
    recordedstates(sol, cp)
    return transienttangent(sol, nothing; factorization, reuse, outputsink, perturbation = cp)
end
