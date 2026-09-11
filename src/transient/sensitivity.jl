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

# which ports the targets are, zero for a component, for the direct
# feedthrough of a target's current into a port wave
targetports(p::TransientProblem, targets) = [t isa Integer ? findfirst(port -> port.number == t, p.circuit.ports) : 0 for t in targets]

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
function transienttangent(sol::TransientSolution, currents::AbstractArray{<:Real};
        targets = porttargets(sol.problem), initialstate = nothing,
        factorization = nothing, reuse = nothing, outputsink = nothing)
    recordedsolution(sol)
    p = sol.problem
    backend = KernelAbstractions.get_backend(sol.finalflux)
    fact = isnothing(factorization) ? transientfactorization(backend) : factorization
    sys = transientsystem(reuse, p, sol.dt, sol.method, backend, fact)
    n, np, nt = length(p), length(p.portimpedances), length(sol.times)
    nq = length(targets)
    ndir = currentshape(currents, nq, nt)
    sys.method isa GaussLegendre && return gaussbatchtangent(sol, currents, targets, initialstate, sys, outputsink)
    h = sys.h
    trapezoidal = sys.method isa Trapezoidal
    injection = devicesparse((sys.Lscale/phi0) .* transientinjection(p, targets), backend)
    ports = tobackend(backend, [q > 0 ? q : 0 for q in targetports(p, targets)])
    # the trapezoidal step reads the grid values of a staged current
    dI = tobackend(backend, reshape(Float64.(collect(ndims(currents) == 4 ? selectdim(currents, 2, 1) : currents)), nq, nt, ndir))
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
        stepmul!(directwork, portmap, view(dI, :, k, :))
        for (s, q) in enumerate((:voltage, :incident, :outgoing))
            cv, cd = coefficients[q]
            out = isnothing(outputsink) ? view((voltage, incident, outgoing)[s], :, k, :) : outwork[s]
            out .= cv .* portwork .+ cd .* directwork
        end
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
            stepmul!(work, injection, view(dI, :, k - 1, :)); rhs .+= work
        end
        stepmul!(work, injection, view(dI, :, k, :)); rhs .+= work
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
    squeeze = a -> isnothing(a) ? nothing : ndims(currents) == 2 ? reshape(a, size(a)[1:end-1]...) : a
    return (; voltage = squeeze(voltage), incident = squeeze(incident), outgoing = squeeze(outgoing),
        finalflux = squeeze(copy(dx)), finalrate = squeeze(copy(dv)))
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
        factorization = nothing, reuse = nothing, sink = nothing, stagesink = nothing)
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
    sys.method isa GaussLegendre && return gaussbatchadjoint(sol, weights, quantity, targets, sys, sink, stagesink)
    h = sys.h
    trapezoidal = sys.method isa Trapezoidal
    w, cv, injectiont, ring, currents = adjointsetup(sol, weights, quantity, targets, sys, sink)
    allocate = (dims...) -> KernelAbstractions.zeros(backend, Float64, dims...)
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
    return (; currents = squeeze(currents), initialflux = squeeze(xbar), initialrate = squeeze(vbar), initialwaves = nothing, initialstates = nothing)
end
