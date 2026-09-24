# `method = Staged()`: source continuation on an adaptively grown harmonic
# grid, for operating points the direct methods cannot reach from a cold
# start.

"""
    StagedStageInfo

One attempted stage of [`stagedhbnlsolve`](@ref). Every attempt is
stored in `solverinfo.stages` of the returned solution in the order it
ran, including stalled steps and growth retreats, so the whole
continuation walk can be examined afterwards.

# Fields
- `label`: `"staged"`.
- `converged`: whether this attempt's inner solve converged.
- `iterations`: total inner Newton iterations of the attempt.
- `grid`: the harmonic truncation the attempt solved on.
- `sfrom`: the last accepted drive fraction before the attempt.
- `starget`: the drive fraction the attempt targeted.
- `ds`: `starget - sfrom` (negative for a growth retreat).
- `action`: `:advance` (drive step on the current grid), `:grow` (a solve
    on a newly grown grid after carrying a converged point up, and the
    drive retreats of a carried point which has not converged on its grid
    yet), or `:final` (a full-drive solve on the finest grid).
- `accepted`: whether the attempt's result was kept as the new operating
    point (a stalled attempt is recorded but not accepted).
- `seconds`: wall time of the attempt, including the stage's system
    assembly.
- `finalresidual`: the residual norm the attempt ended at.
- `inner`: the inner solver's own stage records ([`IterationInfo`](@ref)),
    with their Krylov linear-solve diagnostics when the inner method is a
    [`NewtonKrylov`](@ref).
"""
struct StagedStageInfo <: AbstractStageInfo
    label::String
    converged::Bool
    iterations::Int
    grid::Tuple
    sfrom::Float64
    starget::Float64
    ds::Float64
    action::Symbol
    accepted::Bool
    seconds::Float64
    finalresidual::Float64
    inner::Vector
end

function Base.show(io::IO, ::MIME"text/plain", r::StagedStageInfo)
    print(io, "StagedStageInfo: grid=", r.grid, " s ", round(r.sfrom, digits = 4),
        " -> ", round(r.starget, digits = 4), " ", r.action,
        r.accepted ? " accepted" : " stalled",
        " newton=", r.iterations,
        " |F|=", round(r.finalresidual, sigdigits = 2),
        " (", round(r.seconds, digits = 2), " s)")
end

"""
    defaultgridladder(Nharmonics::NTuple{N,Int})

The default coarse to fine ladder of retained harmonic caps for
[`stagedhbnlsolve`](@ref): `Nharmonics` halved repeatedly in every
dimension down to two harmonics, a tone with two or fewer keeping its
own, coarsest first and `Nharmonics` last.
"""
function defaultgridladder(Nharmonics::NTuple{N,Int}) where {N}
    grids = [Nharmonics]
    g = Nharmonics
    while any(x -> x > 2, g)
        g = map(x -> x > 2 ? max(2, cld(x, 2)) : x, g)
        push!(grids, g)
    end
    return reverse(grids)
end

# The circuit matrices at the drive fraction `s`: the netlist's constant
# current sources scaled with the ports' sources, as `setdrive!` scales the
# whole drive, and `nm` itself when there are none.
function drivenmatrices(nm::CircuitMatrices, componenttypes::Vector{Symbol},
        s::Real)
    any(==(:I), componenttypes) || return nm
    vvn = Any[t === :I && v isa Number ? s*v : v
        for (t, v) in zip(componenttypes, nm.vvn)]
    return CircuitMatrices(nm.Cnm, nm.Gnm, nm.Lb, nm.Lbm, nm.Ljb, nm.Ljbm,
        nm.Mb, nm.invLnm, nm.Rbnm, nm.portindices, nm.portnumbers,
        nm.portimpedances, nm.portenvironmentindices,
        nm.noiseportimpedanceindices, nm.Lmean, vvn)
end

# Embed a converged solution as the initial guess on a larger grid by
# matching mode tuples; modes the small grid did not carry start at zero.
function stagedembed(out, bigmodes)
    small = reshape(Array(out.nodeflux), out.Nmodes, :)
    pos = Dict(m => i for (i, m) in enumerate(bigmodes))
    X = zeros(ComplexF64, length(bigmodes), size(small, 2))
    for (i, m) in enumerate(out.modes)
        haskey(pos, m) && (X[pos[m], :] = small[i, :])
    end
    return vec(X)
end

"""
    stagedhbnlsolve(m::Staged, w::NTuple{N,Float64}, Nharmonics,
        sources::Vector{SourceTuple{N}}, psc::CompiledCircuit,
        circuitdefs::Dict{Any,Any}; kwargs...)

Source continuation on an adaptively grown harmonic grid, reached through
`hbnlsolve(...; method = Staged(...))`; the schedule is the [`Staged`](@ref)
value `m`, validated at its construction. `kwargs` are the keywords of
[`hbnlsolve`](@ref) (`iterations`, `atol`, `Nevaluationharmonics`,
`frequencywindow`, `maxintermodorder`, `dc`, `odd`, `even`,
`keyedarrays`, `sensitivitynames`, `returnoperatingpoint`,
`backend`), which are forwarded to every stage.

Near a critical drive the Newton basin is small and the iteration count
large, so those iterations are spent where they are cheap. The drive is
climbed in warm started steps with only a small set of harmonics retained
as unknowns, and each larger retained set is warm started from the last by
matching mode tuples. A step scales every source of the drive by its
fraction, the ports' and the netlist's constant current sources alike, as
[`setdrive!`](@ref) scales a problem's. The nonlinearity is always evaluated on the full
`Nevaluationharmonics` transform grid, so that every stage sees the same
aliasing of the nonlinear products; the ladder only controls the modes
retained as unknowns, so it is the linear solves which shrink. The mode
set, the circuit matrices, the system with its transforms, the
preconditioner's symbolic factorization and the Krylov workspace of a
grid are built once and rebound to each drive step on it
([`HBReuse`](@ref)), so a step costs its Newton iterations and little
else.

The schedule adapts in both directions, because each truncation has its
own solvability boundary and the boundaries are not monotone in the grid:
a stalled drive step is halved; a stall at the minimum step grows the grid
at the current converged drive; and a point carried to a larger grid
which fails to reconverge there retreats the drive on that grid until it
converges, before the schedule grows past it. Interior points converge
only to `interioratol` under a small iteration budget, since they exist to
keep the iterate inside the basin, and the one expensive solve, the finest
grid at full drive, starts inside the basin with the caller's `atol` and
`iterations`.

Only a stall from a point converged on the finest grid itself is reported
as bracketing a fold, the end of the solution branch (the self
oscillation threshold) between the last converged drive fraction and the
stalled one; the report is what the search saw, not a proof that no
operating point exists at the requested drive. Before reporting it, and
only when the stalled target was below full drive, one further solve at
full drive is attempted from the last converged point with the caller's
own method and tolerance, since a coexisting branch may reach it; failing
that, the solve returns not converged with a warning stating the bracket,
and `solverinfo.sourcefold` holds the last converged drive fraction (it
stays `NaN` for the other ways a schedule ends). No path throws: a
schedule which cannot reach the point (its attempts spent, a carried
point which does not reconverge, a first stage which converges at no
drive down to the minimum step) warns, returns its last attempt marked
not converged, and records the whole walk in `solverinfo.stages`.

# Keywords
- `grids = defaultgridladder(Nharmonics)`: the coarse to fine ladder of
    retained harmonic caps, whose last entry must equal `Nharmonics`.
- `s0 = 0.5`: the first drive fraction attempted.
- `smin = 0.02`: the minimum drive step; a stall below it grows the grid.
- `interioratol = 1e-7`, `interioriterations = 60`: the tolerance and the
    Newton budget of the interior points. The budget is small on purpose: a
    stalled probe is evident within tens of iterations, and interior stalls
    are the overhead of the walk.
- `inner = NewtonKrylov()`: the method of every stage.
- `interiorescalation = false`: whether interior stage solves may escalate
    their preconditioner to the full Jacobian. Off by default, because an
    interior probe exists only to produce a cheap warm start, and at a high
    tone count the full factorization may not fit in memory; a probe which
    fails without escalation is a stall, which the schedule answers with a
    smaller step. The final solve keeps the caller's escalation behavior.
- `maxattempts = 60`: a bound on the total number of stage solves.
- `verbose = false`: print one line per stage solve.
- `warnnotconverged = true`: warn when the schedule ends without the
    requested point, and run the checks of the point it reaches which
    warn, once. The stage solves never warn: a stage which does not
    converge is how the schedule finds its step.
"""
function stagedhbnlsolve(m::Staged, w::NTuple{N,Float64},
    Nharmonics::NTuple{N,Int}, sources::Vector{SourceTuple{N}},
    psc::CompiledCircuit, circuitdefs::Dict{Any,Any};
    kwargs...) where {N}
    # a barrier on the inner method: the schedule holds it as an abstract
    # field, and the stage solves below must see its concrete type, or
    # every stage is called with keywords of unknown type. Through
    # `invokelatest` rather than a plain call, because a plain call with an
    # abstract argument is also inferred for that abstract signature (one
    # method matches), and that inference of the whole continuation and the
    # compiled circuit solve behind it, for an instance which never runs,
    # took longer than compiling the one which does.
    return Base.invokelatest(stagedhbnlsolve, m.inner, m, w, Nharmonics,
        sources, psc, circuitdefs; kwargs...)
end

function stagedhbnlsolve(inner::AbstractHBNonlinearSolver, m::Staged,
    w::NTuple{N,Float64}, Nharmonics::NTuple{N,Int},
    sources::Vector{SourceTuple{N}}, psc::CompiledCircuit,
    circuitdefs::Dict{Any,Any};
    iterations = 1000,
    Nevaluationharmonics::NTuple{N,Int} = map(i -> 2i, Nharmonics),
    frequencywindow = (0, Inf),
    maxintermodorder = Inf, dc::Bool = false, odd::Bool = true,
    even::Bool = false, atol = 1e-8,
    keyedarrays::Bool = true,
    sensitivitynames::AbstractVector = String[],
    returnoperatingpoint::Bool = false, backend = CPU(),
    warnnotconverged::Bool = true) where {N}

    (; s0, smin, interioratol, interioriterations, interiorescalation,
        maxattempts, verbose) = m
    # the ladder of retained harmonic caps: the default for this problem's
    # `Nharmonics`, or the one the schedule states, which must end at it
    grids = isnothing(m.grids) ? defaultgridladder(Nharmonics) :
        Vector{NTuple{N,Int}}(m.grids)
    isempty(grids) && throw(ArgumentError("`grids` must not be empty."))
    grids[end] == Nharmonics || throw(ArgumentError(
        lazy"the finest grid $(grids[end]) must equal `Nharmonics` = $(Nharmonics)."))

    # Every stage uses the full transform grid `Nevaluationharmonics`, so
    # the aliasing of the nonlinear products is the same on every stage;
    # the ladder only sets the modes retained as unknowns.
    #
    # What the stages of one grid share is built once per grid: its mode
    # set with the Fourier indices, the circuit matrices at its mode count,
    # and a reuse object holding the system, the preconditioner's symbolic
    # factorization and the Krylov workspace, which each drive step rebinds
    # to its sources rather than rebuilds (see `HBReuse`). The component
    # values are resolved once for every grid. The first grid's mode set
    # is built before any solve, which checks `Nevaluationharmonics`.
    vvn = componentvaluestonumber(psc.componentvalues, circuitdefs)
    function gridstate(grid)
        freq, indices = pumpmodeset(w, Nharmonics, Nevaluationharmonics;
            dc = dc, odd = odd, even = even,
            maxintermodorder = maxintermodorder,
            frequencywindow = frequencywindow,
            retained = map(min, Nharmonics, grid))
        return (freq = freq, indices = indices,
            nm = numericmatrices(psc, vvn; Nmodes = length(freq.modes)),
            reuse = HBReuse())
    end
    scaled(s) = SourceTuple{N}[(mode = t.mode, port = t.port,
        current = s*t.current) for t in sources]
    solve(state, s, x0, final) = hbnlsolve(w, scaled(s), state.freq,
        state.indices, psc, drivenmatrices(state.nm, psc.componenttypes, s);
        reuse = state.reuse,
        method = (final || interiorescalation) ? inner :
            withescalation(inner, false),
        # a stage solve which does not converge is how the continuation
        # finds its step rather than a failure of the solve the caller
        # asked for, and the outcome of the schedule is warned about and
        # checked here
        warnnotconverged = false,
        # typed here, whatever the loop below inferred for its carried point,
        # so the stage solve is called with keywords of known type
        x0 = initialguess(x0),
        keyedarrays = final ? keyedarrays : false,
        sensitivitynames = final ? sensitivitynames : String[],
        returnoperatingpoint = final ? returnoperatingpoint : false,
        backend = backend,
        atol = final ? atol : interioratol,
        iterations = final ? iterations : interioriterations)

    gi = 1
    state = gridstate(grids[gi])   # what the stages of the current grid share
    s = 0.0             # last converged drive fraction on the current grid
    ds = s0
    x = nothing         # its solution, in the raw nodeflux layout
    out = nothing
    attempts = 0
    stagerecords = AbstractStageInfo[]
    pendinggrow = false
    last = nothing      # the most recent attempt, converged or not
    gaveup = false      # the schedule ended without the requested point
    fold = NaN          # the last converged drive fraction below a fold
    record = function (cand, grid, sfrom, starget, action, accepted, secs)
        si = cand.solverinfo
        push!(stagerecords, StagedStageInfo("staged", si.converged,
            sum(st -> st.iterations, si.stages; init = 0), grid,
            Float64(sfrom), Float64(starget), Float64(starget - sfrom),
            action, accepted, secs, Float64(si.finalresidual), si.stages))
        return nothing
    end
    while true
        attempts += 1
        if attempts > maxattempts
            warnnotconverged && @warn lazy"the staged schedule did not converge in `maxattempts` = $(maxattempts) stage solves; last converged drive fraction $(s) on grid $(grids[gi])."
            gaveup = true
            break
        end
        starget = min(1.0, s + ds)
        final = gi == length(grids) && starget >= 1.0
        t0 = time_ns()
        cand = solve(state, starget, x, final)
        last = cand
        ok = cand.solverinfo.converged
        # the first solve after carrying a full drive point to a larger
        # grid is a growth, not a drive advance
        record(cand, grids[gi], s, starget,
            final ? :final : (pendinggrow ? :grow : :advance), ok,
            (time_ns() - t0)/1e9)
        ok && (pendinggrow = false)
        if verbose
            st = cand.solverinfo.stages[end]
            println("staged: grid=", grids[gi], " s=", round(starget, digits = 4),
                " newton=", st.iterations,
                " |F|=", round(cand.solverinfo.finalresidual, sigdigits = 2),
                ok ? "" : " STALL")
        end
        if ok
            out = cand
            s = starget
            final && break
            if s >= 1.0
                # full drive reached on a coarse grid: grow the grid
                state = gridstate(grids[gi+1])
                x = stagedembed(out, state.freq.modes)
                gi += 1
                pendinggrow = true
            else
                x = vec(reshape(Array(out.nodeflux), out.Nmodes, :))
                ds = min(2*ds, 1.0 - s)
            end
        elseif starget - s > smin
            # halve the effective step rather than `ds`, so that a target
            # capped at full drive is not re-attempted identically
            ds = (starget - s)/2
        elseif isnothing(x)
            # nothing has converged: the first grid stalls at every drive
            # down to the minimum step, and there is no point to grow from
            # or to bracket a fold with
            warnnotconverged && @warn lazy"the first stage converged at no drive fraction down to $(round(starget, sigdigits = 3)) on grid $(grids[gi]); lower `s0` or `smin`."
            gaveup = true
            break
        elseif pendinggrow || gi < length(grids)
            # A point carried from a coarser grid which has not converged on
            # this one retreats its drive here and walks back up, since
            # neither a growth past this grid nor a fold diagnosis may rest
            # on a drive established with a different truncation (`x` is the
            # carried point). A stall of a point converged on this grid,
            # below the finest, is its own solvability boundary: the grid
            # grows at the converged drive, and the point retreats on the
            # new grid if it does not reconverge there.
            factors, what = if pendinggrow
                (0.9, 0.8, 0.65, 0.5), "retreat on "
            else
                state = gridstate(grids[gi+1])
                x = stagedembed(out, state.freq.modes)
                gi += 1
                (1.0, 0.9, 0.8, 0.65, 0.5), "grow -> "
            end
            reconverged = false
            for f in factors
                t0 = time_ns()
                re = solve(state, f*s, x, false)
                last = re
                record(re, grids[gi], s, f*s, :grow,
                    re.solverinfo.converged, (time_ns() - t0)/1e9)
                verbose && println("staged: ", what, grids[gi], " at s=",
                    round(f*s, digits = 4), " |F|=",
                    round(re.solverinfo.finalresidual, sigdigits = 2),
                    re.solverinfo.converged ? "" : " STALL")
                if re.solverinfo.converged
                    out = re
                    s = f*s
                    x = vec(reshape(Array(out.nodeflux), out.Nmodes, :))
                    reconverged = true
                    break
                end
            end
            if !reconverged
                warnnotconverged && @warn lazy"the carried point did not reconverge on grid $(grids[gi]) even at half its drive; the truncation boundaries of the ladder are too far apart. Add an intermediate grid."
                gaveup = true
                break
            end
            pendinggrow = false
            ds = s0/2
        else
            # The continuation branch ends between s and starget. Before
            # declaring the operating point nonexistent, attempt the full
            # drive directly from the last converged point: coexisting
            # branches are common in this regime, and a plain damped solve
            # may reach one. One bounded attempt, recorded either way.
            if starget < 1.0 && !isnothing(x)
                t0 = time_ns()
                jump = solve(state, 1.0, x, true)
                last = jump
                record(jump, grids[gi], s, 1.0, :final,
                    jump.solverinfo.converged, (time_ns() - t0)/1e9)
                verbose && println("staged: branch-end jump to s=1.0 |F|=",
                    round(jump.solverinfo.finalresidual, sigdigits = 2),
                    jump.solverinfo.converged ? "" : " STALL")
                if jump.solverinfo.converged
                    out = jump
                    break
                end
            end
            warnnotconverged && @warn lazy"no harmonic balance solution was found at the requested drive: the source continuation converged at $(round(s, digits = 4)) of the requested amplitudes on the finest grid $(grids[end]) and stalled at $(round(starget, digits = 4)), and a direct solve at full drive from the last converged point also failed. This is what a fold of the solution branch between those two drives (the self oscillation threshold) looks like, but a failed search is not a proof that no operating point exists: a tighter continuation, another initial point or a different method may still reach one."
            fold = s
            gaveup = true
            break
        end
    end
    # a schedule which gave up returns its last attempt, marked not
    # converged whatever that attempt was, with the whole walk recorded;
    # one which did not has its outcome checked, once
    if gaveup
        out = last
    elseif warnnotconverged
        checkjunctioncurrents(out, psc)
    end
    # the diagnostics are the whole walk, one StagedStageInfo per attempt
    # with its inner solver records; every schedule records a first
    # attempt, whose solve holds the initial residual
    si = out.solverinfo
    F0 = stagerecords[1].inner[1].normresidual[1]
    newsi = SolverInfo(stagerecords, F0, si.finalresidual,
        !gaveup && si.converged, fold)
    vals = Any[getfield(out, f) for f in fieldnames(typeof(out))]
    vals[findfirst(==(:solverinfo), fieldnames(typeof(out)))] = newsi
    out = typeof(out)(vals...)
    return out
end
