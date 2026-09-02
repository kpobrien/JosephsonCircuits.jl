# `method = :staged`: source continuation on an adaptively grown harmonic
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
- `action`: `:advance` (drive step on the current grid), `:grow` (first
    solve on a larger grid after carrying a converged point up), or
    `:final` (the full-drive solve on the finest grid).
- `accepted`: whether the attempt's result was kept as the new operating
    point (a stalled attempt is recorded but not accepted).
- `seconds`: wall time of the attempt, including the stage's system
    assembly.
- `finalresidual`: the residual norm the attempt ended at.
- `inner`: the inner solver's own stage records ([`IterationInfo`](@ref)),
    with their Krylov linear-solve diagnostics when the inner method is
    `:newtonkrylov`.
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
dimension down to two harmonics, coarsest first and `Nharmonics` last.
"""
function defaultgridladder(Nharmonics::NTuple{N,Int}) where {N}
    grids = [Nharmonics]
    g = Nharmonics
    while any(x -> x > 2, g)
        g = map(x -> max(2, cld(x, 2)), g)
        g == grids[end] && break
        push!(grids, g)
    end
    return reverse(grids)
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
    stagedhbnlsolve(w, Nharmonics, sources, circuit, circuitdefs;
        grids = defaultgridladder(Nharmonics), s0 = 0.5, smin = 0.02,
        interiorftol = 1e-7, interioriterations = 60,
        innermethod = :newtonkrylov, interiorescalation = false,
        maxattempts = 60, verbose = false, kwargs...)

Source continuation on an adaptively grown harmonic grid, reached through
`hbnlsolve(...; method = :staged, stagedkwargs = (; ...))`. `kwargs` are
the keywords of [`hbnlsolve`](@ref) (`iterations`, `ftol`,
`maxharmonics`, `maxintermodorder`, `dc`, `odd`, `even`, `symfreqvar`,
`sorting`, `keyedarrays`, `sensitivitynames`, `returnoperatingpoint`,
`krylovcouplingmodes`, `krylovrecycle`, `krylovharvest`, `krylovkwargs`,
`factorization`, `backend`, `precision`), which are forwarded to every
stage.

Near a critical drive the Newton basin is small and the iteration count
large, so those iterations are spent where they are cheap. The drive is
climbed in warm started steps with only a small set of harmonics retained
as unknowns, and each larger retained set is warm started from the last by
matching mode tuples. The nonlinearity is always evaluated on the full
`Nharmonics` transform grid, so that every stage sees the same aliasing of
the nonlinear products; the ladder only controls `maxharmonics`, the modes
retained as unknowns, so it is the linear solves which shrink.

The schedule adapts in both directions, because each truncation has its
own solvability boundary and the boundaries are not monotone in the grid:
a stalled drive step is halved; a stall at the minimum step grows the grid
at the current converged drive; and a carried point which fails to
reconverge after growth retreats the drive on the new grid until it
converges. Interior points converge only to `interiorftol` under a small
iteration budget, since they exist to keep the iterate inside the basin,
and the one expensive solve, the finest grid at full drive, starts inside
the basin with the caller's `ftol` and `iterations`.

A point carried to the finest grid which stalls there without ever having
converged on that grid is not diagnosed as a fold: the drive is retreated
on the finest grid and climbed back. Only a stall from a point converged on
the finest grid itself brackets a fold, the end of the solution branch
(the self oscillation threshold) between the last converged drive fraction
and the stalled one. Before reporting it, one plain damped solve at the
full drive is attempted from the last converged point, since a coexisting
branch may reach it; failing that, an error states the bracket.

# Keywords
- `grids = defaultgridladder(Nharmonics)`: the coarse to fine ladder of
    retained harmonic caps, whose last entry must equal `Nharmonics`.
- `s0 = 0.5`: the first drive fraction attempted.
- `smin = 0.02`: the minimum drive step; a stall below it grows the grid.
- `interiorftol = 1e-7`, `interioriterations = 60`: the tolerance and the
    Newton budget of the interior points. The budget is small on purpose: a
    stalled probe is evident within tens of iterations, and interior stalls
    are the overhead of the walk.
- `innermethod = :newtonkrylov`: the solver of every stage.
- `interiorescalation = false`: whether interior stage solves may escalate
    their preconditioner to the full Jacobian. Off by default, because an
    interior probe exists only to produce a cheap warm start, and at a high
    tone count the full factorization may not fit in memory; a probe which
    fails without escalation is a stall, which the schedule answers with a
    smaller step. The final solve keeps the caller's escalation behavior.
- `maxattempts = 60`: a bound on the total number of stage solves.
- `verbose = false`: print one line per stage solve.
"""
function stagedhbnlsolve(w::NTuple{N,Number}, Nharmonics::NTuple{N,Int},
    sources, circuit, circuitdefs;
    grids::Vector{NTuple{N,Int}} = defaultgridladder(Nharmonics),
    s0 = 0.5, smin = 0.02, interiorftol = 1e-7,
    interioriterations::Integer = 60, innermethod::Symbol = :newtonkrylov,
    interiorescalation::Bool = false,
    maxattempts::Integer = 60, verbose::Bool = false,
    iterations = 1000, maxharmonics::NTuple{N,Number} = Nharmonics,
    maxintermodorder = Inf, dc::Bool = false, odd::Bool = true,
    even::Bool = false, ftol = 1e-8, symfreqvar = nothing,
    sorting = :number, keyedarrays::Bool = true,
    sensitivitynames::Vector{String} = String[],
    returnoperatingpoint::Bool = false, krylovcouplingmodes = :none,
    krylovrecycle::Integer = 0, krylovharvest::Integer = 8,
    krylovkwargs::NamedTuple = (;), factorization = nothing,
    backend = CPU(), precision::Type{<:AbstractFloat} = Float64) where {N}

    isempty(grids) && throw(ArgumentError("`grids` must not be empty."))
    grids[end] == Nharmonics || throw(ArgumentError(
        lazy"the finest grid $(grids[end]) must equal `Nharmonics` = $(Nharmonics)."))
    0 < s0 <= 1 || throw(ArgumentError(lazy"`s0` = $(s0) must be in (0, 1]."))
    innermethod === :staged && throw(ArgumentError(
        "`innermethod` must be a non-staged method."))

    # Every stage uses the full transform grid `Nharmonics`, so the
    # aliasing of the nonlinear products is the same on every stage; the
    # ladder only sets `maxharmonics`, the modes retained as unknowns.
    modesof(grid) = removeconjfreqs(truncfreqs(calcfreqsrdft(Nharmonics);
        dc = dc, odd = odd, even = even,
        maxintermodorder = maxintermodorder,
        maxharmonics = map(min, maxharmonics, grid))).modes
    scaled(s) = [(mode = t.mode, port = t.port, current = s*t.current)
        for t in sources]
    solve(grid, s, x0, final) = hbnlsolve(w, Nharmonics, scaled(s), circuit,
        circuitdefs; dc = dc, odd = odd, even = even,
        maxintermodorder = maxintermodorder,
        maxharmonics = map(min, maxharmonics, grid),
        method = innermethod, x0 = x0, symfreqvar = symfreqvar,
        sorting = sorting,
        keyedarrays = final ? keyedarrays : false,
        sensitivitynames = final ? sensitivitynames : String[],
        returnoperatingpoint = final ? returnoperatingpoint : false,
        krylovcouplingmodes = krylovcouplingmodes,
        krylovrecycle = krylovrecycle, krylovharvest = krylovharvest,
        krylovkwargs = (final || interiorescalation) ? krylovkwargs :
            merge((; krylovescalate = typemax(Int)), krylovkwargs),
        factorization = factorization,
        backend = backend, precision = precision,
        ftol = final ? ftol : interiorftol,
        iterations = final ? iterations : interioriterations)

    gi = 1
    s = 0.0             # last converged drive fraction on the current grid
    ds = s0
    x = nothing         # its solution, in the raw nodeflux layout
    out = nothing
    attempts = 0
    stagerecords = AbstractStageInfo[]
    pendinggrow = false
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
        attempts > maxattempts && error(
            lazy"the staged schedule did not converge in `maxattempts` = $(maxattempts) stage solves; last converged drive fraction $(s) on grid $(grids[gi]).")
        starget = min(1.0, s + ds)
        final = gi == length(grids) && starget >= 1.0
        t0 = time_ns()
        cand = solve(grids[gi], starget, x, final)
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
                x = stagedembed(out, modesof(grids[gi+1]))
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
        elseif gi < length(grids)
            # the current grid's own solvability boundary: grow the grid at
            # the converged drive, retreating the drive on the new grid if
            # the carried point does not reconverge there
            isnothing(x) && error(
                "the first stage stalled at its first drive step; lower `s0`.")
            bigmodes = modesof(grids[gi+1])
            bigx = stagedembed(out, bigmodes)
            gi += 1
            reconverged = false
            for f in (1.0, 0.9, 0.8, 0.65, 0.5)
                t0 = time_ns()
                re = solve(grids[gi], f*s, bigx, false)
                record(re, grids[gi], s, f*s, :grow,
                    re.solverinfo.converged, (time_ns() - t0)/1e9)
                verbose && println("staged: grow -> ", grids[gi], " at s=",
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
            reconverged || error(
                lazy"the carried point did not reconverge on grid $(grids[gi]) even at half its drive; the truncation boundaries of the ladder are too far apart. Add an intermediate grid.")
            ds = s0/2
        elseif pendinggrow
            # The stalled point was carried to the finest grid from a
            # coarser one and has never converged here, so a fold diagnosis
            # would rest on a drive established with a different truncation.
            # Retreat the drive on the finest grid and walk back up; only a
            # stall from a point converged on the finest grid brackets a
            # fold.
            isnothing(x) && error(
                "the first stage stalled at its first drive step; lower `s0`.")
            reconverged = false
            for f in (0.9, 0.8, 0.65, 0.5)
                t0 = time_ns()
                re = solve(grids[gi], f*s, x, false)
                record(re, grids[gi], s, f*s, :grow,
                    re.solverinfo.converged, (time_ns() - t0)/1e9)
                verbose && println("staged: retreat on ", grids[gi], " at s=",
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
            reconverged || error(
                lazy"the carried point did not reconverge on the finest grid $(grids[gi]) even at half its drive; the truncation boundaries of the ladder are too far apart. Add an intermediate grid.")
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
                jump = solve(grids[gi], 1.0, x, true)
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
            error(lazy"no harmonic balance solution was found at the requested drive: the source continuation converged at $(round(s, digits = 4)) of the requested amplitudes on the finest grid $(grids[end]) and stalled at $(round(starget, digits = 4)), so the solution branch ends between them (a fold, the self oscillation threshold), and a direct solve at full drive from the last converged point also failed. No operating point at the requested drive is reachable from below.")
        end
    end
    # the diagnostics are the whole walk, one StagedStageInfo per attempt
    # with its inner solver records
    si = out.solverinfo
    F0 = try
        stagerecords[1].inner[1].normresidual[1]
    catch
        si.initialresidual
    end
    newsi = SolverInfo(stagerecords, F0, si.finalresidual, si.converged,
        si.sourcefold)
    vals = Any[getfield(out, f) for f in fieldnames(typeof(out))]
    vals[findfirst(==(:solverinfo), fieldnames(typeof(out)))] = newsi
    out = typeof(out)(vals...)
    return out
end

