# =========================================================================
# method = :staged -- source continuation on an adaptively grown grid.
# =========================================================================

"""
    StagedStageInfo

One attempted stage of [`stagedhbnlsolve`](@ref), stored in
`solverinfo.stages` of the returned solution -- every attempt appears, in
the order it ran, including stalled steps and growth retreats, so the
whole continuation walk can be examined afterwards.

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

The default coarse-to-fine harmonic grid ladder for [`stagedhbnlsolve`](@ref):
repeated halving of every dimension down to two harmonics, finest last.
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

# embed a converged solution into a larger grid's initial guess by matching
# mode tuples; modes the small grid does not carry start at zero
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
    stagedhbnlsolve(w, Nharmonics, sources, circuit, circuitdefs; kwargs...)

Source continuation on an adaptively grown harmonic grid, reached through
`hbnlsolve(...; method = :staged)`.

Near a critical drive the Newton basin is small and the iteration count
large, so those iterations are spent where they are cheap: the drive is
climbed in warm started steps with only a small set of harmonics retained
as unknowns (the nonlinearity is always evaluated on the full transform
grid; see `grids` below), and each larger retained set is warm started
from the last by matching mode tuples. The
schedule is adaptive in both directions because every truncation has its own
solvability boundary, and the boundaries are not monotone in the grid
(measured: the (2,2) truncation of a 64 junction RPM at strong two tone
drive converges at 98.4% of the drive while (4,2) does not): a stalled drive
step is first halved, a stall at the minimum step grows the grid at the
current converged drive, and a carried point that fails to reconverge after
growth retreats the drive on the new grid until it converges. Interior
points converge only to `interiorftol` under a small iteration budget --
they exist to keep the iterate inside the basin -- and the single expensive
solve, the finest grid at full drive, starts inside the basin and gets the
caller's `ftol` and `iterations`.

A point carried to the finest grid at full drive that stalls there without
ever having converged on that grid is not diagnosed as a fold: the drive is
retreated on the finest grid and climbed back, exactly as a mid-ladder
growth retreats, since the stalled drive value was established on different
physics. Only a stall from a point converged on the finest grid itself
brackets a fold: between the
last converged drive fraction and the stalled one the solution branch ends
(the self oscillation threshold), and the error says so, with the bracket.
Distinguishing that from a merely hard problem is something no cold started
method can do. Measured on a 128 junction RPM at equal 1.6 uA two tone
drive -- reported "unreachable" by every direct method -- the walk brackets
the fold at 91.9-93.4% of the drive: the operating point does not exist.

# Keywords (through `stagedkwargs`)
- `grids = defaultgridladder(Nharmonics)`: the coarse-to-fine ladder of
    RETAINED harmonic caps; the finest entry must equal `Nharmonics`.
    Every stage evaluates the nonlinearity on the full `Nharmonics` rdft
    transform grid -- the transform sets the aliasing of the nonlinear
    products and costs only N*log(N) -- and the ladder controls
    `maxharmonics`, which modes are retained as unknowns, so the linear
    solves are what shrink. Stages therefore share one alias-consistent
    grid, and a stage's solvability boundary reflects the physics rather
    than the aliasing artifacts of a small transform.
- `s0 = 0.5`: the first drive fraction attempted.
- `smin = 0.02`: the minimum drive step; a stall below it grows the grid.
- `interiorftol = 1e-7`, `interioriterations = 60`: tolerance and Newton
    budget of the interior points. The budget is deliberately small: a
    stalled probe is evident within tens of iterations, and interior stalls
    are the pure overhead of the walk (measured accepted interior points
    take 4-24 iterations, 51 near a pole).
- `innermethod = :newtonkrylov`: the solver for every stage.
- `interiorescalation = false`: whether interior stage solves may escalate
    their preconditioner to the full Jacobian. Off by default: an interior
    probe exists only to produce a cheap warm start, escalating inside one
    pays a full factorization for a throwaway point, and at high tone
    counts that factorization may not fit in memory at all (measured: a
    527 mode three tone grid ran out of 30 GB the moment an interior probe
    escalated, while the whole staged walk to that point took 19 seconds).
    A probe that fails without escalation is simply a stall, which the
    schedule answers with a smaller step. The final solve keeps the
    caller's escalation behavior.
- `maxattempts = 60`: bound on the total number of stage solves.
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

    # Every stage shares the FULL rdft box `Nharmonics`: the transform
    # grid sets the aliasing of the nonlinear products and costs only
    # N*log(N), so shrinking it with the stage would make each stage's
    # solvability boundary reflect its own aliasing artifacts rather than
    # the physics. The ladder controls `maxharmonics` -- which modes are
    # RETAINED as unknowns -- so the linear solves stay small while every
    # stage evaluates cos(phi(t)) on the same alias-consistent grid.
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
        # the first solve after carrying a full-drive point to a larger
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
                # full drive reached on a coarse grid: grow
                x = stagedembed(out, modesof(grids[gi+1]))
                gi += 1
                pendinggrow = true
            else
                x = vec(reshape(Array(out.nodeflux), out.Nmodes, :))
                ds = min(2*ds, 1.0 - s)
            end
        elseif starget - s > smin
            # the effective step, not `ds`: a target capped at full drive
            # must not be re-attempted identically while `ds` halves above
            # the cap
            ds = (starget - s)/2
        elseif gi < length(grids)
            # the current grid's own pole: grow at the converged drive,
            # retreating the drive on the new grid if the carried point does
            # not reconverge there (the boundaries are not monotone)
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
            # The stalled point was CARRIED to the finest grid from a
            # coarser one and has never converged here, so a fold diagnosis
            # would rest on a drive value established on different physics.
            # Retreat the drive on the finest grid and walk back up, exactly
            # as a mid-ladder growth retreats on its new grid; only a stall
            # from a point converged on the finest grid itself brackets a
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
            # declaring the operating point nonexistent, attempt the
            # requested drive directly from the last converged point:
            # coexisting branches are common in this regime, and the end of
            # the *continuation* branch does not preclude a solution at
            # full drive on another branch reachable by a plain damped
            # solve from here. One bounded attempt; its record appears in
            # the diagnostics either way.
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
    # the returned diagnostics are the whole walk: one StagedStageInfo per
    # attempt, each carrying its inner solver records, in place of just the
    # final solve's stages
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

