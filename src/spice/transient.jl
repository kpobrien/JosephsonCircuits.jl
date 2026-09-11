# WRspice as a back end of the transient: the compiled problem is written
# as a netlist with its drives as sources, run through the wrspice
# executable, and the rawfile read back into the package's solution type,
# so that the demodulation and the I/Q measurements read either solver's
# solution the same way.

"""
    WRspice(; executable = nothing, dphimax = 0.01, jjaccel = true,
        maxdata = 2e9)

WRSPICE as the stepping rule of [`transientsolve`](@ref), for
cross checking the package's own rules against an independent simulator.
The problem's circuit is written as a WRSPICE netlist with each junction
an instance of the `jj` model (see [`exportnetlist`](@ref)) and each
drive a current source — a constant one for a waveform given as a
number, and a piecewise linear one through the waveform's samples on the
grid of `dt` for a callable — and the simulation runs through the
`wrspice` executable in a temporary directory and is read back as a
[`TransientSolution`](@ref) holding the times, the port voltages and the
port waves, and with `record = :phases` the junction phases from the
phase nodes of the `jj` instances, as `(junction, time)` on the
junction axis every solver uses, in the branch orientation of the
package. [`transientdemodulate`](@ref) and [`transientiq`](@ref) then read
the solution as they read one stepped by the package.

A run keeps those traces and no others: the terminals of the ports, and
the phase nodes of the junctions when they are recorded. WRSPICE would
otherwise write every node of the circuit at every print point, which is
a gigabyte of a record nothing reads for every thousand nodes and every
hundred thousand points.

WRSPICE integrates with its own adaptive internal steps: `dt` is the
grid the sources are sampled on and the output is printed on, and
`dphimax`, the largest junction phase change WRSPICE allows per internal
step, controls the accuracy in place of the Newton tolerances, which do
not apply. `jjaccel` selects WRSPICE's accelerated convergence testing
for Josephson circuits, and `maxdata` its limit on the exported data in
kilobytes. `executable` is the path of the `wrspice` executable, or
`nothing` for the one [`wrspice_cmd`](@ref) finds; loading the
XicTools_jll package provides one.

Supported are the circuits of [`transientproblem`](@ref) without
scattering blocks, with the sinusoidal junction relation only, from the
zero initial state under WRSPICE's own initial conditions, with
`record = :ports` or `:phases`. An ideal [`TransmissionLine`](@ref) is
written as the SPICE lossless line element; the wave record `linewaves`
stays empty, since WRSPICE keeps the line histories to itself. The
junctions of every run share one `jj` model, whose subgap loss is set
as small as the model allows but is not zero, where the package's
junctions are lossless. The solution's final state fields are NaN,
since WRSPICE does not hand over a state, and the tangent, the adjoint
and the noise need a solution of the package's own rules.

# Examples
```julia
p = transientproblem(circuit, circuitdefs;
    sources = [TransientSource(1, t -> Ip*sin(wp*t))])
native = transientsolve(p, (0.0, 100e-9); dt = 1e-12)
spice = transientsolve(p, (0.0, 100e-9); dt = 1e-12, method = WRspice())
transientdemodulate(native, 2, wp/(2*pi)), transientdemodulate(spice, 2, wp/(2*pi))
```
"""
struct WRspice <: AbstractTransientIntegrator
    executable::Any
    dphimax::Float64
    jjaccel::Bool
    maxdata::Float64
end

function WRspice(; executable = nothing, dphimax::Real = 0.01,
        jjaccel::Bool = true, maxdata::Real = 2e9)
    isfinite(dphimax) && 0 < dphimax <= pi ||
        throw(ArgumentError("dphimax is the largest junction phase change per internal step, in radians between zero and pi."))
    isfinite(maxdata) && maxdata >= 1e3 ||
        throw(ArgumentError("maxdata is the WRSPICE limit on the exported data in kilobytes, at least 1e3."))
    return WRspice(executable, Float64(dphimax), jjaccel, Float64(maxdata))
end

# The transient through WRspice, behind `transientsolve(p, tspan; dt,
# method = WRspice())`. Everything the run needs is on the problem: the
# compiled circuit and its resolved values for the netlist, the drives
# with their waveforms for the sources, and the port data for the waves.
function wrspicetransient(p::TransientProblem, tspan, method::WRspice; dt,
        saveevery, record, initialstate, backend, linearsolver,
        factorization, reuse, checkpointevery, rtol, atol, maxiters, maxsteps)
    isempty(p.blocks) || throw(ArgumentError(
        "WRSPICE has no element for a scattering block; give the circuit as lumped elements and transmission lines only."))
    isnothing(p.relations) || throw(ArgumentError(
        "the WRSPICE jj model is the sinusoidal junction; a circuit with another current phase relation cannot run through it."))
    backend isa CPU || throw(ArgumentError("a WRSPICE run lives on the host; backend must be CPU()."))
    isnothing(linearsolver) || throw(ArgumentError("WRSPICE solves its own steps; it takes no linearsolver."))
    isnothing(factorization) || throw(ArgumentError("WRSPICE solves its own steps; it takes no factorization."))
    isnothing(reuse) || throw(ArgumentError("a WRSPICE run keeps no workspace to reuse."))
    checkpointevery == 0 || throw(ArgumentError("a WRSPICE run records no checkpoints."))
    record in (:ports, :phases) || throw(ArgumentError(
        "record must be :ports or :phases for WRspice, which does not hand over the package's states."))
    initialstate == transientstate(p) || throw(ArgumentError(
        "WRSPICE integrates from the zero state under its own initial conditions; shape the drive with a ramp instead of an initial state."))
    t0, tf, nsteps, h = transientgrid(tspan, dt, maxsteps, saveevery,
        maxiters, rtol, atol, nothing, nothing)
    nsteps % saveevery == 0 || throw(ArgumentError(
        lazy"the $(nsteps) steps of the grid are not a whole number of saves of every $(saveevery); choose dt or saveevery so they divide."))
    printstep = h*saveevery
    nsaved = nsteps ÷ saveevery

    input, junctions = wrspiceinput(p, method, t0, h, nsteps, printstep, record)
    executable = isnothing(method.executable) ? wrspice_cmd() : method.executable
    out = spice_run(input, executable)
    times, voltage, phases = wrspiceread(out, p, junctions, t0, tf, h,
        saveevery, nsteps, record)

    # the port waves, from the voltages, the terminations the ports own
    # and the drive currents, as `portwaves!` forms them on a step
    nports = length(p.portimpedances)
    current = zeros(nports, length(times))
    for d in p.drives
        d.portindex > 0 || continue
        for (j, t) in enumerate(times)
            current[d.portindex, j] += d.current(t)
        end
    end
    incident = similar(voltage)
    outgoing = similar(voltage)
    for k in 1:nports
        z = p.portimpedances[k]
        for j in eachindex(times)
            v = voltage[k, j]
            i = current[k, j] - p.portconductances[k]*v
            incident[k, j] = (v + z*i)/(2*sqrt(z))
            outgoing[k, j] = (v - z*i)/(2*sqrt(z))
        end
    end

    N = length(p)
    return TransientSolution(p, method, h, times, voltage, incident, outgoing,
        phases, nothing, nothing, nothing, nothing, nothing,
        zeros(N), zeros(N), fill(NaN, N), fill(NaN, N),
        nothing, nothing, nothing, (; steps = length(times) - 1))
end

# the input of a run: the netlist, the sources and the control block
function wrspiceinput(p::TransientProblem, method::WRspice, t0, h, nsteps,
        printstep, record::Symbol = :ports)
    n = exportnetlist(p.circuit, p.matrices.vvn; jj = true)
    names = p.circuit.nodenames
    for jn in n.junctions
        jn.phasenode in names && throw(ArgumentError(
            lazy"the net $(jn.phasenode) collides with the phase node WRSPICE gives a junction; leave the integers from the node count upward free as net names."))
    end
    lines = String[n.netlist]
    push!(lines, "* the drives and the constant sources")
    for (k, d) in enumerate(p.drives)
        a, b = wrspicesourcenodes(p, k)
        if d.waveform isa Number
            push!(lines, "isrcd$(k) $(a) $(b) $(Float64(d.waveform))")
        else
            push!(lines, wrspicepwl("isrcd$(k)", a, b, d.current, t0, h, nsteps))
        end
    end
    for r in 1:p.Nnodal
        iszero(p.constantcurrent[r]) && continue
        push!(lines, "isrcc$(r) 0 $(names[r + 1]) $(p.constantcurrent[r])")
    end
    push!(lines, "")
    push!(lines, ".tran $(printstep) $(h*nsteps) uic")
    push!(lines, "")
    push!(lines, ".control")
    push!(lines, "set maxdata=$(method.maxdata)")
    method.jjaccel && push!(lines, "set jjaccel=1")
    push!(lines, "set dphimax=$(method.dphimax)")
    saved = wrspicesaved(p, n.junctions, record)
    isempty(saved) || push!(lines, "save " * join(saved, " "))
    push!(lines, "run")
    push!(lines, "set filetype=binary")
    push!(lines, "write")
    push!(lines, ".endc")
    push!(lines, "")
    return join(lines, "\n"), n.junctions
end

# The traces a run keeps: the terminals of the ports, which carry the
# voltages and through them the waves, and with `record = :phases` the
# phase node of every junction. WRSPICE otherwise writes every node of
# the circuit at every print point, which on a long line is gigabytes of
# a record nothing reads. Ground is not a trace and repeated terminals
# are named once.
function wrspicesaved(p::TransientProblem, junctions, record::Symbol)
    names = p.circuit.nodenames
    saved = String[]
    for k in eachindex(p.portimpedances), idx in (p.portpositive[k], p.portnegative[k])
        idx == 0 && continue
        trace = "v($(names[idx + 1]))"
        trace in saved || push!(saved, trace)
    end
    if record === :phases
        for jn in junctions
            push!(saved, "v($(jn.phasenode))")
        end
    end
    return saved
end

# The two nodes of drive `k` from its injection column: the current is
# drawn from `a` and injected into `b`, which is how SPICE orients a
# current source, from `a` through the source to `b`. A terminal the
# column does not hold is ground.
function wrspicesourcenodes(p::TransientProblem, k::Int)
    names = p.circuit.nodenames
    a, b = "0", "0"
    for i in nzrange(p.injection, k)
        r = rowvals(p.injection)[i]
        if nonzeros(p.injection)[i] > 0
            b = names[r + 1]
        else
            a = names[r + 1]
        end
    end
    return a, b
end

# a piecewise linear current source through the samples of a waveform on
# the grid, in the simulation's own time, which starts at zero
function wrspicepwl(name, a, b, current, t0, h, nsteps)
    io = IOBuffer()
    print(io, name, " ", a, " ", b, " pwl(")
    for j in 0:nsteps
        j == 0 || print(io, j % 4 == 0 ? "\n+ " : " ")
        print(io, j*h, " ", current(t0 + j*h))
    end
    print(io, ")")
    return String(take!(io))
end

# The sign relating the phase WRSPICE reports for a junction written from
# node `n1` to node `n2`, both one indexed with ground first, to the
# branch orientation of the incidence matrix, whose columns leave ground
# out.
function wrspicephasesign(Rbn, b::Int, n1::Int, n2::Int)
    n1 > 1 && return Rbn[b, n1 - 1] > 0 ? 1.0 : -1.0
    return Rbn[b, n2 - 1] > 0 ? -1.0 : 1.0
end

# the times, the port voltages and the junction phases of a rawfile,
# validated against the requested grid, whose times are formed as the
# stepping rules form theirs so a decimated save holds the times of a
# full one
function wrspiceread(out::SpiceRaw, p::TransientProblem, junctions, t0, tf,
        h, saveevery, nsteps, record::Symbol)
    printstep = h*saveevery
    nsaved = nsteps ÷ saveevery
    haskey(out.values, "S") && haskey(out.values, "V") || throw(ArgumentError(
        "the WRSPICE rawfile holds no transient; expected a time axis and voltages."))
    rawtimes = vec(out.values["S"])
    length(rawtimes) == nsaved + 1 || throw(ArgumentError(
        lazy"WRSPICE returned $(length(rawtimes)) print points where $(nsaved + 1) were requested."))
    deviation = maximum(abs(rawtimes[j] - (j - 1)*printstep) for j in eachindex(rawtimes))
    deviation <= 1e-6*printstep || throw(ArgumentError(
        lazy"the print times of WRSPICE are off the requested grid by up to $(deviation) seconds."))
    times = [j == nsteps ? tf : t0 + j*h for j in 0:saveevery:nsteps]

    # the traces by name; the rawfile names a node voltage v(name)
    rows = Dict{String,Int}()
    for (i, name) in enumerate(out.variables["V"])
        rows[lowercase(name)] = i
    end
    V = out.values["V"]
    function noderow(idx::Int)
        idx == 0 && return 0
        key = lowercase("v($(p.circuit.nodenames[idx + 1]))")
        if !haskey(rows, key)
            held = join(sort!(collect(keys(rows))), ", ")
            throw(ArgumentError(lazy"the WRSPICE output has no trace $(key); it holds $(held)."))
        end
        return rows[key]
    end

    nports = length(p.portimpedances)
    voltage = zeros(nports, length(times))
    for k in 1:nports
        rp = noderow(p.portpositive[k])
        rn = noderow(p.portnegative[k])
        for j in eachindex(times)
            voltage[k, j] = (rp == 0 ? 0.0 : V[rp, j]) - (rn == 0 ? 0.0 : V[rn, j])
        end
    end

    record === :phases || return times, voltage, nothing
    # each phase node onto the junction axis of the solvers, the nonzero
    # entries of the branch inductance vector, in the branch orientation,
    # as `(junction, time)` at the saved times
    Ljb = p.matrices.Ljb
    phases = zeros(length(Ljb.nzind), length(times))
    for jn in junctions
        n1, n2 = p.circuit.nodeindices[1, jn.index], p.circuit.nodeindices[2, jn.index]
        b = p.graph.edge2indexdict[(n1, n2)]
        row = searchsortedfirst(Ljb.nzind, b)
        key = lowercase("v($(jn.phasenode))")
        haskey(rows, key) || throw(ArgumentError(
            lazy"the WRSPICE output has no phase trace $(key) for the junction $(p.circuit.componentnames[jn.index])."))
        r = rows[key]
        s = wrspicephasesign(p.graph.Rbn, b, n1, n2)
        for j in eachindex(times)
            phases[row, j] = s*V[r, j]
        end
    end
    return times, voltage, phases
end
