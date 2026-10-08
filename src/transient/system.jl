# The circuit in time, on the state harmonic balance solves for: the node
# fluxes in units of the reduced flux quantum, the auxiliary branch
# currents of the mutually coupled inductors, and the gauge rows of the
# floating inductive subnetworks, at one mode. The matrices are the ones
# `numericmatrices` builds and the augmentation is the one `hbnlsolve`
# folds into its linear term, so the transient adds integration in time to
# the circuit machinery and no assembly of its own. The equations are
#
#     C phi'' + G phi' + invL phi + Ic sin(phi_b) + A_mna = I(t),
#
# scaled row by row by `Lscale/phi0` as harmonic balance scales them, so
# that a junction enters as `Lscale/Lj` and the drive as `Lscale*I/phi0`.
# The scale is a property of the problem, not of a step, so that a state
# means the same thing at every step size.

const TransientWaveform = FunctionWrapper{Float64, Tuple{Float64}}

transientwaveform(w::Number) = TransientWaveform(Returns(Float64(w)))
transientwaveform(w) = TransientWaveform(w)

"""
    TransientSource(target, current)

A real instantaneous current in Amperes, a number or a callable
`current(t)` of the time in seconds. An integer `target` is a port
number: a Norton source injects the current into the port's first
(positive) terminal, and the port's termination is part of the compiled
circuit already. A string or symbol names a `CurrentSource` component,
whose constant value the waveform replaces, and which drives its current
through itself from its first terminal to its second: it draws the
current from the node at its first terminal and delivers it to the node
at its second, the opposite sense of a port source. Several sources on
one target add. Unlike a harmonic balance source the callable
returns the physical waveform, not a Fourier coefficient, and it must be
deterministic, since the tangent and the adjoint evaluate it again on the
recorded grid.
"""
struct TransientSource{F}
    target::Union{Int,String}
    current::F
end
TransientSource(target::Integer, current) = TransientSource{typeof(current)}(Int(target), current)
TransientSource(target::AbstractString, current) = TransientSource{typeof(current)}(String(target), current)
TransientSource(target::Symbol, current) = TransientSource(String(target), current)

# A bound drive: the column of the injection matrix it scales, the port it
# drives (zero for a named source, whose current does not enter a port
# wave), and its waveform through a wrapper, so that a problem has one type
# whatever its waveforms are. The waveform object is kept for setup code
# that reads it.
struct TransientDrive
    portindex::Int
    current::TransientWaveform
    waveform::Any
end
TransientDrive(portindex::Integer, waveform) =
    TransientDrive(Int(portindex), transientwaveform(waveform), waveform)

"""
    BlockModulation

One modulated output of a pumped block realized in time: the output
matrix `C` over the block's states of the cosine (`quadrature = 1`) or
sine (`quadrature = 2`) filter of the harmonic `k`, whose output the
reflected wave carries multiplied by `2 cos(k wp t)` or `-2 sin(k wp t)`
and by the block's envelope.
"""
struct BlockModulation
    harmonic::Int
    quadrature::Int
    C::Matrix{Float64}
end

"""
    TransientBlock

A scattering parameter block as the transient realizes it: the real
constant scattering matrix `S` at the reference impedances `R` of its
ports, the node of the signal and the reference terminal of each port
(zero for ground), the index of the first of its port current unknowns,
and the block's path. Its ports carry the hybrid constitutive equation
of the linearized solver, `(I - S) R^(-1/2) v - (I + S) R^(1/2) i = 0`,
with `v` the port voltages and `i` the port currents entering the block
through the signal terminals, as auxiliary unknowns whose Kirchhoff
couplings enter the node equations; nothing is inverted, so a short or
an open is stamped exactly, and a scattering entry within roundoff of
one is snapped to it, so a fitted feedthrough on the unit circle
carries the exact zeros of its hybrid coefficients (see
[`snapscattering`](@ref)). A block whose matrix depends on frequency
carries its rational part as states, and a pumped block the states of
every filter of its harmonics, the unconverted output on `C` and the
converted ones modulated (see [`BlockModulation`](@ref)).
"""
struct TransientBlock
    definition::Any
    S::Matrix{Float64}
    R::Vector{Float64}
    signal::Vector{Int}
    ref::Vector{Int}
    auxbase::Int
    path::String
    # the rational part, `S(s) = S + C (s I - A)^(-1) B`, empty for a
    # constant block, and where its states begin in the states of all
    # blocks
    A::Matrix{Float64}
    B::Matrix{Float64}
    C::Matrix{Float64}
    zbase::Int
    # the modulated outputs of a pumped block (see
    # [`BlockModulation`](@ref)), its pump frequency and the envelope of
    # its conversion in time; empty, zero and nothing for a block which
    # does not convert
    modulations::Vector{BlockModulation}
    wp::Float64
    envelope::Any
end

# the weight of a modulated output at the time `t`: the modulation, and
# the envelope of the block's conversion
function modulationweight(b::TransientBlock, m::BlockModulation, t)
    env = isnothing(b.envelope) ? 1.0 : Float64(b.envelope(t))
    return m.quadrature == 1 ? 2env*cos(m.harmonic*b.wp*t) : -2env*sin(m.harmonic*b.wp*t)
end

"""
    TransientLine

An ideal lossless transmission line as the transient realizes it: its
characteristic impedance `Z`, its delay, and the signal and reference
node of each of its two ports (zero for ground). In time the line is the
method of characteristics: at each port the current into the line is
`v/Z - 2 q/sqrt(Z)`, a conductance `1/Z` at the line's own impedance
plus a current from the wave `q` that entered the far port a delay
earlier, and the wave leaving a port, `a = v/sqrt(Z) - q`, is kept as
the port's history for the far port to read a delay later. Nothing is
inverted and nothing resonates in the equations, since the round trips
that make the line's admittance singular in frequency are the recursion
through the history; a mismatch to the connected circuit is the shared
node. The delay must be at least one step, so that every wave a step
reads is accepted history, and the read interpolates the endpoint
history with the centered Lagrange stencil of up to six samples, a
quintic where the delay leaves three accepted samples past the query and
of lower order for a shorter delay.
"""
struct TransientLine
    Z::Float64
    delay::Float64
    signal::NTuple{2,Int}
    ref::NTuple{2,Int}
    path::String
end

"""
    FloatingBalance

The floating subnetworks of a [`TransientProblem`](@ref) which its drives
feed, judged whenever the drives are evaluated (`checkbalance`): a net
current into a subnetwork no element connects to ground has no path back.
`islands` indexes the problem's `floatingcomponents` and `nodes` names
their nodes; `net[k, r]` is the net injection of a unit current of drive
`k` into island `r`, a column for each island holding the drives which
feed it, two islands at most for a drive's two terminals, `constant[r]`
the net constant current into it and `magnitude[r]` the sum of the
magnitudes that net gathers, and the total may differ from zero by
`terms` roundings of the magnitudes it sums.
"""
struct FloatingBalance
    islands::Vector{Int}
    nodes::Vector{String}
    net::SparseMatrixCSC{Float64,Int}
    constant::Vector{Float64}
    magnitude::Vector{Float64}
    terms::Int
end

"""
    TransientProblem

A compiled circuit with its matrices at one mode, its modified nodal
analysis augmentation, its bound drives and its port data, ready to be
integrated in time by [`transientsolve`](@ref). Built by
[`transientproblem`](@ref).

# Fields
- `circuit`, `matrices`: the compiled circuit and its
    [`CircuitMatrices`](@ref) at one mode, with real values.
- `Nnodal`, `Naux`: the node flux unknowns and the auxiliary branch
    currents of the mutually coupled inductors; the state has
    `Nnodal + Naux` entries.
- `Lscale`: the inductance scale of the equations, the mean inductance
    of the circuit as harmonic balance uses at zero frequency, fixed for
    the problem so that a state is independent of the step: an auxiliary
    current `i` is stored as `Lscale*i/phi0`.
- `coupledbranches`, `floatingcomponents`, `gaugeindices`: the branches
    the augmentation promotes, the subnetworks no element connects to
    ground, whose flux offset is free, and the gauge rows fixing it.
- `inertialess`: the directions of the state that carry no inertia, one
    set of state indices per subnetwork of the capacitive graph no
    capacitor connects to ground, including every unknown no capacitor
    touches; the algebraic constraints a state must satisfy at the start
    lie along them.
- `algebraic`: the directions among them that carry no dissipation
    either, one set per subnetwork no capacitor or resistor connects to
    ground, each a union of inertialess ones; along them the equations
    constrain the flux alone, and the rate is what the constraint
    differentiated says.
- `directions`, `constraints`: the algebraic directions as the columns
    of a sparse matrix over the state, in the order of `algebraic`, and
    the constraints along them as the rows of one over the equations, the
    combinations of the equations which no rate enters.
- `injection`: the unscaled injection of a unit current of each drive
    into the node equations, one sparse column per drive, in the
    orientation of the drive.
- `drives`: the bound drives, in the order of `injection`'s columns.
- `constantcurrent`: the unscaled constant node current of the netlist's
    current sources not replaced by a drive.
- `balance`: the floating subnetworks the drives feed, whose net current
    is judged whenever the drives are evaluated, see
    [`FloatingBalance`](@ref).
- `ports`: the compiled ports in the order of their numbers, the order
    of every port axis of the transient, as of harmonic balance's.
- `portpositive`, `portnegative`, `portimpedances`, `portconductances`:
    the node of each port terminal (zero for ground), the reference
    impedance and the conductance of the termination the port owns, in
    the order of `ports`, for the port waves.
- `blocks`: the scattering parameter blocks with a realization in time,
    see [`TransientBlock`](@ref), whose port currents are auxiliary
    unknowns after the coupled inductor currents.
- `lines`: the ideal transmission lines, see [`TransientLine`](@ref),
    whose wave histories the solve keeps.
- `relations`: the current-phase relation of every junction, in the
    order of the junction rows of `RJ`, or `nothing` where every junction
    is sinusoidal.
- `C`, `G`, `L`, `lineE`: the scaled capacitance, conductance and
    augmented inverse inductance the transient steps, with the stamps of
    the blocks and the lines, and the lines' incidence, on the host.
- `RJ`, `lmolj`: the junction rows of the incidence and the
    coefficients `Lscale/Lj`, in the order of `relations`.
"""
struct TransientProblem
    circuit::CompiledCircuit
    matrices::CircuitMatrices
    Nnodal::Int
    Naux::Int
    Lscale::Float64
    coupledbranches::Vector{Int}
    floatingcomponents::Vector{Vector{Int}}
    gaugeindices::Vector{Int}
    inertialess::Vector{Vector{Int}}
    algebraic::Vector{Vector{Int}}
    directions::SparseMatrixCSC{Float64,Int}
    constraints::SparseMatrixCSC{Float64,Int}
    injection::SparseMatrixCSC{Float64,Int}
    drives::Vector{TransientDrive}
    constantcurrent::Vector{Float64}
    balance::FloatingBalance
    ports::Vector{CompiledPort}
    portpositive::Vector{Int}
    portnegative::Vector{Int}
    portimpedances::Vector{Float64}
    portconductances::Vector{Float64}
    blocks::Vector{TransientBlock}
    lines::Vector{TransientLine}
    # the current-phase relation of every junction, in the order of the
    # nonzero entries of `matrices.Ljb`, which is the order of the junction
    # rows of `RJ` and of `lmolj`. `nothing` when every one of them is the
    # sinusoidal Josephson relation.
    relations::Union{Nothing,JunctionRelations{Matrix{Float64},Vector{Int}}}
    # the scaled matrices and the junction incidence, built once for the
    # classification and for every system of the problem
    C::SparseMatrixCSC{Float64,Int}
    G::SparseMatrixCSC{Float64,Int}
    L::SparseMatrixCSC{Float64,Int}
    lineE::SparseMatrixCSC{Float64,Int}
    RJ::SparseMatrixCSC{Float64,Int}
    lmolj::Vector{Float64}
end

Base.length(p::TransientProblem) = p.Nnodal + p.Naux

# The problem `p` with its drives or its blocks replaced and every other
# field shared, the one place which lists the fields for such a copy.
function TransientProblem(p::TransientProblem; injection = p.injection, drives = p.drives,
        constantcurrent = p.constantcurrent, balance = p.balance, blocks = p.blocks)
    return TransientProblem(p.circuit, p.matrices, p.Nnodal, p.Naux, p.Lscale, p.coupledbranches,
        p.floatingcomponents, p.gaugeindices, p.inertialess, p.algebraic, p.directions, p.constraints,
        injection, drives, constantcurrent, balance, p.ports, p.portpositive, p.portnegative, p.portimpedances,
        p.portconductances, blocks, p.lines, p.relations, p.C, p.G, p.L, p.lineE, p.RJ, p.lmolj)
end

# a real value of a component, or an argument error naming it. The values
# of the table are already numbers, from `numericvalues`; a port impedance
# is checked here too and may not be
function transientreal(value, name)
    (value isa Number && !checkissymbolic(value)) || throw(ArgumentError(
        lazy"the component $(name) has the value $(value); the transient needs real, frequency independent values."))
    iszero(imag(value)) || throw(ArgumentError(
        lazy"the component $(name) has the complex value $(value); dissipation in time needs a real resistor, not an imaginary part."))
    v = Float64(real(value))
    (isfinite(v) || v == Inf) || throw(ArgumentError(
        lazy"the component $(name) has the nonfinite value $(value)."))
    return v
end

"""
    transientproblem(circuit, circuitdefs = Dict(); sources = ())

Compile a circuit for integration in time: the same compiler and
[`numericmatrices`](@ref) as harmonic balance, at one mode, with the
mutually coupled inductors promoted to auxiliary branch currents and the
floating inductive subnetworks gauge fixed as [`hbnlsolve`](@ref) does.
`sources` is a tuple or vector of [`TransientSource`](@ref)s; the
netlist's constant `CurrentSource` components keep their constant values
unless a source names them. The circuit may be a typed [`Circuit`](@ref)
or a compiled circuit.

Supported are real, constant resistors, capacitors, inductors, mutual
inductors, Josephson junctions and nonlinear inductors with a polynomial
current-phase relation (see [`PolynomialCPR`](@ref)), current sources and
ports. Frequency dependent or complex values are rejected, since they
need a causal realization in time, which a block may carry: a
[`ScatteringParameters`](@ref) block with a constant real matrix is
realized as it is, a [`RationalScattering`](@ref) block and a pumped
block fitted for time with `RationalScattering(block, npoles)` by their
states, see [`TransientBlock`](@ref), and an ideal
[`TransmissionLine`](@ref) by its delay, see [`TransientLine`](@ref);
circuits with blocks or lines step under [`GaussLegendre`](@ref). Any
other block is rejected.
"""
Base.@nospecializeinfer function transientproblem(circuit::CompilableCircuit,
        circuitdefs::AbstractDict = Dict{Symbol,Any}();
        sources = ())
    # compiled once for every kind of circuit and every waveform: the body
    # reads the compiled circuit and binds the drives through wrappers
    @nospecialize circuit sources
    psc = compile(circuit)::CompiledCircuit
    vvn = numericvalues(psc, circuitdefs)
    for k in eachindex(vvn)
        v = transientreal(vvn[k], psc.componentnames[k])
        psc.componenttypes[k] == :R || isfinite(v) || throw(ArgumentError(
            lazy"the component $(psc.componentnames[k]) has an infinite value; only a resistor may be an open."))
        vvn[k] = v
    end
    checkstaticstiffnessvalues(psc.componenttypes, vvn)
    nm = numericmatrices(psc, vvn; Nmodes = 1)
    for (name, A) in (("capacitance", nm.Cnm), ("conductance", nm.Gnm),
            ("inverse inductance", nm.invLnm))
        all(isfinite, nonzeros(A)) || throw(ArgumentError(
            lazy"the $(name) matrix has a nonfinite entry."))
    end
    # every capacitance nonnegative, so that the capacitance matrix, a
    # weighted Laplacian with the capacitances to ground on its diagonal,
    # is positive semidefinite, as the time stepping needs it
    for k in eachindex(vvn)
        psc.componenttypes[k] == :C && vvn[k] < 0 && throw(ArgumentError(
            lazy"the capacitor $(psc.componentnames[k]) has the negative value $(vvn[k]); a negative capacitance has no meaning in time."))
    end

    Nnodal = psc.Nnodes - 1
    coupledbranches = mnacoupledbranches(nm.Mb)
    # the auxiliary unknowns: the coupled inductor currents, then the
    # port currents of the scattering blocks
    blocks = transientblocks(psc, Nnodal + length(coupledbranches))
    lines = transientlines(psc)
    Naux = length(coupledbranches) + sum(b -> length(b.signal), blocks; init = 0)
    floatingcomponents = transientfloatingcomponents(psc, vvn)
    gaugeindices = calcdcgaugeindices(floatingcomponents, [0.0], 1)
    # the scale of the equations, which no step size enters: the package's
    # at zero frequency, the mean inductance, and for a circuit without
    # inductance the impedance scale times the impedance and mean
    # capacitance, or a picosecond
    Lscale = transientscale(psc, vvn, nm)
    C, G, L, lineE = transientlinearmatrices(nm, coupledbranches, psc.topology.Rbn, gaugeindices, blocks, lines, Lscale, Nnodal, Naux)
    inertialess, algebraic, directions, constraints =
        transientclassification(psc, vvn, G, L, Nnodal, length(coupledbranches), blocks)
    RJ, lmolj = junctionincidence(nm, Naux, Lscale)

    # the port terminals and terminations, for the port waves, in the
    # order of the ports' numbers
    ports = orderedports(psc)
    np = length(ports)
    portpositive = [p.positivenode - 1 for p in ports]
    portnegative = [p.negativenode - 1 for p in ports]
    portimpedances = [transientreal(nm.portimpedances[findfirst(==(ports[k].number), nm.portnumbers)],
        "the impedance of port $(ports[k].number)") for k in 1:np]
    all(>(0), portimpedances) || throw(ArgumentError("port reference impedances must be positive."))
    portconductances = [p.environment == 0 ? 0.0 : 1/vvn[p.environment] for p in ports]

    drives, injection, constantcurrent, balance = bindsources(psc, vvn, ports, portpositive, portnegative, Nnodal + Naux,
        floatingcomponents, sources)
    return TransientProblem(psc, nm, Nnodal, Naux, Lscale, coupledbranches,
        floatingcomponents, gaugeindices, inertialess, algebraic, directions, constraints,
        injection, drives, constantcurrent, balance, ports, portpositive, portnegative, portimpedances, portconductances, blocks, lines,
        calcjunctionrelations(psc.componenttypes, psc.nodeindices,
            psc.junctioncprs, psc.topology.edge2indexdict, nm.Ljb), C, G, L, lineE, RJ, lmolj)
end

# the junction rows of the incidence, over the state with its `Naux`
# auxiliary unknowns, and the coefficients `Lscale/Lj`, in the order of
# the junctions of the matrices
function junctionincidence(nm::CircuitMatrices, Naux::Int, Lscale::Float64)
    Rbnm = hcat(nm.Rbnm, spzeros(eltype(nm.Rbnm), size(nm.Rbnm, 1), Naux))
    Ljb = nm.Ljb
    return SparseMatrixCSC{Float64,Int}(Rbnm[Ljb.nzind, :]), Float64[Lscale/Ljb.nzval[i] for i in eachindex(Ljb.nzval)]
end

# The drives of `sources` bound to a compiled circuit with `n` unknowns
# and the ports `ports` with their terminals: a unit current of each
# source as a node injection, one column per source, a port's into its
# positive terminal and a named current source's through itself from its
# first terminal to its second, the constant current of the netlist's
# current sources no source replaced, from the values `vvn`, and the
# balance of the floating subnetworks. A subnetwork of `floating`, the
# node sets no element connects to ground, has no path back for a net
# current driven into it, which would flow through its gauge row, an
# inductor of the problem's scale to ground, and set its flux by that
# scale, and harmonic balance refuses the subnetwork itself. What counts
# is the sources' total: sources whose net currents into a subnetwork
# cancel, two drives of one waveform into it and out of it say, give it
# the forcing of one source across it. A subnetwork no drive feeds has its
# constant currents judged here, beyond the rounding of their sum; one a
# drive feeds is judged with its drives whenever they are evaluated
# (`checkbalance`). Compiled once for every collection of sources, whose
# waveforms the drives wrap.
Base.@nospecializeinfer function bindsources(psc::CompiledCircuit, vvn::Vector, ports::Vector{CompiledPort}, portpositive,
        portnegative, n::Int, floating::Vector{Vector{Int}}, @nospecialize(sources))
    drives = TransientDrive[]
    rows, cols, vals = Int[], Int[], Float64[]
    replaced = Set{Int}()
    for (k, source) in enumerate(sources)
        source isa TransientSource || throw(ArgumentError("sources must contain TransientSource objects."))
        if source.target isa Int
            q = portindex(ports, source.target)
            n1, n2 = portpositive[q], portnegative[q]
            push!(drives, TransientDrive(q, source.current))
        else
            c = get(psc.componentnamedict, source.target, 0)
            (c > 0 && psc.componenttypes[c] == :I) || throw(ArgumentError(
                lazy"$(source.target) does not name a CurrentSource of the circuit."))
            push!(replaced, c)
            n1, n2 = psc.nodeindices[2, c] - 1, psc.nodeindices[1, c] - 1
            push!(drives, TransientDrive(0, source.current))
        end
        n1 > 0 && (push!(rows, n1); push!(cols, k); push!(vals, 1.0))
        n2 > 0 && (push!(rows, n2); push!(cols, k); push!(vals, -1.0))
    end
    injection = sparse(rows, cols, vals, n, length(drives))
    # the constant current into each node and the magnitudes it sums
    constantcurrent, magnitude = zeros(n), zeros(n)
    for c in psc.currentsources
        c in replaced && continue
        n1, n2 = psc.nodeindices[2, c] - 1, psc.nodeindices[1, c] - 1
        n1 > 0 && (constantcurrent[n1] += vvn[c]; magnitude[n1] += abs(vvn[c]))
        n2 > 0 && (constantcurrent[n2] -= vvn[c]; magnitude[n2] += abs(vvn[c]))
    end
    # the net injection of each drive into each subnetwork, gathered from
    # the drives' terminals, a column for each subnetwork; a drive with
    # both terminals in one subnetwork feeds it nothing
    islandof = zeros(Int, n)
    for (r, island) in enumerate(floating), i in island
        islandof[i - 1] = r
    end
    ks, rs, vs = Int[], Int[], Float64[]
    for k in axes(injection, 2), q in nzrange(injection, k)
        r = islandof[rowvals(injection)[q]]
        r > 0 && (push!(ks, k); push!(rs, r); push!(vs, nonzeros(injection)[q]))
    end
    nets = dropzeros!(sparse(ks, rs, vs, length(drives), length(floating)))
    fed = Int[]; constants = Float64[]; magnitudes = Float64[]
    for (r, island) in enumerate(floating)
        net = sum(constantcurrent[i - 1] for i in island)
        scale = sum(i -> magnitude[i - 1], island)
        if !isempty(nzrange(nets, r))
            push!(fed, r); push!(constants, net); push!(magnitudes, scale)
        else
            nodes = join(psc.nodenames[island], ", ")
            abs(net) <= 2length(psc.currentsources)*eps(Float64)*scale || throw(ArgumentError(
                lazy"the constant current sources drive a net current into the nodes ($(nodes)), which no element connects to ground, so the current has no path back; connect them to ground (a resistor or a capacitor will do)."))
        end
    end
    balance = FloatingBalance(fed, [join(psc.nodenames[floating[r]], ", ") for r in fed],
        nets[:, fed], constants, magnitudes, 2*(length(psc.currentsources) + length(drives)))
    return drives, injection, constantcurrent, balance
end

# the net current of the drives, at the `values` they take at `t`, and of
# the constant sources into each floating subnetwork a drive feeds,
# refused beyond the rounding of its sum; each subnetwork reads the drives
# which feed it
function checkbalance(p::TransientProblem, values::AbstractVector, t)
    b = p.balance
    drives, nets = rowvals(b.net), nonzeros(b.net)
    for r in eachindex(b.islands)
        net, scale = b.constant[r], b.magnitude[r]
        for q in nzrange(b.net, r)
            x = nets[q]*values[drives[q]]
            net += x
            scale += abs(x)
        end
        abs(net) <= b.terms*eps(Float64)*scale || throw(ArgumentError(
            lazy"the sources drive a net current of $(net) A into the nodes ($(b.nodes[r])) at t = $(t) s, which no element connects to ground, so the current has no path back; balance the sources, or connect the nodes to ground (a resistor or a capacitor will do)."))
    end
    return nothing
end

# the index of the port numbered `number` among `ports`, the row of its
# trace
function portindex(ports::Vector{CompiledPort}, number)
    number isa Integer || throw(ArgumentError(lazy"a port is named by its number, not $(number)."))
    q = findfirst(port -> port.number == number, ports)
    isnothing(q) && throw(ArgumentError(lazy"there is no port $(number)."))
    return q
end

# the ideal lines of a compiled circuit, in compiled order
function transientlines(psc::CompiledCircuit)
    lines = TransientLine[]
    for cb in psc.scatteringblocks
        cb.definition isa LinearizedScattering && continue
        provider = cb.definition.provider
        provider isa TransmissionLineProvider || continue
        (isfinite(provider.Z0) && provider.Z0 > 0 && isfinite(provider.delay) && provider.delay >= 0) || throw(ArgumentError(
            lazy"the transmission line at $(cb.path) needs a finite positive impedance and a finite delay."))
        provider.delay > 0 || throw(ArgumentError(
            lazy"the transmission line at $(cb.path) has no length; write it as a through block, ScatteringParameters([0 1; 1 0]; zref)."))
        push!(lines, TransientLine(provider.Z0, provider.delay, (cb.signalnodes[1] - 1, cb.signalnodes[2] - 1),
            (cb.refnodes[1] - 1, cb.refnodes[2] - 1), cb.path))
    end
    return lines
end


"""
    snapscattering(S::AbstractMatrix)

The real scattering matrix `S` with every entry within `1e-12` of a
perfect open, short or isolation snapped to the exact `1`, `-1` or `0`,
which is how a block realized in time stamps it. The hybrid
coefficients `I - S` and `I + S` of a port whose feedthrough reaches
the unit circle then hold exact zeros where the algebra has them, so
the endpoint's rate system sees a zero row rather than the roundoff
residue of one, which its balancing would otherwise scale up into an
equation whose inverted singular value multiplies the residual by the
reciprocal of machine epsilon at every step. The snap moves an entry
by less than `1e-12`, far below any scattering the data resolves.
"""
function snapscattering(S::AbstractMatrix)
    snap = x -> abs(x) < 1e-12 ? zero(x) :
        abs(x - 1) < 1e-12 ? one(x) :
        abs(x + 1) < 1e-12 ? -one(x) : x
    return snap.(Matrix{Float64}(real.(S)))
end

# the blocks of a compiled circuit the transient realizes, in compiled
# order, their port currents laid out from `offset`
function transientblocks(psc::CompiledCircuit, offset::Int)
    blocks = TransientBlock[]
    zbase = 0
    for cb in psc.scatteringblocks
        def = cb.definition
        if def isa LinearizedScattering
            # a pumped block fitted for time: the states of every filter
            # in one realization, the unconverted output on the constant
            # rows, and the converted outputs modulated, its pump phase
            # folded into them
            realizedintime(def) || throw(ArgumentError(
                lazy"the pumped scattering block at $(cb.path) has no realization in time; fit it with RationalScattering(block, npoles)."))
            n = def.nports
            As, Bs = Matrix{Float64}[], Matrix{Float64}[]
            modulations = BlockModulation[]
            p0 = def.providers[1]
            push!(As, p0.A); push!(Bs, p0.B)
            outputs = Tuple{Int,Int,Matrix{Float64}}[]
            for (j, k) in enumerate(def.harmonics)
                j == 1 && continue
                p = def.providers[j]
                for (q, part) in enumerate((p.cosine, p.sine))
                    push!(As, part.A); push!(Bs, part.B)
                    push!(outputs, (k, q, part.C))
                end
            end
            nz = sum(size(A, 1) for A in As)
            A = zeros(nz, nz)
            B = zeros(nz, n)
            z = 0
            for (Ab, Bb) in zip(As, Bs)
                m = size(Ab, 1)
                A[z + 1:z + m, z + 1:z + m] .= Ab
                B[z + 1:z + m, :] .= Bb
                z += m
            end
            C = zeros(n, nz)
            C[:, 1:size(p0.A, 1)] .= p0.C
            z = size(p0.A, 1)
            for (k, q, Cpart) in outputs
                m = size(Cpart, 2)
                Cm = zeros(n, nz)
                Cm[:, z + 1:z + m] .= Cpart
                push!(modulations, BlockModulation(k, q, Cm))
                z += m
            end
            # the pump phase turns each harmonic's two outputs into each
            # other, `exp(i k phase) (G_c + i G_s)`, as harmonic balance
            # turns the harmonic (see readharmonics!)
            if !iszero(def.phase)
                for j in 1:2:length(modulations)
                    mc, ms = modulations[j], modulations[j + 1]
                    c, s = cos(mc.harmonic*def.phase), sin(mc.harmonic*def.phase)
                    modulations[j] = BlockModulation(mc.harmonic, 1, c .* mc.C .- s .* ms.C)
                    modulations[j + 1] = BlockModulation(ms.harmonic, 2, s .* mc.C .+ c .* ms.C)
                end
            end
            push!(blocks, TransientBlock(def, snapscattering(p0.D), Float64.(def.zref), cb.signalnodes .- 1, cb.refnodes .- 1,
                offset, cb.path, A, B, C, zbase, modulations, def.wp, def.envelope))
            zbase += nz
            offset += n
            continue
        end
        provider = def.provider
        provider isa TransmissionLineProvider && continue
        if provider isa RationalScatteringProvider
            push!(blocks, TransientBlock(def, snapscattering(provider.D), Float64.(def.zref), cb.signalnodes .- 1, cb.refnodes .- 1,
                offset, cb.path, copy(provider.A), copy(provider.B), copy(provider.C), zbase, BlockModulation[], 0.0, nothing))
            zbase += size(provider.A, 1)
        else
            provider isa ConstantMatrixProvider || throw(ArgumentError(
                lazy"the scattering block at $(cb.path) depends on frequency without a realization in time; give it as RationalScattering or as a TransmissionLine."))
            S = provider.A
            all(x -> isreal(x) && isfinite(x), S) || throw(ArgumentError(
                lazy"the scattering block at $(cb.path) has a complex or nonfinite matrix; a block realized in time is real."))
            push!(blocks, TransientBlock(def, snapscattering(S), Float64.(def.zref), cb.signalnodes .- 1, cb.refnodes .- 1,
                offset, cb.path, zeros(0, 0), zeros(0, def.nports), zeros(def.nports, 0), zbase, BlockModulation[], 0.0, nothing))
        end
        offset += def.nports
    end
    return blocks
end

# the number of states of the rational blocks of a problem
blockstates(p::TransientProblem) = sum(b -> size(b.A, 1), p.blocks; init = 0)

# the inductance scale of a problem: the mean inductance of the circuit,
# or without one `Z0^2 * Cmean`, or `Z0 * 1 ps`
function transientscale(psc::CompiledCircuit, vvn::Vector, nm::CircuitMatrices)
    Lscale = real(calcsolverscale((0.0,), psc.componenttypes, vvn, nm.portimpedances, nm.Lmean))
    isfinite(Lscale) && Lscale > 0 && any(t -> t in (:L, :Lj), psc.componenttypes) && return Lscale
    Z0 = real(calcsolverscale((1.0,), psc.componenttypes, vvn, nm.portimpedances, 1.0))
    caps = [Float64(vvn[k]) for k in eachindex(psc.componenttypes) if psc.componenttypes[k] == :C && vvn[k] > 0]
    return isempty(caps) ? Z0*1e-12 : Z0^2*exp(sum(log, caps)/length(caps))
end

# The directions of the state without inertia, and among them those the
# equations constrain instead of driving. The capacitance matrix is the
# Laplacian of the capacitive graph, so its left null space is spanned by
# the indicators of the subnetworks no capacitor connects to ground: a
# node without a capacitor is such a subnetwork of one node, and so is
# every auxiliary current. Along each of them the equations are
# algebraic, `z' (G v + L x + J(x)) = z' b(t)`, and nothing integrates
# them; a state must satisfy them at the start, and the solver checks
# that it does. Which of them determine a rate and which constrain the
# flux is read off the rate system of those equations: with `Z0` the
# indicators of the capacitor free node subnetworks and `Ea` the block
# current rows, the equations along `Q = [Z0'; Ea']` are linear in the
# rates `alpha` along `Z0` and the block currents `u` through
# `K = [Z0' G Z0  Z0' L Ea; Ea' G Z0  Ea' L Ea]`, the conductances of
# the resistors and the lines and the blocks' hybrid rows, which relate
# a port's current to the rate across it or, where `I + S` is singular,
# constrain the rates. A right null vector of `K` is a flux direction no
# equation's rate determines, an algebraic direction the rule projects;
# the left null vector with it is the combination of the equations, the
# node rows and the block rows, that constrains it, in which the block
# currents cancel; and the rest of `K` determines the other rates and
# the currents, which the endpoint reads. So an open block port adds
# no conductance and leaves a junction's node algebraic, a short or a
# through joins the nodes it ties, and a resistive port grounds one, as
# the equations say rather than as a graph of the ports would guess. A
# coupled inductor current is an algebraic direction of its own, its row
# constraining the flux. Returns the inertialess subnetworks, the
# supports of the algebraic directions, the directions as the columns of
# a sparse matrix, and the constraints as its rows over every equation.
function transientclassification(psc::CompiledCircuit, vvn::Vector, G::SparseMatrixCSC, L::SparseMatrixCSC,
        Nnodal::Int, Naux::Int, blocks = TransientBlock[])
    n = size(G, 1)
    islands = transientsubnetworks(psc, vvn, (:C,))
    inertialess = copy(islands)
    for k in 1:Naux
        push!(inertialess, [Nnodal + k])
    end
    for b in blocks, q in eachindex(b.signal)
        push!(inertialess, [b.auxbase + q])
    end
    k0 = length(islands)
    Z0 = sparse(foldl(append!, islands; init = Int[]), [c for (c, z) in enumerate(islands) for _ in z],
        ones(sum(length, islands; init = 0)), n, k0)
    auxrows = [b.auxbase + q for b in blocks for q in eachindex(b.signal)]
    Ea = sparse(auxrows, 1:length(auxrows), ones(length(auxrows)), n, length(auxrows))
    rs = ratesystem(G, L, Z0, Ea)
    na = length(auxrows)
    # The algebraic directions are the null vectors without a current: a
    # null vector with one, a through's current between two nodes with
    # capacitance, is a current the differential equations determine, of
    # no concern to the projection; a null vector with both a flux and a
    # current part would be a flux direction tied to an undetermined
    # current, which the circuit does not support. A block of the rate
    # system has null vectors of its own (see ratesystem), so its
    # directions are formed within it.
    # the null vectors are orthonormal, so a part below 1e-8 is roundoff
    vr, vc, vv = Int[], Int[], Float64[]
    d = 0
    for b in rs.blocks
        size(b.rightnull, 2) == 0 && continue
        flux, current = findall(<=(k0), b.cols), findall(>(k0), b.cols)
        Nalpha = b.rightnull[flux, :]
        cu = isempty(current) ? Matrix(1.0I, size(Nalpha, 2), size(Nalpha, 2)) :
            nullspace(b.rightnull[current, :]; atol = 1e-8)
        V = orthonormalcolumns(Nalpha*cu)
        rank(Nalpha; atol = 1e-8) == size(V, 2) || throw(ArgumentError(
            "a flux direction without capacitance is tied to a scattering block's port current that no equation determines; the circuit is singular in time."))
        for j in axes(V, 2), (i, c) in enumerate(b.cols[flux])
            push!(vr, c); push!(vc, d + j); push!(vv, V[i, j])
        end
        d += size(V, 2)
    end
    Valpha = sparse(vr, vc, vv, k0, d)
    directions = Z0*Valpha
    # the constraints are the combinations of the equations without any
    # rate: the block currents cancel in every left null vector and the
    # rates along the islands too, and the combinations without the rate
    # of any other node are kept, as many as the directions. The left
    # null vectors over every equation, one row each, a node row taking
    # its island's coefficient; those whose leaks reach no node in common
    # are combined apart
    lr, lc, lv = Int[], Int[], Float64[]
    e = 0
    for b in rs.blocks, j in axes(b.leftnull, 2)
        e += 1
        for (i, r) in enumerate(b.rows)
            x = b.leftnull[i, j]
            for node in (r <= k0 ? islands[r] : (auxrows[r - k0],))
                push!(lr, e); push!(lc, node); push!(lv, x)
            end
        end
    end
    rowst = sparse(lc, lr, lv, n, e)
    leak = sparse(transpose(rowst))*G
    leakt = sparse(transpose(leak))
    atol = 1e-8*max(norm(G, Inf), floatmin(Float64))
    cr, cc, cv = Int[], Int[], Float64[]
    f = 0
    at = zeros(Int, n)
    for (R, C) in blockcomponents(leak)
        isempty(R) && continue
        cl = isempty(C) ? Matrix(1.0I, length(R), length(R)) : nullspace(densesub!(at, leakt, C, R); atol)
        # the group's constraints on the nodes its combinations span
        span = sort!(unique!([rowvals(rowst)[q] for r in R for q in nzrange(rowst, r)]))
        D = densesub!(at, rowst, span, R)*cl
        for j in axes(D, 2), i in axes(D, 1)
            iszero(D[i, j]) && continue
            push!(cr, span[i]); push!(cc, f + j); push!(cv, D[i, j])
        end
        f += size(cl, 2)
    end
    constraints = sparse(cc, cr, cv, f, n)
    size(constraints, 1) == d || throw(ArgumentError(
        "a scattering block ties the rate of a node with capacitance to a constraint on a node without one; that coupling is not supported in time."))
    algebraic = [Int[] for _ in 1:d]
    for j in 1:d, q in nzrange(Valpha, j)
        abs(nonzeros(Valpha)[q]) > 1e-8 && append!(algebraic[j], islands[rowvals(Valpha)[q]])
    end
    foreach(sort!, algebraic)
    # the coupled inductor currents, each its own direction
    coupled = sparse(Nnodal .+ (1:Naux), 1:Naux, ones(Naux), n, Naux)
    directions = hcat(directions, coupled)
    constraints = vcat(constraints, sparse(transpose(coupled)))
    append!(algebraic, [[Nnodal + k] for k in 1:Naux])
    order = sortperm(algebraic; by = z -> (first(z), length(z)))
    return inertialess, algebraic[order], directions[:, order], constraints[order, :]
end

# A block of the rate system, its rows and its columns, which no other
# block shares, and its right and left null vectors as orthonormal columns
struct RateBlock
    rows::Vector{Int}
    cols::Vector{Int}
    rightnull::Matrix{Float64}
    leftnull::Matrix{Float64}
end

# A block of the rate system too large to decompose densely, factorized by
# the sparse QR factorization of its balanced matrix, which sets last the
# columns whose remaining norm falls below the rank tolerance and so
# reveals the block's rank as the singular values reveal a small block's:
# its rows and columns, the factorization's `Q`, its rank and the leading
# triangle of its `R`, its row and column permutations, `Q R` being the
# balanced block permuted, the balancing of the block's rows and columns,
# and the null vectors of the balanced block as orthonormal columns, which
# the minimum norm solution leaves out
struct RateFactor
    rows::Vector{Int}
    cols::Vector{Int}
    Q::SparseArrays.SPQR.QRSparseQ{Float64,Int}
    rank::Int
    R11::SparseMatrixCSC{Float64,Int}
    prow::Vector{Int}
    pcol::Vector{Int}
    dr::Vector{Float64}
    dc::Vector{Float64}
    null::Matrix{Float64}
end

# The pseudoinverse of a rate system on its range, from its rows to its
# columns: the small blocks' as one sparse matrix, and the factorized
# blocks' applied through their factorizations (see ratesolve!)
struct RatePseudoinverse
    dense::SparseMatrixCSC{Float64,Int}
    factors::Vector{RateFactor}
end

# The rate system along the capacitor free islands `Z0` and the block
# current rows `Ea`, balanced by its columns' and rows' magnitudes and
# decomposed: its rank to a relative tolerance, the right and left null
# vectors of the unbalanced system as orthonormal columns, and its
# pseudoinverse on the range, which reads the rates and the currents at
# an endpoint and leaves the null directions alone. The balancing keeps
# a conductance small in the equations' scale from being taken for zero
# next to a block's rows, while an exact dependence, a through's two
# rows, stays one. The islands and the block rows couple only through
# the nodes they share, so the system falls into blocks which share no
# row or column, the connected components of its pattern; each is
# decomposed alone, its singular values held against the largest of them
# all, which is the decomposition of the whole system, and its
# pseudoinverse and null vectors are its own. A block of more than
# `maxdense` rows or columns, a cascade of blocks through nodes without
# capacitance say, is factorized by sparse QR instead (see RateFactor),
# whose cost grows with its pattern rather than as the cube of its size,
# at the same tolerance, the largest singular value bounded there by the
# norms of its balanced matrix. Returns the blocks, the pseudoinverse and
# the size of the system.
function ratesystem(G::SparseMatrixCSC, L::SparseMatrixCSC, Z0::SparseMatrixCSC, Ea::SparseMatrixCSC;
        rtol = 1e-8, maxdense::Integer = 64)
    Q = sparse(transpose(hcat(Z0, Ea)))
    K = dropzeros!(hcat(Q*G*Z0, Q*L*Ea))
    m = size(K, 1)
    dc = ones(m)
    for j in 1:m
        x = maximum(abs, view(nonzeros(K), nzrange(K, j)); init = 0.0)
        x > 0 && (dc[j] = 1/x)
    end
    Kb = K*Diagonal(dc)
    rmax = zeros(m)
    for (i, x) in zip(rowvals(Kb), nonzeros(Kb))
        rmax[i] = max(rmax[i], abs(x))
    end
    dr = [x > 0 ? 1/x : 1.0 for x in rmax]
    Kb = Diagonal(dr)*Kb
    parts = blockcomponents(Kb)
    large = ((R, C),) -> max(length(R), length(C)) > maxdense
    at = zeros(Int, m)
    F = [isempty(R) || isempty(C) || large((R, C)) ? nothing : svd(densesub!(at, Kb, R, C); full = true) for (R, C) in parts]
    smax = maximum((first(f.S) for f in F if !isnothing(f)); init = 0.0)
    for (R, C) in Iterators.filter(large, parts)
        B = Kb[R, C]
        smax = max(smax, sqrt(opnorm(B, 1)*opnorm(B, Inf)))
    end
    blocks, factors = RateBlock[], RateFactor[]
    pr, pc, pv = Int[], Int[], Float64[]
    for ((R, C), f) in zip(parts, F)
        if large((R, C))
            q = qr(Kb[R, C]; tol = rtol*smax)
            r = rank(q)
            R11 = q.R[1:r, 1:r]
            prow, pcol = q.prow, q.pcol
            # the null vectors: each dependent column a unit step, the
            # independent ones solved for through the triangle
            X = Matrix(q.R[1:r, r + 1:end])
            ldiv!(UpperTriangular(R11), X)
            N = zeros(length(C), length(C) - r)
            N[pcol, :] = vcat(-X, Matrix(1.0I, length(C) - r, length(C) - r))
            Lb = zeros(length(R), length(R) - r)
            Lb[prow, :] = q.Q*vcat(zeros(r, length(R) - r), Matrix(1.0I, length(R) - r, length(R) - r))
            push!(blocks, RateBlock(R, C, orthonormalcolumns(dc[C] .* N), orthonormalcolumns(dr[R] .* Lb)))
            push!(factors, RateFactor(R, C, q.Q, r, R11, prow, pcol, dr[R], dc[C], orthonormalcolumns(N)))
        elseif isnothing(f)
            # a row or a column without an entry, its own null vector
            push!(blocks, RateBlock(R, C, Matrix(1.0I, length(C), length(C)), Matrix(1.0I, length(R), length(R))))
        else
            r = count(s -> s > rtol*smax, f.S)
            pinv = (dc[C] .* f.V[:, 1:r])*Diagonal(1 ./ f.S[1:r])*transpose(f.U[:, 1:r] .* dr[R])
            for (jj, row) in enumerate(R), (ii, col) in enumerate(C)
                push!(pr, col); push!(pc, row); push!(pv, pinv[ii, jj])
            end
            push!(blocks, RateBlock(R, C, orthonormalcolumns(dc[C] .* f.V[:, r + 1:end]),
                orthonormalcolumns(dr[R] .* f.U[:, r + 1:end])))
        end
    end
    return (; blocks, pinv = RatePseudoinverse(sparse(pr, pc, pv, m, m), factors), m)
end

# `theta = P g`, or with `transposed` `theta = P' g`, for the
# pseudoinverse `P` of a rate system, over the columns of `g`: the small
# blocks through their sparse matrix and each factorized block through its
# factorization, with `work` the buffers of `ratework`
function ratesolve!(theta, P::RatePseudoinverse, g, work; transposed::Bool = false)
    transposed ? mul!(theta, transpose(P.dense), g) : mul!(theta, P.dense, g)
    for (f, w) in zip(P.factors, work)
        transposed ? factorsolvetranspose!(theta, f, g, w) : factorsolve!(theta, f, g, w)
    end
    return theta
end

# the buffers of the factorized blocks' solves over `m` columns: for each
# block a column over its rows, a column over its columns, and the
# coefficients along its null vectors
ratework(P::RatePseudoinverse, m::Integer) =
    [(zeros(length(f.rows), m), zeros(length(f.cols), m), zeros(size(f.null, 2), m)) for f in P.factors]

# The minimum norm solution through a factorized block,
# `theta[cols] = dc .* (B^+ (dr .* g[rows]))` for its balanced matrix `B`:
# the right hand side along `Q'`, the triangle solved for the independent
# columns, the dependent ones zero, and the null vectors taken out
function factorsolve!(theta, f::RateFactor, g, (b, x, c))
    r = f.rank
    for (k, i) in enumerate(f.prow), j in axes(g, 2)
        b[k, j] = f.dr[i]*g[f.rows[i], j]
    end
    lmul!(adjoint(f.Q), b)
    ldiv!(UpperTriangular(f.R11), view(b, 1:r, :))
    fill!(x, 0)
    for q in 1:r, j in axes(g, 2)
        x[f.pcol[q], j] = b[q, j]
    end
    mul!(c, transpose(f.null), x)
    mul!(x, f.null, c, -1.0, 1.0)
    for (k, col) in enumerate(f.cols), j in axes(g, 2)
        theta[col, j] = f.dc[k]*x[k, j]
    end
    return theta
end

# Its transpose, `theta[rows] = dr .* (B^+' (dc .* g[cols]))`: the right
# hand side with the null vectors taken out, the transposed triangle solved
# along the independent columns, and the solution along `Q`
function factorsolvetranspose!(theta, f::RateFactor, g, (b, x, c))
    r = f.rank
    for (k, col) in enumerate(f.cols), j in axes(g, 2)
        x[k, j] = f.dc[k]*g[col, j]
    end
    mul!(c, transpose(f.null), x)
    mul!(x, f.null, c, -1.0, 1.0)
    fill!(b, 0)
    for q in 1:r, j in axes(g, 2)
        b[q, j] = x[f.pcol[q], j]
    end
    ldiv!(transpose(UpperTriangular(f.R11)), view(b, 1:r, :))
    lmul!(f.Q, b)
    for (k, i) in enumerate(f.prow), j in axes(g, 2)
        theta[f.rows[i], j] = f.dr[i]*b[k, j]
    end
    return theta
end

# The dense submatrix of `A` on the rows `R` and the columns `C`, read from
# the entries of its columns through `at`, a vector over its rows which
# holds zeros and is left so: work in proportion to the entries read,
# whatever the size of `A`.
function densesub!(at::Vector{Int}, A::SparseMatrixCSC, R, C)
    for (i, r) in enumerate(R)
        at[r] = i
    end
    B = zeros(length(R), length(C))
    for (j, c) in enumerate(C), q in nzrange(A, c)
        i = at[rowvals(A)[q]]
        i > 0 && (B[i, j] = nonzeros(A)[q])
    end
    for r in R
        at[r] = 0
    end
    return B
end

# The connected components of the pattern of `A`, its rows and its
# columns joined by its entries: blocks which share no row or column,
# each as its sorted rows and columns, in the order of their first
# column, and a row without an entry a block of its own after them.
function blockcomponents(A::SparseMatrixCSC)
    m, n = size(A)
    At = sparse(transpose(A))
    rowblock, colblock = falses(m), falses(n)
    blocks = Tuple{Vector{Int},Vector{Int}}[]
    queue = Int[]
    for j0 in 1:n
        colblock[j0] && continue
        rows, cols = Int[], Int[]
        colblock[j0] = true
        push!(queue, j0)
        while !isempty(queue)
            j = pop!(queue)
            push!(cols, j)
            for q in nzrange(A, j)
                i = rowvals(A)[q]
                rowblock[i] && continue
                rowblock[i] = true
                push!(rows, i)
                for p in nzrange(At, i)
                    k = rowvals(At)[p]
                    colblock[k] && continue
                    colblock[k] = true
                    push!(queue, k)
                end
            end
        end
        push!(blocks, (sort!(rows), sort!(cols)))
    end
    for i in 1:m
        rowblock[i] || push!(blocks, ([i], Int[]))
    end
    return blocks
end

# an orthonormal basis of the span of the columns of `M`
function orthonormalcolumns(M::AbstractMatrix; rtol = 1e-10)
    size(M, 2) == 0 && return zeros(size(M, 1), 0)
    F = svd(M)
    r = count(s -> s > rtol*max(F.S[1], floatmin(Float64)), F.S)
    return F.U[:, 1:r]
end

# the branches of the elements of the given types which have a finite
# nonzero value, and with `blocks` the ports of every scattering block or
# line, which connect their two terminals as an element does
function transientedges(psc::CompiledCircuit, vvn::Vector, types;
        blocks::Bool = false)
    edges = Tuple{Int,Int}[]
    for k in eachindex(psc.componenttypes)
        psc.componenttypes[k] in types || continue
        v = vvn[k]
        (v isa Number && isfinite(v) && !iszero(v)) || continue
        push!(edges, (psc.nodeindices[1, k], psc.nodeindices[2, k]))
    end
    blocks && for cb in psc.scatteringblocks, q in eachindex(cb.signalnodes)
        push!(edges, (cb.signalnodes[q], cb.refnodes[q]))
    end
    return edges
end

# the subnetworks of the nodes that no element of the given types with a
# finite nonzero value connects to ground, as sorted lists of state indices
# in the order of their first node
function transientsubnetworks(psc::CompiledCircuit, vvn::Vector, types)
    nodes = nodecomponents(psc.Nnodes, transientedges(psc, vvn, types))
    return [n .- 1 for n in nodes]
end

# The floating subnetworks of the circuit in time. Harmonic balance gauge
# fixes the components of the static flux stiffness graph, whose edges are
# the inductors and junctions, because at zero frequency a capacitor and a
# resistor are open and the flux of a node reached through them alone is
# undetermined. In time that flux is the integral of the node's voltage and
# is determined; only a subnetwork no element connects to ground has a free
# flux offset, and it is those the gauge rows fix. The same union-find as
# `calcstaticfluxcomponents`, over every two terminal element with a finite
# nonzero value.
transientfloatingcomponents(psc::CompiledCircuit, vvn::Vector) =
    nodecomponents(psc.Nnodes,
        transientedges(psc, vvn, (:C, :R, :L, :Lj); blocks = true))

"""
    TransientState

The state a transient starts from or ends at, in the solver's units: the
scaled node fluxes `flux/phi0` with the auxiliary currents of the coupled
inductors and the scattering blocks appended in the problem's fixed units
`Lscale*i/phi0`, their rates, the history of the wave leaving each port
of each transmission line in sqrt(W), two per line in compiled order, as
the columns of `waves` at the spacing `wavesdt` ending at the start, or a
single column the lines carry unchanged before the start, and the states
of the rational blocks. Built by [`transientstate`](@ref), from the
physical node fluxes and voltages of a problem or from the end of a
solution, and given to [`transientsolve`](@ref) as `initialstate`, which
takes nothing else, since a bare pair of arrays does not say which units
it is in. A state does not depend on the step: a history recorded at one
step is read at another through the same interpolation the lines read
their history with.
"""
struct TransientState
    # host vectors of `Float64`, to which the constructor converts the
    # vectors it is given, a device's among them: the one form the step
    # drivers take, so that they compile once whatever form a state is
    # given in
    flux::Vector{Float64}
    rate::Vector{Float64}
    waves::Matrix{Float64}
    wavesdt::Float64
    blockstates::Vector{Float64}
end

Base.:(==)(a::TransientState, b::TransientState) = a.flux == b.flux && a.rate == b.rate && a.waves == b.waves &&
    a.wavesdt == b.wavesdt && a.blockstates == b.blockstates

"""
    transientstate(problem; flux = zeros(...), voltage = zeros(...),
        linecurrents = zeros(...))
    transientstate(solution)

The initial state of a transient, a [`TransientState`](@ref), from the
node fluxes in Weber and the node voltages in Volts, in the compiled node
order without ground: the scaled fluxes `flux/phi0`, augmented with the
auxiliary currents the coupled inductors' constitutive equations imply
and normalized into the gauge of the floating subnetworks, and the scaled
flux rates `voltage/phi0` with the auxiliary rates the same equations
imply. The default is the zero state. [`transientsolve`](@ref) checks
that a state satisfies the algebraic equations of the circuit at the
start; it does not project one that does not. With transmission lines
the waves leaving the ports of each line come from the port voltages and
`linecurrents`, the direct current into the first port of each line,
zero by default; before the start the lines carry those waves unchanged.

From a solution, the state at its end, to start another solve from: its
final fluxes and rates, the waves leaving each line port over the delay
window before the end, which is what the lines read after the start, so
a continuation is the uninterrupted solve to the solver's tolerance at
the same step and to the interpolation of the history at another, and
the final states of the rational blocks, which every solve of the
package's rules keeps whatever it records.
"""
function transientstate(p::TransientProblem; flux = zeros(p.Nnodal), voltage = zeros(p.Nnodal),
        linecurrents = zeros(length(p.lines)))
    length(flux) == p.Nnodal && length(voltage) == p.Nnodal || throw(DimensionMismatch(
        lazy"the circuit has $(p.Nnodal) nodes without ground; give one flux and one voltage per node."))
    all(isfinite, flux) && all(isfinite, voltage) || throw(ArgumentError("the initial state must be finite."))
    x = vcat(Float64.(flux) ./ phi0, zeros(p.Naux))
    v = vcat(Float64.(voltage) ./ phi0, zeros(p.Naux))
    # the gauge and the auxiliary currents, as hbnlsolve normalizes an
    # initial guess; the constitutive rows are algebraic, so the auxiliary
    # rates follow the node rates the same way
    xc = complex(x)
    mnagaugenormalize!(xc, p.floatingcomponents, [0.0], 1)
    mnainitialauxind!(xc, p.coupledbranches, p.matrices.Lb, p.matrices.Mb,
        p.circuit.topology.Rbn, 1, p.Nnodal, p.Lscale)
    vc = complex(v)
    mnagaugenormalize!(vc, p.floatingcomponents, [0.0], 1)
    mnainitialauxind!(vc, p.coupledbranches, p.matrices.Lb, p.matrices.Mb,
        p.circuit.topology.Rbn, 1, p.Nnodal, p.Lscale)
    xr = real.(xc)
    # the port currents of the blocks from their constitutive rows at the
    # port voltages, where the rows determine them
    for b in p.blocks
        vE = [(b.signal[q] > 0 ? voltage[b.signal[q]] : 0.0) - (b.ref[q] > 0 ? voltage[b.ref[q]] : 0.0) for q in eachindex(b.signal)]
        Cb = (I + b.S) .* transpose(sqrt.(b.R))
        rank(Cb) == size(Cb, 1) || continue
        Bb = (I - b.S) .* transpose(1 ./ sqrt.(b.R))
        xr[b.auxbase + 1:b.auxbase + length(b.signal)] .= Cb \ (p.Lscale .* (Bb*vE) ./ phi0)
    end
    # the states of the rational blocks at rest under the port voltages,
    # `z = -A^(-1) B a` with the incident waves `a` of the direct current
    # solution of each block's own equations at zero frequency
    zs = zeros(blockstates(p))
    for b in p.blocks
        nz = size(b.A, 1)
        nz == 0 && continue
        vE = [(b.signal[q] > 0 ? voltage[b.signal[q]] : 0.0) - (b.ref[q] > 0 ? voltage[b.ref[q]] : 0.0) for q in eachindex(b.signal)]
        S0 = b.S .- b.C*(b.A \ b.B)
        Cb = (I + S0) .* transpose(sqrt.(b.R))
        rank(Cb) == size(Cb, 1) || continue
        i0 = Cb \ (((I - S0) .* transpose(1 ./ sqrt.(b.R)))*vE)
        a0 = (vE ./ sqrt.(b.R) .+ sqrt.(b.R) .* i0) ./ 2
        zs[b.zbase + 1:b.zbase + nz] .= -(b.A \ (b.B*a0))
        xr[b.auxbase + 1:b.auxbase + length(b.signal)] .= p.Lscale .* i0 ./ phi0
    end
    length(linecurrents) == length(p.lines) || throw(DimensionMismatch(lazy"give one direct current per transmission line ($(length(p.lines)))."))
    waves = zeros(2length(p.lines))
    for (l, line) in enumerate(p.lines)
        vp = [(line.signal[q] > 0 ? voltage[line.signal[q]] : 0.0) - (line.ref[q] > 0 ? voltage[line.ref[q]] : 0.0) for q in 1:2]
        waves[2l - 1] = (vp[1]/sqrt(line.Z) + sqrt(line.Z)*linecurrents[l])/2
        waves[2l] = (vp[2]/sqrt(line.Z) - sqrt(line.Z)*linecurrents[l])/2
    end
    return TransientState(xr, real.(vc), reshape(waves, :, 1), 0.0, zs)
end

# The drive of a harmonic balance solution in time, one source for each
# of its sources: at a nonzero mode of the frequency `f` the current
# `2 real(current*cis(f t))`, and at the zero mode the current itself.
# Harmonic balance drives a port along its branch, into the branch's
# destination whichever order the port's terminals are written in, and
# the transient into the port's positive terminal, so the current is
# turned where the two differ.
function orbitsources(psc::CompiledCircuit, nonlinear::NonlinearHB)
    w = collect(Float64, nonlinear.w)
    _, destination = branchendpoints(psc.topology.Rbn, psc.topology.Nbranches)
    return map(nonlinear.sources) do s
        port = psc.ports[findfirst(q -> q.number == s.port, psc.ports)]
        b = psc.topology.edge2indexdict[(port.positivenode, port.negativenode)]
        current = destination[b] == port.positivenode ? s.current : -s.current
        all(iszero, s.mode) && return TransientSource(s.port, real(current))
        f = sum(s.mode .* w)
        return TransientSource(s.port, t -> 2real(current*cis(f*t)))
    end
end

# The problem with the conversion of every pumped block always on: the
# periodic device harmonic balance and the pole analysis describe, which
# the envelope of a block gates only for a record that starts it.
function alwayson(p::TransientProblem)
    all(b -> isnothing(b.envelope), p.blocks) && return p
    blocks = [TransientBlock(b.definition, b.S, b.R, b.signal, b.ref, b.auxbase, b.path, b.A, b.B, b.C, b.zbase,
        b.modulations, b.wp, nothing) for b in p.blocks]
    return TransientProblem(p; blocks)
end

# The state of a transient on a harmonic balance orbit at time zero, the
# origin of its phases: the node fluxes and voltages of the orbit's
# Fourier series, with the direct voltage the solution holds apart, the
# port currents of the blocks as the solution determined them with the
# whole circuit, the states of each block's filters in the steady state
# of every mode, under the incident waves the port voltages and currents
# make, and the history of the wave leaving each line port,
# `(v + Z i)/(2 sqrt(Z))` from the port's voltage and the current into the
# line the solution determined, over the prehistory a solve at the step
# `dt` reads, at that step. A pumped block's filters are driven by its
# incident waves alone, its modulations acting on their outputs, so the
# states of every mode are formed as an unpumped block's are.
function orbitstate(p::TransientProblem, nonlinear::NonlinearHB; dt::Union{Nothing,Real} = nothing)
    isempty(p.lines) || !isnothing(dt) || throw(ArgumentError(
        "a transient on a harmonic balance orbit through a transmission line starts from the line's history, sampled at the step: give dt."))
    modes = nonlinear.modes
    F = reshape(initialguess(nonlinear.nodeflux), length(modes), :)
    size(F, 2) == p.Nnodal || throw(DimensionMismatch(
        "the harmonic balance solution is of a circuit with another number of nodes."))
    w = collect(Float64, nonlinear.w)
    f = [sum(mode .* w) for mode in modes]
    # the real signal: the zero mode once, every other mode twice its
    # real part
    c = [all(iszero, mode) ? 1.0 : 2.0 for mode in modes]
    dc = isnothing(nonlinear.dcnodevoltage) ? zeros(p.Nnodal) : real.(initialguess(nonlinear.dcnodevoltage))
    flux = phi0 .* vec(sum(real.(c .* F); dims = 1))
    voltage = dc .+ phi0 .* vec(sum(real.(c .* (im .* f) .* F); dims = 1))
    state = transientstate(p; flux, voltage)
    x, v, zs = copy(state.flux), copy(state.rate), copy(state.blockstates)
    modevoltage(k) = all(iszero, modes[k]) ? complex(dc) : phi0 .* (im*f[k]) .* F[k, :]
    block = Dict(cb.path => j for (j, cb) in enumerate(p.circuit.scatteringblocks))
    for b in p.blocks
        n, nz = length(b.signal), size(b.A, 1)
        currents = nonlinear.blockcurrents[block[b.path]]
        current, rate, states = zeros(n), zeros(n), zeros(nz)
        for k in eachindex(modes)
            V = modevoltage(k)
            vk = [(b.signal[q] > 0 ? V[b.signal[q]] : zero(eltype(V))) - (b.ref[q] > 0 ? V[b.ref[q]] : zero(eltype(V))) for q in 1:n]
            ik = currents[k, :]
            all(iszero, vk) && all(iszero, ik) && continue
            current .+= real.(c[k] .* ik)
            rate .+= real.(c[k] .* (im*f[k]) .* ik)
            ak = (vk ./ sqrt.(b.R) .+ sqrt.(b.R) .* ik) ./ 2
            nz > 0 && (states .+= real.(c[k] .* ((im*f[k]*I - b.A) \ (b.B*ak))))
        end
        x[b.auxbase + 1:b.auxbase + n] .= p.Lscale .* current ./ phi0
        v[b.auxbase + 1:b.auxbase + n] .= p.Lscale .* rate ./ phi0
        zs[b.zbase + 1:b.zbase + nz] .= states
    end
    isempty(p.lines) && return TransientState(x, v, state.waves, state.wavesdt, zs)
    # the columns end at zero, as a solve's history ends at its start
    npre = lineprehistory(p, dt)
    times = -(npre - 1:-1:0) .* dt
    waves = zeros(2length(p.lines), npre)
    for k in eachindex(modes)
        V = modevoltage(k)
        weight = all(iszero, modes[k]) ? 1.0 : 2.0
        for (l, line) in enumerate(p.lines), e in 1:2
            vk = (line.signal[e] > 0 ? V[line.signal[e]] : zero(eltype(V))) - (line.ref[e] > 0 ? V[line.ref[e]] : zero(eltype(V)))
            ak = (vk + line.Z*nonlinear.blockcurrents[block[line.path]][k, e])/(2sqrt(line.Z))
            iszero(ak) || (view(waves, 2(l - 1) + e, :) .+= weight .* real.(cis.(f[k] .* times) .* ak))
        end
    end
    return TransientState(x, v, waves, Float64(dt), zs)
end

# the initial states of the rational blocks of a state, zero when the
# state has none
initialblockstates(state::TransientState, p::TransientProblem) = isempty(state.blockstates) ? zeros(blockstates(p)) : Float64.(collect(state.blockstates))

# The history of the waves leaving the line ports before the start, the
# `npre` columns at the step `h` a solve reads, from a state: a single
# column carried unchanged, a history at the same step as it is, its
# earliest column repeated before it begins, and a history at another
# step read at the columns' times through the same stencil the lines
# read their history with. Zero for a circuit without lines.
function initialwaves(state::TransientState, p::TransientProblem, h, npre::Int)
    nl = 2length(p.lines)
    out = zeros(nl, npre)
    (nl == 0 || isempty(state.waves)) && return out
    w = state.waves
    size(w, 1) == nl || throw(DimensionMismatch(lazy"the state needs $(nl) line waves; use transientstate."))
    nh = size(w, 2)
    if nh == 1 || state.wavesdt <= 0
        out .= view(w, :, nh)
    elseif state.wavesdt == h
        for c in 1:npre
            out[:, c] .= view(w, :, max(nh - (npre - c), 1))
        end
    else
        weights = zeros(6)
        tpre = -(nh - 1)*state.wavesdt
        for c in 1:npre
            first, nst = linestencil!(weights, -(npre - c)*h, tpre, state.wavesdt, nh)
            for k in 1:nst
                out[:, c] .+= weights[k] .* view(w, :, first + k)
            end
        end
    end
    return out
end

# the current each line port forces into its node from the far port's
# wave a delay earlier, `2 q / sqrt(Z)` in Amperes, for the waves `q`
# arriving at the ports
lineforcing(p::TransientProblem, q) = [2q[k]/sqrt(p.lines[(k + 1) ÷ 2].Z) for k in eachindex(q)]
# the waves arriving at the ports before the start: the far ports'
# initial waves, which the lines carry unchanged
arrivingwaves(p::TransientProblem, waves) = [waves[isodd(k) ? k + 1 : k - 1] for k in eachindex(waves)]

# The waves arriving at the line ports at the start of a solve, at the
# time `t` after it, from a state's own history at its own step: the far
# port's wave a delay earlier, read with the same stencil the lines read
# theirs with, and the far port's single wave for a constant history.
# The check of the algebraic equations reads these, so a state is judged
# against the history it carries; the interpolation which resamples that
# history onto another step is an error of the continuation and not an
# inconsistency of the state.
function initialarrivals(state::TransientState, p::TransientProblem, t = 0.0)
    nl = 2length(p.lines)
    q = zeros(nl)
    (nl == 0 || isempty(state.waves)) && return q
    w = state.waves
    nh = size(w, 2)
    (nh == 1 || state.wavesdt <= 0) && return arrivingwaves(p, view(w, :, nh))
    weights = zeros(6)
    tpre = -(nh - 1)*state.wavesdt
    for (l, line) in enumerate(p.lines), e in 1:2
        row = 2(l - 1) + e
        far = isodd(row) ? row + 1 : row - 1
        first, nst = linestencil!(weights, t - line.delay, tpre, state.wavesdt, nh)
        for k in 1:nst
            q[row] += weights[k]*w[far, first + k]
        end
    end
    return q
end

# The currents the lines force at the start of a solve from a state's own
# history, `2 q / sqrt(Z)` of the arriving waves, and their rate by the
# central difference of `delta` a stepper's reading takes of its history
# (see `readstepper!`), so that a state taken from the end of a solve
# meets the check of the next start as it met the reading there
function stateforcing(state::TransientState, p::TransientProblem, delta)
    forced = lineforcing(p, initialarrivals(state, p))
    rate = (lineforcing(p, initialarrivals(state, p, delta)) .- lineforcing(p, initialarrivals(state, p, -delta))) ./ (2delta)
    return forced, rate
end
