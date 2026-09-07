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
transientwaveform(w::TransientWaveform) = w

"""
    TransientSource(target, current)

A real instantaneous current in Amperes, a number or a callable
`current(t)` of the time in seconds. An integer `target` is a port
number: positive current is injected into the port's positive terminal by
a Norton source, and the port's termination is part of the compiled
circuit already. A string or symbol names a `CurrentSource` component,
whose constant value the waveform replaces; positive current then flows
out of the component's first terminal and into its second. Several
sources on one target add. Unlike a harmonic balance source the callable
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
an open is stamped exactly. A block whose matrix depends on frequency
needs a causal realization, which the transient does not yet have.
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
history with a cubic, centered when the delay allows.
"""
struct TransientLine
    Z::Float64
    delay::Float64
    signal::NTuple{2,Int}
    ref::NTuple{2,Int}
    path::String
end

"""
    TransientProblem

A compiled circuit with its matrices at one mode, its modified nodal
analysis augmentation, its bound drives and its port data, ready to be
integrated in time by [`transientsolve`](@ref). Built by
[`transientproblem`](@ref).

# Fields
- `circuit`, `graph`, `matrices`: the compiled circuit, its graph and its
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
- `injection`: the unscaled injection of a unit current of each drive
    into the node equations, one sparse column per drive, in the
    orientation of the drive.
- `drives`: the bound drives, in the order of `injection`'s columns.
- `constantcurrent`: the unscaled constant node current of the netlist's
    current sources not replaced by a drive.
- `portpositive`, `portnegative`, `portimpedances`, `portconductances`:
    the node of each port terminal (zero for ground), the reference
    impedance and the conductance of the termination the port owns, in
    compiled port order, for the port waves.
- `blocks`: the scattering parameter blocks with a realization in time,
    see [`TransientBlock`](@ref), whose port currents are auxiliary
    unknowns after the coupled inductor currents.
- `lines`: the ideal transmission lines, see [`TransientLine`](@ref),
    whose wave histories the solve keeps.
"""
struct TransientProblem
    circuit::CompiledCircuit
    graph::CircuitGraph
    matrices::CircuitMatrices
    Nnodal::Int
    Naux::Int
    Lscale::Float64
    coupledbranches::Vector{Int}
    floatingcomponents::Vector{Vector{Int}}
    gaugeindices::Vector{Int}
    inertialess::Vector{Vector{Int}}
    algebraic::Vector{Vector{Int}}
    directions::Matrix{Float64}
    constraints::Matrix{Float64}
    rateextraction::Matrix{Float64}
    injection::SparseMatrixCSC{Float64,Int}
    drives::Vector{TransientDrive}
    constantcurrent::Vector{Float64}
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
    relations::Union{Nothing,JunctionRelations{Matrix{Float64},Vector{Bool}}}
end

Base.length(p::TransientProblem) = p.Nnodal + p.Naux

# a real, finite value of a component, or an argument error naming it
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
    transientproblem(circuit, circuitdefs = Dict(); sources = (),
        sorting = defaultsorting(circuit))

Compile a circuit for integration in time: the same compiler and
[`numericmatrices`](@ref) as harmonic balance, at one mode, with the
mutually coupled inductors promoted to auxiliary branch currents and the
floating inductive subnetworks gauge fixed as [`hbnlsolve`](@ref) does.
`sources` is a tuple or vector of [`TransientSource`](@ref)s; the
netlist's constant `CurrentSource` components keep their constant values
unless a source names them. The circuit may be a typed [`Circuit`](@ref),
a compiled circuit or a legacy netlist.

Supported are real, constant resistors, capacitors, inductors, mutual
inductors, sinusoidal Josephson junctions, current sources and ports.
Frequency dependent or complex values are rejected, since they need a
causal realization in time. A [`ScatteringParameters`](@ref) block with a
constant real matrix is realized as it is, see [`TransientBlock`](@ref);
any other block is rejected for the same reason.
"""
function transientproblem(circuit, circuitdefs = Dict{Symbol,Any}();
        sources = (), sorting::Symbol = defaultsorting(circuit))
    psc = compile(circuit; sorting)
    cg = calccircuitgraph(psc; loops = false)
    vvn = componentvaluestonumber(psc.componentvalues, circuitdefs)
    checkcomponentvaluesdefined(psc.componentnames, vvn, circuitdefs)
    for k in eachindex(vvn)
        v = transientreal(vvn[k], psc.componentnames[k])
        psc.componenttypes[k] == :R || isfinite(v) || throw(ArgumentError(
            lazy"the component $(psc.componentnames[k]) has an infinite value; only a resistor may be an open."))
        vvn[k] = v
    end
    checkstaticstiffnessvalues(psc.componenttypes, vvn)
    nm = numericmatrices(psc, cg, vvn; Nmodes = 1)
    for (name, A) in (("capacitance", nm.Cnm), ("conductance", nm.Gnm),
            ("inverse inductance", nm.invLnm))
        all(isfinite, nonzeros(A)) || throw(ArgumentError(
            lazy"the $(name) matrix has a nonfinite entry."))
    end
    all(>(0), nm.Cnm.nzval[[k for j in axes(nm.Cnm, 2) for k in nzrange(nm.Cnm, j) if rowvals(nm.Cnm)[k] == j]]) ||
        throw(ArgumentError("every node needs a nonnegative capacitance to ground; a negative capacitance has no meaning in time."))

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
    _, Gs, Ls, _ = transientlinearmatrices(nm, coupledbranches, cg.Rbn, gaugeindices, blocks, lines, Lscale, Nnodal, Naux)
    inertialess, algebraic, directions, constraints, rateextraction =
        transientclassification(psc, vvn, Gs, Ls, Nnodal, length(coupledbranches), blocks)

    # the port terminals and terminations, for the port waves
    np = length(psc.ports)
    portpositive = [p.positivenode - 1 for p in psc.ports]
    portnegative = [p.negativenode - 1 for p in psc.ports]
    portimpedances = [transientreal(nm.portimpedances[k], "the impedance of port $(psc.ports[k].number)") for k in 1:np]
    all(>(0), portimpedances) || throw(ArgumentError("port reference impedances must be positive."))
    portconductances = [p.environment == 0 ? 0.0 : 1/vvn[p.environment] for p in psc.ports]

    # the drives: a unit current of each source as a node injection, one
    # column per source, and the constant current of the netlist sources
    # no source replaced
    drives = TransientDrive[]
    rows, cols, vals = Int[], Int[], Float64[]
    replaced = Set{Int}()
    for (k, source) in enumerate(sources)
        source isa TransientSource || throw(ArgumentError("sources must contain TransientSource objects."))
        if source.target isa Int
            p = findfirst(port -> port.number == source.target, psc.ports)
            isnothing(p) && throw(ArgumentError(lazy"there is no port $(source.target)."))
            n1, n2 = portpositive[p], portnegative[p]
            push!(drives, TransientDrive(p, source.current))
        else
            c = get(psc.componentnamedict, source.target, 0)
            (c > 0 && psc.componenttypes[c] == :I) || throw(ArgumentError(
                lazy"$(source.target) does not name a CurrentSource of the circuit."))
            push!(replaced, c)
            # out of the first terminal, into the second
            n1, n2 = psc.nodeindices[2, c] - 1, psc.nodeindices[1, c] - 1
            push!(drives, TransientDrive(0, source.current))
        end
        n1 > 0 && (push!(rows, n1); push!(cols, k); push!(vals, 1.0))
        n2 > 0 && (push!(rows, n2); push!(cols, k); push!(vals, -1.0))
    end
    injection = sparse(rows, cols, vals, Nnodal + Naux, length(drives))
    constantcurrent = zeros(Nnodal + Naux)
    for c in psc.currentsources
        c in replaced && continue
        n1, n2 = psc.nodeindices[2, c] - 1, psc.nodeindices[1, c] - 1
        n1 > 0 && (constantcurrent[n1] += vvn[c])
        n2 > 0 && (constantcurrent[n2] -= vvn[c])
    end
    return TransientProblem(psc, cg, nm, Nnodal, Naux, Lscale, coupledbranches,
        floatingcomponents, gaugeindices, inertialess, algebraic, directions, constraints, rateextraction,
        injection, drives, constantcurrent, portpositive, portnegative, portimpedances, portconductances, blocks, lines,
        calcjunctionrelations(psc.componenttypes, psc.nodeindices,
            psc.junctioncprs, cg.edge2indexdict, nm.Ljb))
end

# the ideal lines of a compiled circuit, in compiled order
function transientlines(psc::CompiledCircuit)
    lines = TransientLine[]
    for cb in psc.scatteringblocks
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


# the blocks of a compiled circuit the transient realizes, in compiled
# order, their port currents laid out from `offset`
function transientblocks(psc::CompiledCircuit, offset::Int)
    blocks = TransientBlock[]
    zbase = 0
    for cb in psc.scatteringblocks
        def = cb.definition
        provider = def.provider
        provider isa TransmissionLineProvider && continue
        if provider isa RationalScatteringProvider
            push!(blocks, TransientBlock(def, copy(provider.D), Float64.(def.zref), cb.signalnodes .- 1, cb.refnodes .- 1,
                offset, cb.path, copy(provider.A), copy(provider.B), copy(provider.C), zbase))
            zbase += size(provider.A, 1)
        else
            provider isa ConstantMatrixProvider || throw(ArgumentError(
                lazy"the scattering block at $(cb.path) depends on frequency without a realization in time; give it as RationalScattering or as a TransmissionLine."))
            S = provider.A
            all(x -> isreal(x) && isfinite(x), S) || throw(ArgumentError(
                lazy"the scattering block at $(cb.path) has a complex or nonfinite matrix; a block realized in time is real."))
            push!(blocks, TransientBlock(def, Matrix{Float64}(real.(S)), Float64.(def.zref), cb.signalnodes .- 1, cb.refnodes .- 1,
                offset, cb.path, zeros(0, 0), zeros(0, def.nports), zeros(def.nports, 0), zbase))
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
# supports of the algebraic directions, the directions as columns, the
# constraints as rows over every equation, and the rows that extract the
# rate along each direction from a rate of the state.
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
    Z0 = sparse(reduce(vcat, islands; init = Int[]), reduce(vcat, [fill(c, length(z)) for (c, z) in enumerate(islands)]; init = Int[]),
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
    # current, which the circuit does not support.
    # the null vectors are orthonormal, so a part below 1e-8 is roundoff
    Nalpha, Nu = rs.rightnull[1:k0, :], rs.rightnull[k0 + 1:k0 + na, :]
    cu = na == 0 ? Matrix(1.0I, size(Nalpha, 2), size(Nalpha, 2)) : nullspace(Nu; atol = 1e-8)
    Valpha = orthonormalcolumns(Nalpha*cu)
    d = size(Valpha, 2)
    rank(Nalpha; atol = 1e-8) == d || throw(ArgumentError(
        "a flux direction without capacitance is tied to a scattering block's port current that no equation determines; the circuit is singular in time."))
    directions = Matrix(Z0*Valpha)
    # the constraints are the combinations of the equations without any
    # rate: the block currents cancel in every left null vector and the
    # rates along the islands too, and the combinations without the rate
    # of any other node are kept, as many as the directions
    rows = Matrix(transpose(Z0*rs.leftnull[1:k0, :] .+ Ea*rs.leftnull[k0 + 1:k0 + na, :]))
    leak = rows*G
    cl = size(rows, 1) == 0 ? zeros(0, 0) : nullspace(Matrix(transpose(leak)); atol = 1e-8*max(norm(G, Inf), floatmin(Float64)))
    constraints = size(rows, 1) == 0 ? zeros(0, n) : Matrix(transpose(cl)*rows)
    size(constraints, 1) == d || throw(ArgumentError(
        "a scattering block ties the rate of a node with capacitance to a constraint on a node without one; that coupling is not supported in time."))
    D0 = [1.0/length(z) for z in islands]
    rateextraction = Matrix(transpose(Valpha)*(D0 .* transpose(Z0)))
    algebraic = [sort!([node for (c, z) in enumerate(islands) if abs(Valpha[c, j]) > 1e-8 for node in z]) for j in 1:d]
    # the coupled inductor currents, each its own direction
    for k in 1:Naux
        e = zeros(n); e[Nnodal + k] = 1.0
        directions = hcat(directions, e)
        constraints = vcat(constraints, transpose(e))
        rateextraction = vcat(rateextraction, transpose(e))
        push!(algebraic, [Nnodal + k])
    end
    order = sortperm(algebraic; by = z -> (first(z), length(z)))
    return inertialess, algebraic[order], directions[:, order], constraints[order, :], rateextraction[order, :]
end

# The rate system along the capacitor free islands `Z0` and the block
# current rows `Ea`, balanced by its columns' and rows' magnitudes and
# decomposed: its rank to a relative tolerance, the right and left null
# vectors of the unbalanced system as orthonormal columns, and its
# pseudoinverse on the range, which reads the rates and the currents at
# an endpoint and leaves the null directions alone. The balancing keeps
# a conductance small in the equations' scale from being taken for zero
# next to a block's rows, while an exact dependence, a through's two
# rows, stays one.
function ratesystem(G::SparseMatrixCSC, L::SparseMatrixCSC, Z0::SparseMatrixCSC, Ea::SparseMatrixCSC; rtol = 1e-8)
    Q = sparse(transpose(hcat(Z0, Ea)))
    K = Matrix(hcat(Q*G*Z0, Q*L*Ea))
    m = size(K, 1)
    m == 0 && return (; K, Minv = zeros(0, 0), rightnull = zeros(0, 0), leftnull = zeros(0, 0))
    dc = [(x = maximum(abs, view(K, :, j)); x > 0 ? 1/x : 1.0) for j in 1:m]
    Kb = K .* transpose(dc)
    dr = [(x = maximum(abs, view(Kb, i, :)); x > 0 ? 1/x : 1.0) for i in 1:m]
    Kb = dr .* Kb
    F = svd(Kb; full = true)
    smax = F.S[1]
    r = count(s -> s > rtol*smax, F.S)
    Minv = (dc .* F.V[:, 1:r])*Diagonal(1 ./ F.S[1:r])*transpose(F.U[:, 1:r] .* dr)
    rightnull = orthonormalcolumns(dc .* F.V[:, r + 1:m])
    leftnull = orthonormalcolumns(dr .* F.U[:, r + 1:m])
    return (; K, Minv, rightnull, leftnull)
end

# an orthonormal basis of the span of the columns of `M`
function orthonormalcolumns(M::AbstractMatrix; rtol = 1e-10)
    size(M, 2) == 0 && return zeros(size(M, 1), 0)
    F = svd(M)
    r = count(s -> s > rtol*max(F.S[1], floatmin(Float64)), F.S)
    return F.U[:, 1:r]
end

# the subnetworks of the nodes that no element of the given types with a
# finite nonzero value, nor any of the extra edges, connects to ground, as
# sorted lists of state indices in the order of their first node
function transientsubnetworks(psc::CompiledCircuit, vvn::Vector, types, edges = Tuple{Int,Int}[])
    Nnodes = psc.Nnodes
    parent = collect(1:Nnodes)
    function findroot(i::Int)
        while parent[i] != i
            parent[i] = parent[parent[i]]
            i = parent[i]
        end
        return i
    end
    for k in eachindex(psc.componenttypes)
        psc.componenttypes[k] in types || continue
        v = vvn[k]
        (v isa Number && isfinite(v) && !iszero(v)) || continue
        a, b = findroot(psc.nodeindices[1, k]), findroot(psc.nodeindices[2, k])
        a == b || (parent[max(a, b)] = min(a, b))
    end
    for (i, j) in edges
        a, b = findroot(i), findroot(j)
        a == b || (parent[max(a, b)] = min(a, b))
    end
    components = Dict{Int,Vector{Int}}()
    for node in 2:Nnodes
        root = findroot(node)
        root == 1 && continue
        push!(get!(components, root, Int[]), node - 1)
    end
    directions = sort!(collect(values(components)); by = first)
    foreach(sort!, directions)
    return directions
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
function transientfloatingcomponents(psc::CompiledCircuit, vvn::Vector)
    Nnodes = psc.Nnodes
    parent = collect(1:Nnodes)
    function findroot(i::Int)
        while parent[i] != i
            parent[i] = parent[parent[i]]
            i = parent[i]
        end
        return i
    end
    for k in eachindex(psc.componenttypes)
        t = psc.componenttypes[k]
        t in (:C, :R, :L, :Lj) || continue
        v = vvn[k]
        (v isa Number && isfinite(v) && !iszero(v)) || continue
        a, b = findroot(psc.nodeindices[1, k]), findroot(psc.nodeindices[2, k])
        a == b || (parent[max(a, b)] = min(a, b))
    end
    # a port of a scattering block or a line connects its two terminals
    # as an element does
    for cb in psc.scatteringblocks, q in eachindex(cb.signalnodes)
        a, b = findroot(cb.signalnodes[q]), findroot(cb.refnodes[q])
        a == b || (parent[max(a, b)] = min(a, b))
    end
    components = Dict{Int,Vector{Int}}()
    for node in 2:Nnodes
        root = findroot(node)
        root == 1 && continue
        push!(get!(components, root, Int[]), node)
    end
    floating = sort!(collect(values(components)); by = first)
    foreach(sort!, floating)
    return floating
end

"""
    transientstate(problem; flux = zeros(...), voltage = zeros(...))

The initial state of a transient from the node fluxes in Weber and the
node voltages in Volts, in the compiled node order without ground: the
pair `(x, v)` of the scaled fluxes `flux/phi0`, augmented with the
auxiliary currents the coupled inductors' constitutive equations imply,
in the problem's fixed units `Lscale*i/phi0`, and normalized into the
gauge of the floating subnetworks, and of the scaled flux rates
`voltage/phi0` with the auxiliary rates the same equations imply. The
default is the zero state. A state does not depend on the step, so the
final state of one solve starts another at any step.
[`transientsolve`](@ref) checks that a state satisfies the algebraic
equations of the circuit at the start; it does not project one that does
not. With transmission lines the state has a third member, `waves`, the
wave leaving each port of each line in sqrt(W), two per line in compiled
order, from the port voltages and `linecurrents`, the direct current
into the first port of each line, zero by default; before the start the
lines carry those waves unchanged.
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
        p.graph.Rbn, 1, p.Nnodal, p.Lscale)
    vc = complex(v)
    mnagaugenormalize!(vc, p.floatingcomponents, [0.0], 1)
    mnainitialauxind!(vc, p.coupledbranches, p.matrices.Lb, p.matrices.Mb,
        p.graph.Rbn, 1, p.Nnodal, p.Lscale)
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
    return (x = xr, v = real.(vc), waves = waves, states = zs)
end

# the initial states of the rational blocks of a state, zero when the
# state has none
initialblockstates(state, p::TransientProblem) = length(state) >= 4 ? Float64.(collect(state[4])) : zeros(blockstates(p))

# the initial waves of a state, zero for a circuit without lines and
# for a state without them
initialwaves(state, p::TransientProblem) = length(state) >= 3 ? Float64.(collect(state[3])) : zeros(2length(p.lines))

# the current each line port forces into its node from the far port's
# wave a delay earlier, `2 q / sqrt(Z)` in Amperes, for the waves `q`
# arriving at the ports
lineforcing(p::TransientProblem, q) = [2q[k]/sqrt(p.lines[(k + 1) ÷ 2].Z) for k in eachindex(q)]
# the waves arriving at the ports before the start: the far ports'
# initial waves, which the lines carry unchanged
arrivingwaves(p::TransientProblem, waves) = [waves[isodd(k) ? k + 1 : k - 1] for k in eachindex(waves)]
