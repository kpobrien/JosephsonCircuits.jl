# Binding values, and assembling from a fixed sparsity pattern.
#
# Everything which runs in a loop in this package is already pattern fixed:
# `freqsubst` shares the input's `colptr` and `rowval` and replaces only
# `nzval`, `spaddkeepzeros` keeps structural zeros so a pattern never depends
# on a value, and the Jacobian assembly writes into a structure built once.
# The front end was the exception -- every matrix was rebuilt from coordinate
# triples on every call, which is why changing one inductance rebuilt every
# sparse matrix in the circuit.
#
# A plan closes that gap. The pattern and the destination of every stamp are
# computed once from the compiled topology; assembly is a scatter-add into a
# preallocated `nzval` with no allocation, no sort and no scan for component
# types. The contributions are accumulated in the order the coordinate form
# would have summed duplicates, so the result is bit for bit what
# `calcnodematrix` produced.

"""
    BoundCircuit

A [`CompiledCircuit`](@ref) with its component values resolved to numbers,
grouped by kind.

Each group is a concrete vector in the compiled group order, so
`capacitors[k]` is the value of the component at flat index
`circuit.capacitors[k]`. Element types are chosen per group exactly as the
matrix assembly chooses them, so a lossy capacitance does not make the
resistances complex.

"""
struct BoundCircuit{TC,TR,TL,TJ,TN,TI,TK}
    circuit::CompiledCircuit
    capacitors::Vector{TC}
    resistors::Vector{TR}
    inductors::Vector{TL}
    junctions::Vector{TJ}
    nonlinearinductors::Vector{TN}
    currentsources::Vector{TI}
    mutualinductors::Vector{TK}
    values::Vector          # the flat component table, resolved
end

function Base.show(io::IO, b::BoundCircuit)
    print(io, "BoundCircuit(", ncomponents(b.circuit), " components: ",
        length(b.capacitors), " C, ", length(b.resistors), " R, ",
        length(b.inductors), " L, ", length(b.junctions), " Lj, ",
        length(b.circuit.ports), " ports)")
end

# The element type the assembly chooses for a group. This is `calcvaluetype`
# restricted to the group's own components rather than scanning the whole
# table for them: the compiled circuit already knows which components those
# are, and the scans are seven passes over every component in the circuit.
# The promotion logic is reproduced exactly, quirks included, because the
# assembled matrices must keep their element types.
function grouptype(values, idx, checkinverse::Bool)
    seen = Dict{DataType,Nothing}()
    valuetype = Nothing
    if !isempty(idx)
        v = values[first(idx)]
        valuetype = typeof(v)
        seen[valuetype] = nothing
        checkinverse && (valuetype = promote_type(typeof(1/v), valuetype))
    end
    for i in idx
        v = values[i]
        if typeof(v) != valuetype
            valuetype = promote_type(typeof(v), valuetype)
            if !haskey(seen, valuetype)
                seen[valuetype] = nothing
                checkinverse && (valuetype = promote_type(typeof(1/v), valuetype))
            end
        end
    end
    return valuetype
end

gather(::Type{T}, values, idx) where {T} = T[values[i] for i in idx]

"""
    bind(c::CompiledCircuit, circuitdefs = Dict{Symbol,Any}())

Resolve the component values of a compiled circuit into concrete per group
arrays.

Symbolic values are substituted from `circuitdefs`; values which depend on
frequency are resolved per mode later and are carried through unevaluated.

A bound circuit records the assumptions the compiled structure rests on --
which inductances are finite, which values are complex, which reference
impedances are positive -- so that rebinding at new values can tell a change
which only moves numbers from one which changes the topology. See
[`structuralkey`](@ref).
"""
bind(c::CompiledCircuit, circuitdefs = Dict{Symbol,Any}()) =
    bindvalues(c, componentvaluestonumber(c.componentvalues, circuitdefs))

"""
    bindvalues(c::CompiledCircuit, values)

A [`BoundCircuit`](@ref) from an already resolved flat value table.

This is the entry point for rebinding: the topology, the groups and the
assembly plans are unchanged when only the numbers move, so a sweep resolves
its values once per point and gathers them into the groups without touching
anything structural. Check [`structuralkey`](@ref) before reusing a plan
across a rebind.
"""
function bindvalues(c::CompiledCircuit, values)
    TC = grouptype(values, c.capacitors, true)
    TR = grouptype(values, c.resistors, true)
    TL = grouptype(values, c.inductors, true)
    TJ = grouptype(values, c.junctions, true)
    TN = grouptype(values, c.nonlinearinductors, true)
    TI = grouptype(values, c.currentsources, true)
    TK = grouptype(values, c.mutualinductors, true)
    return BoundCircuit(c,
        gather(TC, values, c.capacitors),
        gather(TR, values, c.resistors),
        gather(TL, values, c.inductors),
        gather(TJ, values, c.junctions),
        gather(TN, values, c.nonlinearinductors),
        gather(TI, values, c.currentsources),
        gather(TK, values, c.mutualinductors), values)
end

"""
    structuralkey(b::BoundCircuit)

The structural facts a compiled plan depends on, as a comparable value.

A plan may be reused at new component values only while this is unchanged.
The facts are the ones which move a component between behaviors rather than
along a range: an infinite inductance is an open circuit and drops an edge
from the static flux graph, a zero one is a short, a complex capacitance is a
noise channel where a real one is not, and a unit mutual coupling makes the
inverse inductance matrix singular. A value crossing any of these changes the
sparsity pattern, the branch set, or the noise classification, so the plan
must be rebuilt rather than refilled.
"""
function structuralkey(b::BoundCircuit)
    finite(v) = v isa Number ? isfinite(abs(v)) : true
    nonzero(v) = v isa Number ? !iszero(v) : true
    complexvalued(v) = v isa Number && !iszero(imag(v))
    return (
        inductoropen = Bool[!finite(v) for v in b.inductors],
        inductorshort = Bool[!nonzero(v) for v in b.inductors],
        junctionopen = Bool[!finite(v) for v in b.junctions],
        junctionshort = Bool[!nonzero(v) for v in b.junctions],
        resistoropen = Bool[!finite(v) for v in b.resistors],
        resistorshort = Bool[!nonzero(v) for v in b.resistors],
        capacitorlossy = Bool[complexvalued(v) for v in b.capacitors],
        inductorlossy = Bool[complexvalued(v) for v in b.inductors],
        unitcoupling = Bool[v isa Number && isone(abs(v))
            for v in b.mutualinductors],
    )
end

# === nodal stamp plans ===

"""
    NodalStampPlan

The fixed pattern and stamp destinations of a nodal matrix.

`dest[k]` is the position in `nzval` which contribution `k` accumulates into,
`src[k]` names the group value it comes from, and `negate[k]` marks the off
diagonal terms. `invert` inverts the value after negating, which is the order
[`calcnodematrix`](@ref) uses and therefore the one that reproduces its
arithmetic exactly.

The plan is built once per topology; [`assemblenodal!`](@ref) is a
scatter-add which allocates nothing.
"""
struct NodalStampPlan{Ti<:Integer}
    colptr::Vector{Ti}
    rowval::Vector{Ti}
    dest::Vector{Ti}
    src::Vector{Ti}
    negate::Vector{Bool}
    invert::Bool
    n::Int
end

"""
    nodalstampplan(c::CompiledCircuit, group, Nnodes; invert = false)

Build the [`NodalStampPlan`](@ref) of a two terminal group.

`group` is a vector of flat component indices, so the same function plans the
capacitance from the capacitors and the conductance from the resistors and
the port environments alike.
"""
function nodalstampplan(c::CompiledCircuit, group::Vector{Int}, Nnodes::Int;
        invert::Bool = false)

    n = Nnodes - 1
    # the coordinate form the old assembly built, in the order it built it:
    # one diagonal entry for a grounded component, two diagonal and two off
    # diagonal for a floating one
    I = Int[]; J = Int[]; S = Int[]; G = Bool[]
    for (k, i) in enumerate(group)
        n1, n2 = c.nodeindices[1, i], c.nodeindices[2, i]
        if n1 == 1
            push!(I, n2-1); push!(J, n2-1); push!(S, k); push!(G, false)
        elseif n2 == 1
            push!(I, n1-1); push!(J, n1-1); push!(S, k); push!(G, false)
        else
            push!(I, n1-1); push!(J, n1-1); push!(S, k); push!(G, false)
            push!(I, n2-1); push!(J, n2-1); push!(S, k); push!(G, false)
            push!(I, n1-1); push!(J, n2-1); push!(S, k); push!(G, true)
            push!(I, n2-1); push!(J, n1-1); push!(S, k); push!(G, true)
        end
    end

    # the pattern, and where each contribution lands in it. Building it from
    # the same coordinate triples the old assembly used is what makes the two
    # patterns identical, including any structural zero.
    pattern = sparse(I, J, ones(Int, length(I)), n, n)
    dest = Vector{Int}(undef, length(I))
    for k in eachindex(I)
        # I and J are already indices into the grounded-node-removed matrix
        col = J[k]
        r = searchsortedfirst(view(pattern.rowval,
            pattern.colptr[col]:(pattern.colptr[col+1]-1)), I[k])
        dest[k] = pattern.colptr[col] + r - 1
    end
    return NodalStampPlan(pattern.colptr, pattern.rowval, dest, S, G, invert, n)
end

"""
    assemblenodal!(nzval, seen, plan::NodalStampPlan, values)

Accumulate the stamps of `values` into `nzval` against a fixed pattern.

Contributions are summed in the order the coordinate form would have combined
duplicates, so two components on one node pair add in the same order and the
result is bit for bit the old assembly's.
"""
function assemblenodal!(nzval::Vector, seen::Vector{Bool},
        plan::NodalStampPlan, values)
    # The first contribution to a position is assigned and later ones are
    # added to it, rather than accumulating onto a zero. That is what the
    # coordinate form did -- `sparse` combines duplicates starting from the
    # first value -- and it is why this works for element types which have
    # no zero: an empty group, whose element type is Nothing, and a circuit
    # whose values are still symbolic, which has to reach the diagnostic
    # which names the undefined one rather than failing here.
    fill!(seen, false)
    @inbounds for k in eachindex(plan.dest)
        v = values[plan.src[k]]
        plan.negate[k] && (v = -v)
        plan.invert && (v = 1/v)
        d = plan.dest[k]
        nzval[d] = seen[d] ? nzval[d] + v : v
        seen[d] = true
    end
    return nzval
end

"""
    assemblenodal(::Type{T}, plan::NodalStampPlan, values, Nmodes)

The nodal matrix of `values`, repeated along the diagonal for `Nmodes`.
"""
function assemblenodal(::Type{T}, plan::NodalStampPlan, values,
        Nmodes::Integer) where {T}
    nzval = Vector{T}(undef, length(plan.rowval))
    assemblenodal!(nzval, Vector{Bool}(undef, length(plan.rowval)), plan,
        values)
    A = SparseMatrixCSC(plan.n, plan.n, plan.colptr, plan.rowval, nzval)
    return Nmodes == 1 ? A : diagrepeat(A, Nmodes)
end

# === branch stamp plans ===

"""
    BranchStampPlan

The branch each component of a group occupies, and where it lands in the
branch vector.

`nzind` is the sorted list of branches the group touches and `dest[k]` is the
position of group member `k` within it. Two components on one branch share a
destination and are folded together in the order they appear, which is the
order the coordinate form combined them.

The branch of a component is a dictionary lookup on its node pair, and there
is one per component per assembly; doing it once is most of what this plan
saves.
"""
struct BranchStampPlan{Ti<:Integer}
    nzind::Vector{Ti}
    dest::Vector{Ti}
    n::Int
end

"""
    branchstampplan(c::CompiledCircuit, group, edge2indexdict, Nbranches)

Build the [`BranchStampPlan`](@ref) of a two terminal group.
"""
function branchstampplan(c::CompiledCircuit, group::Vector{Int},
        edge2indexdict::Dict, Nbranches::Int)
    branch = [edge2indexdict[(c.nodeindices[1,i], c.nodeindices[2,i])]
              for i in group]
    nzind = sort(unique(branch))
    position = Dict(b => k for (k, b) in enumerate(nzind))
    dest = [position[b] for b in branch]
    return BranchStampPlan(nzind, dest, Nbranches)
end

"""
    assemblebranch!(nzval, seen, plan::BranchStampPlan, values, combine)

Fold `values` into the branch vector `nzval` against a fixed set of branches.

`combine` is applied to the running value and the new one in the order the
components appear, matching `sparsevec`'s combination of duplicate indices:
two inductors on one branch combine as a parallel inductance, and two
junctions raise the error that says to separate them.
"""
function assemblebranch!(nzval::Vector, seen::Vector{Bool},
        plan::BranchStampPlan, values, combine::F) where {F}
    fill!(seen, false)
    @inbounds for k in eachindex(plan.dest)
        d = plan.dest[k]
        v = values[k]
        nzval[d] = seen[d] ? combine(nzval[d], v) : v
        seen[d] = true
    end
    return nzval
end

"""
    assemblebranch(::Type{T}, plan::BranchStampPlan, values, combine, Nmodes)

The branch vector of `values`, repeated along the diagonal for `Nmodes`.
"""
function assemblebranch(::Type{T}, plan::BranchStampPlan, values,
        combine::F, Nmodes::Integer) where {T,F}
    nzval = Vector{T}(undef, length(plan.nzind))
    seen = Vector{Bool}(undef, length(plan.nzind))
    assemblebranch!(nzval, seen, plan, values, combine)
    v = SparseVector(plan.n, plan.nzind, nzval)
    return Nmodes == 1 ? v : diagrepeat(v, Nmodes)
end

# === inverse nodal inductance ===
#
# The solvers represent mutually coupled inductor branches by auxiliary MNA
# branch currents with their un-inverted branch inductance matrix as explicit
# constitutive equations, so those branches are dropped here and no
# inductance matrix is ever inverted. What is left,
#
#     transpose(Rbn) * diagm(1/L) * Rbn
#
# over the retained inductive branches, is the same stamp shape as a nodal
# capacitance: a branch deposits +1/L on the diagonal of each of its nodes
# and -1/L on the two off diagonal entries. So it is planned with the same
# machinery, reading the node pair and the signs out of the incidence matrix
# and taking the reciprocals of the branch inductances as its values. That
# is also the order the triple product forms them in -- reciprocal first,
# then the sign of the incidence product -- so the arithmetic matches.

"""
    InverseInductancePlan

The fixed pattern of the inverse nodal inductance matrix, and which branch
inductances feed it.

`positions` selects the retained branches out of the assembled branch
inductance vector: the coupled ones are dropped, because the solvers carry
them as auxiliary MNA currents instead.
"""
struct InverseInductancePlan{Ti<:Integer}
    stamp::NodalStampPlan{Ti}
    positions::Vector{Int}
end

"""
    inverseinductanceplan(c, cg, Lb, coupled)

Build the [`InverseInductancePlan`](@ref) from the incidence matrix, the
branches which carry an inductance, and the branches which are mutually
coupled.
"""
function inverseinductanceplan(c::CompiledCircuit, cg::CircuitGraph,
        Lb::SparseVector, coupled)
    coupledset = Set(coupled)
    positions = [k for (k, b) in enumerate(Lb.nzind) if !(b in coupledset)]
    branches = Lb.nzind[positions]
    Rbn = cg.Rbn
    n = size(Rbn, 2)
    I = Int[]; J = Int[]; S = Int[]; G = Bool[]
    Rt = sparse(transpose(Rbn))    # columns are branches, so a row is one lookup
    for (k, b) in enumerate(branches)
        nodes = Rt.rowval[Rt.colptr[b]:(Rt.colptr[b+1]-1)]
        signs = Rt.nzval[Rt.colptr[b]:(Rt.colptr[b+1]-1)]
        for (bi, sb) in zip(nodes, signs), (ai, sa) in zip(nodes, signs)
            push!(I, ai); push!(J, bi); push!(S, k); push!(G, sa*sb < 0)
        end
    end
    pattern = sparse(I, J, ones(Int, length(I)), n, n)
    dest = Vector{Int}(undef, length(I))
    for k in eachindex(I)
        col = J[k]
        r = searchsortedfirst(view(pattern.rowval,
            pattern.colptr[col]:(pattern.colptr[col+1]-1)), I[k])
        dest[k] = pattern.colptr[col] + r - 1
    end
    stamp = NodalStampPlan(pattern.colptr, pattern.rowval, dest, S, G, false, n)
    return InverseInductancePlan(stamp, positions)
end

"""
    assembleinvinductance(::Type{T}, plan, Lb, Nmodes)

The inverse nodal inductance matrix of the branch inductances `Lb`.
"""
function assembleinvinductance(::Type{T}, plan::InverseInductancePlan,
        Lb::SparseVector, Nmodes::Integer) where {T}
    isempty(plan.positions) &&
        return spzeros(T, Nmodes*plan.stamp.n, Nmodes*plan.stamp.n)
    invL = T[1/Lb.nzval[p] for p in plan.positions]
    return assemblenodal(T, plan.stamp, invL, Nmodes)
end

# === the circuit matrices ===

"""
    CircuitMatrixPlan

Everything about a circuit's matrices which depends on its topology but not
on its values, for one mode count.

Holds the nodal and branch stamp plans and the mode expanded incidence
matrix. Rebinding at new component values reuses all of it; only the values
are refilled. See [`circuitmatrixplan`](@ref) and
[`assemblematrices`](@ref).
"""
struct CircuitMatrixPlan{Ti<:Integer}
    circuit::CompiledCircuit
    graph::CircuitGraph
    Nmodes::Int
    capacitance::NodalStampPlan{Ti}
    conductance::NodalStampPlan{Ti}
    inductance::BranchStampPlan{Ti}
    junction::BranchStampPlan{Ti}
    invinductance::InverseInductancePlan{Ti}
    Rbnm::SparseMatrixCSC{Int,Int}
end

"""
    circuitmatrixplan(c::CompiledCircuit, cg::CircuitGraph; Nmodes = 1)

Build the [`CircuitMatrixPlan`](@ref) of a compiled circuit.

The conductance plan covers the resistors and the port owned environments
together, because at this stage an environment is realized as an ordinary
resistor; when ports become direct boundary stamps it gains its own plan.
"""
function circuitmatrixplan(c::CompiledCircuit, cg::CircuitGraph,
        b::BoundCircuit; Nmodes::Int = 1)
    inductance = branchstampplan(c, c.inductors, cg.edge2indexdict,
        cg.Nbranches)
    # which branches are mutually coupled is a structural fact, but it is
    # read off the mutual inductance matrix, so the plan is built from a
    # bound circuit and is valid while `structuralkey` is unchanged
    Mb = calcMb(c.componenttypes, c.nodeindices, b.values,
        c.componentnamedict, c.mutualinductorbranchnames, cg.edge2indexdict,
        1, cg.Nbranches)
    Lb = assemblebranch(eltype(b.inductors), inductance, b.inductors,
        combine_reciprocal_sum, 1)
    return CircuitMatrixPlan(c, cg, Nmodes,
        nodalstampplan(c, c.capacitors, c.Nnodes),
        nodalstampplan(c, c.resistors, c.Nnodes; invert = true),
        inductance,
        branchstampplan(c, c.junctions, cg.edge2indexdict, cg.Nbranches),
        inverseinductanceplan(c, cg, Lb, mnacoupledbranches(Mb)),
        diagrepeat(cg.Rbn, Nmodes))
end

"""
    assemblematrices(plan::CircuitMatrixPlan, b::BoundCircuit)

The [`CircuitMatrices`](@ref) of a bound circuit, assembled against the fixed
patterns of `plan`.

Numerically identical to [`numericmatrices`](@ref) on the same circuit, entry
for entry. The parts which are cheap to rebuild and depend on the values in
ways a stamp plan does not express -- the mutual inductance matrix, the
inverse inductance triple product, the solver scale and the port index
lists -- are still computed the same way; the plan covers the four which
dominate.
"""
function assemblematrices(plan::CircuitMatrixPlan, b::BoundCircuit)
    c = plan.circuit
    cg = plan.graph
    vvn = b.values
    Nmodes = plan.Nmodes
    ct, ni = c.componenttypes, c.nodeindices

    TC = eltype(b.capacitors)
    TR = eltype(b.resistors)
    TL = eltype(b.inductors)
    TJ = eltype(b.junctions)

    Cnm = assemblenodal(TC, plan.capacitance, b.capacitors, Nmodes)
    Gnm = assemblenodal(TR, plan.conductance, b.resistors, Nmodes)
    Lb = assemblebranch(TL, plan.inductance, b.inductors,
        combine_reciprocal_sum, 1)
    Lbm = assemblebranch(TL, plan.inductance, b.inductors,
        combine_reciprocal_sum, Nmodes)
    Ljb = assemblebranch(TJ, plan.junction, b.junctions, combine_error, 1)
    Ljbm = assemblebranch(TJ, plan.junction, b.junctions, combine_error,
        Nmodes)

    Mb = calcMb(ct, ni, vvn, c.componentnamedict,
        c.mutualinductorbranchnames, cg.edge2indexdict, 1, cg.Nbranches)
    checkcoupledbranchinductors(c.componentnames, ct, ni, cg.edge2indexdict,
        Mb)
    invLnm = assembleinvinductance(TL, plan.invinductance, Lb, Nmodes)

    Lmean = calcLmean(ct, vvn)
    # by role, not by geometry: the port states its reference impedance and
    # what realizes it, and nothing here looks at what shares its branch
    portindices, portnumbers = portindicesnumbers(c)
    portimpedances = portreferenceimpedances(c, vvn)
    noiseportimpedanceindices = noiseindices(c, vvn)

    return CircuitMatrices(Cnm, Gnm, Lb, Lbm, Ljb, Ljbm, Mb, invLnm,
        plan.Rbnm, portindices, portnumbers, portimpedances,
        portenvironmentindices(c), noiseportimpedanceindices, Lmean, vvn)
end

# === the solver front end ===

"""
    compilable(circuit)

Whether a circuit can take the compiled path. A netlist of tuples and an
already parsed circuit cannot, and take the original one.
"""
compilable(::Union{Circuit,ElaboratedCircuit}) = true
compilable(::AbstractVector) = true
compilable(x) = false

"""
    preparecircuit(circuit, circuitdefs; sorting = :name, Nmodes = 1)

Everything the solvers need from a circuit: the compiled circuit, its graph,
and its numeric matrices.

Every circuit takes this path. A netlist of tuples becomes a typed circuit
first, which parses to the same thing at the same sorting and refuses the
same inputs for the same reasons.
"""
function preparecircuit(circuit, circuitdefs; sorting::Symbol = :name,
        Nmodes::Int = 1)
    c = compile(circuit; sorting = sorting)
    # the loop enumeration is quadratic in the number of inductive loops and
    # nothing reads it
    cg = calccircuitgraph(c; loops = false)
    b = bind(c, circuitdefs)
    nm = assemblematrices(circuitmatrixplan(c, cg, b; Nmodes = Nmodes), b)
    return c, cg, nm
end

"""
    assemblegrid(compiled, cg, circuitdefs, Nmodes)

The numeric matrices of a compiled circuit at a different mode count.

The pump and the signal sweep of one analysis use different mode grids, so
each assembles its own; the topology and the values behind them are the
same.
"""
function assemblegrid(c::CompiledCircuit, cg::CircuitGraph, circuitdefs,
        Nmodes::Integer)
    b = circuitdefs isa AbstractVector ? bindvalues(c, circuitdefs) :
        bind(c, circuitdefs)
    return assemblematrices(circuitmatrixplan(c, cg, b; Nmodes = Nmodes), b)
end

"""
    scatteringstampsystem(blocks::Vector{CompiledScatteringBlock}, Nmodes;
        auxoffset, Ntotal, scale = 1.0)

The stamp system of the compiled scattering blocks of a circuit.

A compiled block is one instance carrying its own terminal map, so this
needs no regrouping and none of the checks which the per port form does: a
block cannot be missing a port, cannot repeat one, and cannot be confused
with another instance of the same definition.
"""
function scatteringstampsystem(blocks::Vector{CompiledScatteringBlock},
    Nmodes::Integer; auxoffset::Integer, Ntotal::Integer,
    scale::Real = 1.0)

    isempty(blocks) && return nothing
    stamped = StampedScatteringBlock[]
    auxbase = auxoffset
    for b in blocks
        n = b.definition.nports
        push!(stamped, StampedScatteringBlock(b.definition,
            copy(b.signalnodes), copy(b.refnodes), auxbase,
            b.path*"/port1"))
        auxbase += n*Nmodes
    end
    return scatteringstampsystem(stamped, Nmodes, Ntotal, scale)
end
