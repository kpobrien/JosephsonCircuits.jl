# Binding values, and assembling the circuit matrices from a fixed
# sparsity pattern.
#
# The circuit matrices are assembled on a plan: the pattern and the
# destination of every stamp are computed once from the compiled topology,
# and an assembly is a scatter-add into a preallocated `nzval` with no
# allocation, no sort and no scan for component types, for
# `numericmatrices` and the solvers alike. This matches the rest of the
# solve, where `freqsubst` replaces only `nzval`, `spaddkeepzeros` keeps
# structural zeros so a pattern never depends on a value, and the Jacobian
# assembly writes into a structure built once. The contributions of parallel
# components are summed in the order the netlist lists them. A plan holds
# the patterns and the resolved couplings and ports; the scratch of a refill
# is a solve's own.

"""
    FrequencyDependentComponents

The capacitors, resistors and inductors of a bound circuit whose values
are [`FrequencyDependent`](@ref): the distinct functions of their values,
as `===` tells them apart, in `providers`, and for each group the position
of each such component in it and the index of its function.
"""
struct FrequencyDependentComponents
    providers::Vector{Any}
    capacitors::Vector{Tuple{Int,Int}}
    resistors::Vector{Tuple{Int,Int}}
    inductors::Vector{Tuple{Int,Int}}
end

"""
    BoundCircuit

A [`CompiledCircuit`](@ref) with its component values resolved to numbers,
grouped by kind.

Each group is a concrete vector in the compiled group order, so
`capacitors[k]` is the value of the component at flat index
`circuit.capacitors[k]`. Element types are chosen per group exactly as the
matrix assembly chooses them, so a lossy capacitance does not make the
resistances complex.

A component whose value depends on the frequency holds a value in its
group which leaves it open, a zero capacitance or an infinite resistance
or inductance, so that the matrices assembled from the groups describe the
circuit without it; `frequencydependent` records it, and the terms of
[`FrequencyDependentTerms`](@ref) add it at each mode frequency.

# Fields
- `circuit`: the [`CompiledCircuit`](@ref).
- `capacitors`, `resistors`, `inductors`, `junctions`, `mutualinductors`:
    the values of each group, in the compiled group order, every nonlinear
    inductor's among the junctions. The current sources are read from
    `values`, so they have no group here.
- `values`: the flat component table, resolved to numbers.
- `frequencydependent`: the [`FrequencyDependentComponents`](@ref), or
    `nothing` when no value depends on the frequency.
"""
struct BoundCircuit{TC,TR,TL,TJ,TK}
    circuit::CompiledCircuit
    capacitors::Vector{TC}
    resistors::Vector{TR}
    inductors::Vector{TL}
    junctions::Vector{TJ}
    mutualinductors::Vector{TK}
    values::Vector          # the flat component table, resolved
    frequencydependent::Union{Nothing,FrequencyDependentComponents}
end

# the frequency dependent capacitors, resistors and inductors of the flat
# value table `values`, or `nothing` when there are none; a frequency
# dependent junction, mutual coupling or coupled inductor is refused
function frequencydependentcomponents(c::CompiledCircuit, values)
    any(v -> v isa FrequencyDependent, values) || return nothing
    checkfrequencyindependent(c, values)
    providers = Any[]
    index = IdDict{Any,Int}()
    function members(group)
        out = Tuple{Int,Int}[]
        for (k, i) in enumerate(group)
            v = values[i]
            v isa FrequencyDependent || continue
            j = get!(index, v.f) do
                push!(providers, v.f)
                length(providers)
            end
            push!(out, (k, j))
        end
        return out
    end
    C, R, L = members(c.capacitors), members(c.resistors), members(c.inductors)
    isempty(providers) && return nothing
    return FrequencyDependentComponents(providers, C, R, L)
end

function Base.show(io::IO, b::BoundCircuit)
    print(io, "BoundCircuit(", ncomponents(b.circuit), " components: ",
        length(b.capacitors), " C, ", length(b.resistors), " R, ",
        length(b.inductors), " L, ", length(b.junctions), " Lj, ",
        length(b.circuit.ports), " ports)")
end

# The element type of a value group. Numbers are held as `Float64` however
# the netlist wrote them (`1` or `1.0f0` for an inductance, `50` for a
# resistance, `50 + 0im` in a table of complex numbers), empty groups
# included, and as `ComplexF64` only where a value has an imaginary part,
# so that a circuit has one matrix type per real or complex distinction
# of its values rather than one per way of writing them; every solver
# method downstream is compiled once per matrix type, and a circuit whose
# values are real runs on the code compiled for real values whatever
# table they came in. A frequency dependent value is held as the value
# which leaves its component open (see `gather`) and adds no type. A group
# with a symbolic value keeps the promotion of the types present and of
# their inverses; a value which is neither a number nor symbolic
# (`nothing`, a vector) is refused there, by the component `names` gives
# it.
function grouptype(values, idx, names = eachindex(values))
    isempty(idx) && return Float64
    complex = false
    plain = true
    for i in idx
        v = values[i]
        v isa FrequencyDependent && continue
        if v isa Complex && plainnumber(real(v))
            complex |= !iszero(imag(v))
        elseif !plainnumber(v)
            plain = false
            break
        end
    end
    plain && return complex ? ComplexF64 : Float64
    valuetype = Union{}
    for i in idx
        v = values[i]
        v isa FrequencyDependent && continue
        v isa Number || checkissymbolic(v) || throw(ArgumentError(
            lazy"The component $(names[i]) has the value $(repr(v)), which is not a number, a name or a frequency dependent value."))
        valuetype = promote_type(valuetype, typeof(v))
        # the inverse too: the conductance of a symbolic resistance is a
        # different expression type from the resistance itself
        valuetype = promote_type(valuetype, typeof(1/v))
    end
    return valuetype
end

plainnumber(v) = v isa Union{Integer,AbstractFloat,Rational,Irrational}

# the values of a group, a frequency dependent one as the value which
# leaves its component open: a zero capacitance, an infinite resistance or
# inductance
gather(::Type{T}, values, idx, open) where {T} =
    T[values[i] isa FrequencyDependent ? open : values[i] for i in idx]

# the mean of the linear and Josephson inductances, zero when the circuit
# has none. A frequency dependent inductance has no one value and is left
# out. An empty group is left out of the sum rather than started from its
# own zero, so the mean is the sum's own type and a value type with no
# zero of its own, such as a symbolic expression, is still summed.
function inductancemean(b::BoundCircuit)
    fd = b.frequencydependent
    inductors = isnothing(fd) || isempty(fd.inductors) ? b.inductors :
        b.inductors[setdiff(eachindex(b.inductors), first.(fd.inductors))]
    junctions = b.junctions
    n = length(inductors) + length(junctions)
    n == 0 && return 0.0
    total = if isempty(junctions)
        sum(inductors)
    elseif isempty(inductors)
        sum(junctions)
    else
        sum(inductors) + sum(junctions)
    end
    return total/n
end

# The inductance of a junction, a mutual coupling and the inductances it
# couples are read once, by the nonlinear term and the mutual stamps,
# rather than at each mode frequency, so none of them may depend on
# frequency.
function checkfrequencyindependent(c::CompiledCircuit, values)
    for i in c.junctions
        values[i] isa FrequencyDependent && throw(ArgumentError(
            lazy"The junction $(c.componentnames[i]) has the frequency dependent value $(values[i]); the inductance of a junction cannot depend on frequency."))
    end
    for (k, i, j) in c.couplings, m in (k, i, j)
        values[m] isa FrequencyDependent || continue
        kind = m == k ? "mutual coupling" : "coupled inductor"
        throw(ArgumentError(lazy"The $(kind) $(c.componentnames[m]) has the frequency dependent value $(values[m]); a mutual coupling and the inductances it couples cannot depend on frequency."))
    end
    return nothing
end

"""
    numericvalues(c::CompiledCircuit, definitions)

The flat value table of a compiled circuit resolved at `definitions`, with
every entry a number. A value which still depends on a name the
definitions do not give, and one which is frequency dependent rather than
a number, are both refused naming their component, so a caller which can
work with neither fails with the cause rather than downstream. The
definitions may be any dictionary.
"""
function numericvalues(c::CompiledCircuit, definitions)
    vvn = componentvaluestonumber(c.componentvalues,
        definitiontable(definitions))
    checkcomponentvaluesdefined(c.componentnames, vvn)
    for (i, v) in enumerate(vvn)
        v isa Number || throw(ArgumentError(
            lazy"the component $(c.componentnames[i]) has the value $(v) at these definitions, which is not a number."))
    end
    return vvn
end

"""
    bindvalues(c::CompiledCircuit, values)

A [`BoundCircuit`](@ref) from the flat value table `values`, the component
values resolved with [`componentvaluestonumber`](@ref): the values gathered
into the per group arrays, each in the element type its group promotes to.
A name the definitions did not give, a value which is neither a number nor
symbolic, and a frequency dependent junction, mutual coupling or coupled
inductor are refused naming their components.

The topology, the groups and the assembly plans do not depend on the
values, so a sweep binds its values at each point and refills the matrices
on the same plan, as long as each group keeps its element type.
"""
function bindvalues(c::CompiledCircuit, values)
    length(values) == ncomponents(c) || throw(DimensionMismatch(
        "componenttypes and componentvalues should have the same length"))
    names = c.componentnames
    # a name the definitions did not give is held by no group, port or
    # source: reported here with every other undefined value
    any(v -> v isa Union{Symbol,AbstractString}, values) &&
        checkcomponentvaluesdefined(names, values)
    TC = grouptype(values, c.capacitors, names)
    TL = grouptype(values, c.inductors, names)
    TR = grouptype(values, c.resistors, names)
    TJ = grouptype(values, c.junctions, names)
    TK = grouptype(values, c.mutualinductors, names)
    return BoundCircuit(c,
        gather(TC, values, c.capacitors, zero(TC)),
        gather(TR, values, c.resistors, convert(TR, Inf)),
        gather(TL, values, c.inductors, convert(TL, Inf)),
        gather(TJ, values, c.junctions, convert(TJ, Inf)),
        gather(TK, values, c.mutualinductors, zero(TK)), values,
        frequencydependentcomponents(c, values))
end

# === nodal stamp plans ===

"""
    NodalStampPlan

The fixed pattern and stamp destinations of a nodal matrix.

`dest[k]` is the position in `nzval` which contribution `k` accumulates into,
`src[k]` names the group value it comes from, and `weights[k]` holds its
incidence coefficient, `1` on the diagonal and `-1` off it. `invert`
inverts the value after applying that coefficient, so a conductance is
`-1/R` off the diagonal.

The plan is built once per topology; [`assemblenodal!`](@ref) is a
scatter-add which allocates nothing.
"""
struct NodalStampPlan{Ti<:Integer}
    colptr::Vector{Ti}
    rowval::Vector{Ti}
    dest::Vector{Ti}
    src::Vector{Ti}
    weights::Vector{Int8}
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
    # one diagonal entry for a grounded component, two diagonal and two off
    # diagonal for a floating one, and none for a component whose terminals
    # are one node, which carries no current
    I = Int[]; J = Int[]; S = Int[]; G = Int8[]
    for (k, i) in enumerate(group)
        n1, n2 = c.nodeindices[1, i], c.nodeindices[2, i]
        n1 == n2 && continue
        if n1 == 1
            push!(I, n2-1); push!(J, n2-1); push!(S, k); push!(G, 1)
        elseif n2 == 1
            push!(I, n1-1); push!(J, n1-1); push!(S, k); push!(G, 1)
        else
            push!(I, n1-1); push!(J, n1-1); push!(S, k); push!(G, 1)
            push!(I, n2-1); push!(J, n2-1); push!(S, k); push!(G, 1)
            push!(I, n1-1); push!(J, n2-1); push!(S, k); push!(G, -1)
            push!(I, n2-1); push!(J, n1-1); push!(S, k); push!(G, -1)
        end
    end

    # the pattern, and where each contribution lands in it, structural
    # zeros included
    return nodalstampplan(I, J, S, G, invert, n)
end

# the pattern of a set of coordinate contributions and where each lands in
# it, which the nodal, the inverse inductance and the mutual stamps share
function nodalstampplan(I, J, S, G, invert::Bool, n::Int)
    pattern = sparse(I, J, ones(Int, length(I)), n, n)
    dest = Vector{Int}(undef, length(I))
    for k in eachindex(I)
        dest[k] = nzposition(pattern, I[k], J[k])
    end
    return NodalStampPlan(pattern.colptr, pattern.rowval, dest, S, G, invert, n)
end

"""
    assemblenodal!(nzval, seen, plan::NodalStampPlan, values)

Accumulate the stamps of `values` into `nzval` against a fixed pattern.

Contributions to one position are summed in the order the components
appear, the order in which `sparse` combines duplicate coordinates.
"""
function assemblenodal!(nzval::Vector, seen::Vector{Bool},
        plan::NodalStampPlan, values)
    # The first contribution to a position is assigned and later ones are
    # added to it, rather than accumulating onto a zero, as `sparse` does
    # when it combines duplicates. This is what makes it work for element
    # types with no zero, a circuit whose values are still symbolic, which
    # must reach the diagnostic naming the undefined value rather than fail
    # here.
    fill!(seen, false)
    @inbounds for k in eachindex(plan.dest)
        v = values[plan.src[k]]
        plan.weights[k] == -1 && (v = -v)
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
destination and are folded together in the order they appear.

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

# the plan of the junctions, which refuses two on one branch by name: they
# cannot be combined into one element
function junctionstampplan(c::CompiledCircuit, edge2indexdict::Dict,
        Nbranches::Int)
    plan = branchstampplan(c, c.junctions, edge2indexdict, Nbranches)
    first = zeros(Int, length(plan.nzind))
    for (k, d) in enumerate(plan.dest)
        if !iszero(first[d])
            a = c.componentnames[c.junctions[first[d]]]
            b = c.componentnames[c.junctions[k]]
            throw(ArgumentError(lazy"The Josephson junctions $(a) and $(b) are on one branch and cannot be combined into a single element. Place them between different nodes."))
        end
        first[d] = k
    end
    return plan
end

"""
    assemblebranch!(nzval, seen, plan::BranchStampPlan, values, combine)

Fold `values` into the branch vector `nzval` against a fixed set of branches.

`combine` is applied to the running value and the new one in the order the
components appear, matching `sparsevec`'s combination of duplicate indices:
two inductors on one branch combine as a parallel inductance. With
`combine = nothing` every branch holds one component, as the junction plan
guarantees (see `junctionstampplan`), and each value is written to its
branch, with no `seen` to keep.
"""
function assemblebranch!(nzval::Vector, seen, plan::BranchStampPlan,
        values, combine::F) where {F}
    if isnothing(combine)
        @inbounds for k in eachindex(plan.dest)
            nzval[plan.dest[k]] = values[k]
        end
        return nzval
    end
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
    seen = isnothing(combine) ? nothing : Vector{Bool}(undef, length(plan.nzind))
    assemblebranch!(nzval, seen, plan, values, combine)
    v = SparseVector(plan.n, plan.nzind, nzval)
    return Nmodes == 1 ? v : diagrepeat(v, Nmodes)
end

# === inverse nodal inductance ===
#
# The solvers represent mutually coupled inductor branches by auxiliary
# branch currents with their (uninverted) branch inductance matrix as
# explicit constitutive equations, so those branches are dropped here and
# no inductance matrix is inverted. What is left,
#
#     transpose(Rbn) * diagm(1/L) * Rbn
#
# over the retained inductive branches, is the same stamp shape as a nodal
# capacitance: a branch deposits +1/L on the diagonal of each of its nodes
# and -1/L on the two off diagonal entries. So it is planned with the same
# machinery, reading the node pair and the signs out of the incidence
# matrix and taking the reciprocals of the branch inductances as values,
# formed in the order the triple product forms them (reciprocal first,
# then the sign of the incidence product) so the arithmetic matches.

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
    inverseinductanceplan(topology, inductivebranches, coupled)

Build the [`InverseInductancePlan`](@ref) from the incidence matrix of
`topology`, the branches which carry an inductance, and the branches which
are mutually coupled.
"""
function inverseinductanceplan(topology::CircuitTopology,
        inductivebranches::Vector{Int}, coupled)
    coupledset = Set(coupled)
    positions = [k for (k, b) in enumerate(inductivebranches) if !(b in coupledset)]
    branches = inductivebranches[positions]
    Rbn = topology.Rbn
    n = size(Rbn, 2)
    I = Int[]; J = Int[]; S = Int[]; G = Int8[]
    Rt = sparse(transpose(Rbn))    # columns are branches, so a row is one lookup
    for (k, b) in enumerate(branches)
        nodes = Rt.rowval[Rt.colptr[b]:(Rt.colptr[b+1]-1)]
        signs = Rt.nzval[Rt.colptr[b]:(Rt.colptr[b+1]-1)]
        for (bi, sb) in zip(nodes, signs), (ai, sa) in zip(nodes, signs)
            push!(I, ai); push!(J, bi); push!(S, k); push!(G, sa*sb)
        end
    end
    stamp = nodalstampplan(I, J, S, G, false, n)
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

"""
    MutualStampPlan

The mutual inductance stamps: every coupling resolved once against the
graph to the two branches it couples and the sign the netlist's terminal
order gives it, so a refill reads the three values of a coupling and
nothing structural. `components` holds, for each coupling, the position of
its coefficient among the bound circuit's mutual inductors and of its two
inductances among its inductors, so the values are read from those typed
groups. `sharedbranches` lists the coupled branches on which more than one
inductor sits, which a refill refuses when the coupling is nonzero.
"""
struct MutualStampPlan{Ti<:Integer}
    stamp::NodalStampPlan{Ti}
    components::Vector{NTuple{3,Int}}
    sharedbranches::Vector{Pair{Int,Vector{Int}}}
end

"""
    mutualstampplan(c::CompiledCircuit)

Build the [`MutualStampPlan`](@ref) of a compiled circuit and its graph, and
return it with the sorted branches the couplings touch. The sign of each
coupling's stamps carries it from the terminal order the netlist declared
to the orientation the graph gave its two branches. A coupling names two
inductors, which the parse of the circuit holds it to.
"""
function mutualstampplan(c::CompiledCircuit)
    topology = c.topology
    I = Int[]; J = Int[]; S = Int[]; G = Int8[]
    from, to = isempty(c.couplings) ? (Int[], Int[]) :
        branchendpoints(topology.Rbn, topology.Nbranches)
    for (n, (_, i, j)) in enumerate(c.couplings)
        e1 = (c.nodeindices[1,i], c.nodeindices[2,i])
        e2 = (c.nodeindices[1,j], c.nodeindices[2,j])
        b1, b2 = topology.edge2indexdict[e1], topology.edge2indexdict[e2]
        negative = ((from[b1], to[b1]) == e1) != ((from[b2], to[b2]) == e2)
        sign = Int8(negative ? -1 : 1)
        append!(I, (b1, b2)); append!(J, (b2, b1))
        append!(S, (n, n)); append!(G, (sign, sign))
    end
    coupled = sort!(unique(I))
    members = Dict(b => Int[] for b in coupled)
    if !isempty(coupled)
        for i in c.inductors
            b = topology.edge2indexdict[(c.nodeindices[1,i], c.nodeindices[2,i])]
            haskey(members, b) && push!(members[b], i)
        end
    end
    shared = [b => members[b] for b in coupled if length(members[b]) > 1]
    # each coupling's positions in the typed value groups
    kpos = Dict(k => n for (n, k) in enumerate(c.mutualinductors))
    lpos = Dict(i => n for (n, i) in enumerate(c.inductors))
    components = NTuple{3,Int}[(kpos[k], lpos[i], lpos[j])
        for (k, i, j) in c.couplings]
    return MutualStampPlan(nodalstampplan(I, J, S, G, false, topology.Nbranches),
        components, shared), coupled
end

# the mutual inductance of each coupling, `K*sqrt(L1*L2)`, from the
# coupling coefficients `K` and the inductances `L` of a bound circuit
function mutualvalues!(out, plan::MutualStampPlan, K::AbstractVector,
        L::AbstractVector)
    for (n, (k, i, j)) in enumerate(plan.components)
        out[n] = K[k]*sqrt(L[i]*L[j])
    end
    return out
end

function assemblemutual(::Type{T}, plan::MutualStampPlan,
        b::BoundCircuit) where {T}
    mv = mutualvalues!(Vector{T}(undef, length(plan.components)), plan,
        b.mutualinductors, b.inductors)
    return assemblenodal(T, plan.stamp, mv, 1)
end

function checkmutualbranches(c::CompiledCircuit, plan::MutualStampPlan, Mb)
    checkcoupleddiagonal(Mb)
    for (branch, components) in plan.sharedbranches
        any(k -> !iszero(nonzeros(Mb)[k]), nzrange(Mb, branch)) || continue
        throwsharedinductors(c.componentnames[components])
    end
    return nothing
end

# === the circuit matrices ===

"""
    CircuitMatrixPlan

Everything about a circuit's matrices which depends on its topology but not
on its values, for one mode count.

Holds the nodal, branch and mutual stamp plans, the ports and the mode
expanded incidence matrix. Rebinding at new component values reuses all of
it; only the values are refilled. See [`circuitmatrixplan`](@ref) and
[`assemblematrices`](@ref).

A plan depends on the compiled circuit alone; the scratch of a refill is a
[`CircuitMatrixWorkspace`](@ref), one per solve.
"""
struct CircuitMatrixPlan{Ti<:Integer}
    circuit::CompiledCircuit
    Nmodes::Int
    capacitance::NodalStampPlan{Ti}
    conductance::NodalStampPlan{Ti}
    inductance::BranchStampPlan{Ti}
    junction::BranchStampPlan{Ti}
    invinductance::InverseInductancePlan{Ti}
    mutual::MutualStampPlan{Ti}
    ports::Vector{CompiledPort}
    Rbnm::SparseMatrixCSC{Int,Int}
end

"""
    circuitmatrixplan(c::CompiledCircuit; Nmodes = 1)

Build the [`CircuitMatrixPlan`](@ref) of a compiled circuit.
The plan depends on neither the values nor their types.

The conductance plan covers the resistors and the port owned environments
together, because an environment is realized as an ordinary resistor.
"""
@noinline function circuitmatrixplan(c::CompiledCircuit; Nmodes::Int = 1)
    topology = c.topology
    size(c.nodeindices,2) == ncomponents(c) || throw(DimensionMismatch(
        "componenttypes, nodeindices, and componentvalues should have the same length"))
    size(c.nodeindices,1) == 2 || throw(DimensionMismatch(
        "nodeindices should have a first dimension size of 2."))
    inductance = branchstampplan(c, c.inductors, topology.edge2indexdict,
        topology.Nbranches)
    mutual, coupled = mutualstampplan(c)
    return CircuitMatrixPlan(c, Nmodes,
        nodalstampplan(c, c.capacitors, c.Nnodes),
        nodalstampplan(c, c.resistors, c.Nnodes; invert = true),
        inductance,
        junctionstampplan(c, topology.edge2indexdict, topology.Nbranches),
        inverseinductanceplan(topology, inductance.nzind, coupled),
        mutual, orderedports(c), diagrepeat(topology.Rbn, Nmodes))
end

"""
    frequencydependentterms(plan::CircuitMatrixPlan,
        fd::FrequencyDependentComponents)

The [`FrequencyDependentTerms`](@ref) of the frequency dependent
components `fd` of a circuit, from the stamp plans: each contribution of
a capacitor or resistor, and each of the inverse inductance of the branch
of an inductor, which the coupled branches never are. `nothing` for
`nothing`.
"""
frequencydependentterms(plan::CircuitMatrixPlan, ::Nothing) = nothing
function frequencydependentterms(plan::CircuitMatrixPlan,
        fd::FrequencyDependentComponents)
    kind, row, col, sign, provider =
        Int8[], Int[], Int[], Int8[], Int[]
    # the terms of the sources `members`, `(source, function)` pairs, of the
    # stamp plan `stamp`: the plan's entries are gathered by source in one
    # pass, so that a component's are found without a pass of their own
    function addterms!(stamp::NodalStampPlan, k::Integer, members)
        entries = Dict{Int,Vector{Int}}(src => Int[] for (src, _) in members)
        for t in eachindex(stamp.src)
            e = get(entries, stamp.src[t], nothing)
            isnothing(e) || push!(e, t)
        end
        for (src, j) in members, t in entries[src]
            d = stamp.dest[t]
            push!(kind, k); push!(row, stamp.rowval[d])
            push!(col, searchsortedlast(stamp.colptr, d))
            push!(sign, stamp.weights[t]); push!(provider, j)
        end
    end
    addterms!(plan.capacitance, TERM_CAPACITANCE, fd.capacitors)
    addterms!(plan.conductance, TERM_CONDUCTANCE, fd.resistors)
    # an inductor's branch, and that branch's place among the retained ones
    # the inverse inductance stamps
    retained = Dict(p => r for (r, p) in enumerate(plan.invinductance.positions))
    addterms!(plan.invinductance.stamp, TERM_INVINDUCTANCE,
        [(retained[plan.inductance.dest[g]], j) for (g, j) in fd.inductors])
    return FrequencyDependentTerms(fd.providers, kind, row, col, sign,
        provider)
end

"""
    assemblematrices(plan::CircuitMatrixPlan, b::BoundCircuit)

The [`CircuitMatrices`](@ref) of a bound circuit, assembled against the fixed
patterns of `plan`.

The assembly of [`numericmatrices`](@ref) and of the solvers: every stamp,
the mutual inductances included, on the plan's patterns and resolved
couplings and ports, so that only the values and the solver scale are
read at an assembly. The numeric types are the bound circuit's.
"""
function assemblematrices(plan::CircuitMatrixPlan, b::BoundCircuit)
    c = plan.circuit
    vvn = b.values
    Nmodes = plan.Nmodes

    TC = eltype(b.capacitors)
    TR = eltype(b.resistors)
    TL = eltype(b.inductors)
    TJ = eltype(b.junctions)

    Cnm = assemblenodal(TC, plan.capacitance, b.capacitors, Nmodes)
    Gnm = assemblenodal(TR, plan.conductance, b.resistors, Nmodes)
    Lb = assemblebranch(TL, plan.inductance, b.inductors,
        combine_reciprocal_sum, 1)
    Ljb = assemblebranch(TJ, plan.junction, b.junctions, nothing, 1)
    Ljbm = Nmodes == 1 ? copy(Ljb) : diagrepeat(Ljb, Nmodes)

    TM = promote_type(eltype(b.inductors), eltype(b.mutualinductors))
    Mb = assemblemutual(TM, plan.mutual, b)
    checkmutualbranches(c, plan.mutual, Mb)
    invLnm = assembleinvinductance(TL, plan.invinductance, Lb, Nmodes)

    Lmean = inductancemean(b)
    # a port states its reference impedance and what realizes it; nothing
    # here looks at what shares its branch
    portindices = [p.component for p in plan.ports]
    portnumbers = [p.number for p in plan.ports]
    portimpedances = portreferenceimpedances(plan.ports, vvn)

    return CircuitMatrices(Cnm, Gnm, Lb, Ljb, Ljbm, Mb, invLnm,
        plan.Rbnm, portindices, portnumbers, portimpedances,
        [p.environment for p in plan.ports], Lmean, vvn,
        frequencydependentterms(plan, b.frequencydependent))
end

# === the same matrices at new values ===

# the stored values of `diagrepeat(A, Nmodes)` are those of `A` repeated
# once per mode within each column, in column order; see `diagrepeat!`
function repeatvalues!(out::AbstractVector, base::AbstractVector,
        colptr::AbstractVector, Nmodes::Integer)
    q = 1
    @inbounds for j in 1:(length(colptr) - 1), _ in 1:Nmodes,
            t in colptr[j]:(colptr[j+1] - 1)
        out[q] = base[t]
        q += 1
    end
    q == length(out) + 1 || throw(DimensionMismatch(
        "the repeated matrix does not have the pattern its base repeats."))
    return out
end

function repeatvalues!(out::AbstractVector, base::AbstractVector,
        Nmodes::Integer)
    @inbounds for i in eachindex(base), k in 1:Nmodes
        out[(i-1)*Nmodes + k] = base[i]
    end
    return out
end

# the scratch of one stamp's assembly, a solve's own rather than the
# plan's, so that solves on one plan do not write over each other
struct StampWorkspace{T}
    values::Vector{T}
    seen::Vector{Bool}
end
StampWorkspace(::Type{T}, n::Int) where {T} =
    StampWorkspace(Vector{T}(undef, n), Vector{Bool}(undef, n))

"""
    CircuitMatrixWorkspace(plan::CircuitMatrixPlan, nm::CircuitMatrices)

Scratch for refilling matrices built from `plan`. Reuse it with the same
patterns and numeric types, and keep separate workspaces for concurrent solves.
"""
struct CircuitMatrixWorkspace{TC,TR,TL,TM}
    capacitance::StampWorkspace{TC}
    conductance::StampWorkspace{TR}
    invinductance::StampWorkspace{TL}
    mutual::StampWorkspace{TM}
    inductor_seen::Vector{Bool}
    invL::Vector{TL}
    mutualvalues::Vector{TM}
end

function CircuitMatrixWorkspace(plan::CircuitMatrixPlan, nm::CircuitMatrices)
    return CircuitMatrixWorkspace(
        StampWorkspace(eltype(nm.Cnm), length(plan.capacitance.rowval)),
        StampWorkspace(eltype(nm.Gnm), length(plan.conductance.rowval)),
        StampWorkspace(eltype(nm.invLnm), length(plan.invinductance.stamp.rowval)),
        StampWorkspace(eltype(nm.Mb), length(plan.mutual.stamp.rowval)),
        Vector{Bool}(undef, length(plan.inductance.nzind)),
        Vector{eltype(nm.invLnm)}(undef, length(plan.invinductance.positions)),
        Vector{eltype(nm.Mb)}(undef, length(plan.mutual.components)))
end

function refillnodal!(A, plan::NodalStampPlan, values, work::StampWorkspace, Nmodes)
    assemblenodal!(work.values, work.seen, plan, values)
    Nmodes == 1 ? copyto!(nonzeros(A), work.values) :
        repeatvalues!(nonzeros(A), work.values, plan.colptr, Nmodes)
    return A
end

function refillbranch!(v, vm, plan::BranchStampPlan, values, seen, combine, Nmodes)
    assemblebranch!(v.nzval, seen, plan, values, combine)
    Nmodes == 1 ? copyto!(vm.nzval, v.nzval) :
        repeatvalues!(vm.nzval, v.nzval, Nmodes)
    return v
end

"""
    assemblematrices!(nm::CircuitMatrices, plan::CircuitMatrixPlan,
        b::BoundCircuit, work = CircuitMatrixWorkspace(plan, nm))

The matrices of `nm` at the values of `b`, written into the storage `nm`
already has. The patterns are the plan's and do not move, so only the
stored values are rewritten: the nodal matrices and the mutual inductance
matrix through their stamp plans into a scratch the size of one mode and
then repeated, and the branch vectors directly. Each group of `b` must
have the element type its matrix was built with, which the solver cache
checks before it refills. Pass a `CircuitMatrixWorkspace(plan, nm)` to
reuse scratch across calls; the default allocates scratch for this call.
The returned [`CircuitMatrices`](@ref) shares its matrix storage with
`nm`; value-dependent metadata is rebuilt.

This is the sweep's assembly: [`hbsolve!`](@ref) calls it at every point.
"""
function assemblematrices!(nm::CircuitMatrices, plan::CircuitMatrixPlan,
        b::BoundCircuit, work::CircuitMatrixWorkspace = CircuitMatrixWorkspace(plan, nm))
    c = plan.circuit
    vvn = b.values
    Nmodes = plan.Nmodes
    refillnodal!(nm.Cnm, plan.capacitance, b.capacitors, work.capacitance, Nmodes)
    refillnodal!(nm.Gnm, plan.conductance, b.resistors, work.conductance, Nmodes)
    assemblebranch!(nm.Lb.nzval, work.inductor_seen, plan.inductance,
        b.inductors, combine_reciprocal_sum)
    refillbranch!(nm.Ljb, nm.Ljbm, plan.junction, b.junctions, nothing,
        nothing, Nmodes)
    mutualvalues!(work.mutualvalues, plan.mutual, b.mutualinductors,
        b.inductors)
    refillnodal!(nm.Mb, plan.mutual.stamp, work.mutualvalues, work.mutual, 1)
    checkmutualbranches(c, plan.mutual, nm.Mb)
    for (k, p) in enumerate(plan.invinductance.positions)
        work.invL[k] = 1/nm.Lb.nzval[p]
    end
    refillnodal!(nm.invLnm, plan.invinductance.stamp, work.invL,
        work.invinductance, Nmodes)
    Lmean = inductancemean(b)
    return CircuitMatrices(nm.Cnm, nm.Gnm, nm.Lb, nm.Ljb, nm.Ljbm,
        nm.Mb, nm.invLnm, nm.Rbnm, nm.portindices, nm.portnumbers,
        portreferenceimpedances(plan.ports, vvn), nm.portenvironmentindices,
        Lmean, vvn, frequencydependentterms(plan, b.frequencydependent))
end

# === the scattering blocks of a compiled circuit ===

"""
    scatteringstampsystem(blocks::Vector{CompiledScatteringBlock}, Nmodes;
        auxoffset, Ntotal, scale = 1.0, modeoffsets, iscale = 1.0)

The stamp system of the compiled scattering blocks of a circuit.

A compiled block is one instance carrying its own terminal map, so the
blocks need no regrouping and no check: a block cannot be missing a port,
cannot repeat one, and cannot be confused with another instance of the
same definition. The auxiliary port currents of each block follow those
of the block before it, from `auxoffset`. `modeoffsets` are the frequency
offsets of the modes, which a pumped block's coupling between them reads,
and `iscale` the scale of the auxiliary port current unknowns; see
[`scatteringstampsystem`](@ref).
"""
function scatteringstampsystem(blocks::Vector{CompiledScatteringBlock},
    Nmodes::Integer; auxoffset::Integer, Ntotal::Integer,
    scale::Real = 1.0, modeoffsets, iscale::Real = 1.0)

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
    return scatteringstampsystem(stamped, Nmodes, Ntotal, scale;
        modeoffsets = modeoffsets, iscale = iscale)
end
