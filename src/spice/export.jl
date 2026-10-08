# The export of a compiled circuit as a WRSPICE netlist: the two terminal
# entries merged by branch into one table, and each branch of the table
# written as one line.

"""
    sumvalues(type::Symbol, value1, value2)

Sum together two values in different ways depending on the circuit component
type: capacitances and coupling coefficients add, and inductances and
resistances combine in parallel.

# Examples
```jldoctest
julia> JosephsonCircuits.sumvalues(:L, 1.0, 4.0)
0.8

julia> JosephsonCircuits.sumvalues(:Lj, 1.0, 4.0)
0.8

julia> JosephsonCircuits.sumvalues(:R, 1.0, 4.0)
0.8

julia> JosephsonCircuits.sumvalues(:C, 1.0, 4.0)
5.0

julia> JosephsonCircuits.sumvalues(:K, 1.0, 4.0)
5.0
```
"""
function sumvalues(type::Symbol, value1, value2)
    if type == :C || type == :K
        return value1+value2
    elseif type == :Lj || type == :L || type == :R
        return 1/(1/value1+1/value2)
    else
        error(lazy"unknown component type in sumvalues")
    end
end

# the real value a SPICE element takes: one with an imaginary part, a lossy
# element, or one which is not a number, a frequency dependent value or an
# undefined parameter, has no SPICE element
function spicevalue(value, name)
    value isa Number || throw(ArgumentError(
        lazy"the value of $(name) is $(value), which is not a number; the netlist export needs every parameter defined and no frequency dependent value."))
    iszero(imag(value)) || throw(ArgumentError(
        lazy"the value of $(name) is complex, $(value); a SPICE element takes a real value, so a lossy element has no line in the netlist."))
    return real(value)
end

"""
    spicebranches(psc::CompiledCircuit, componentvalues::AbstractVector)

The entries of a compiled circuit merged by branch, as the netlist export
writes them. Entries of one kind between the same two nodes are one
branch, whatever their order in the table, with their values combined by
[`sumvalues`](@ref); the couplings between the same two inductors are one
branch too, and an inductor which a mutual inductor couples is a branch of
its own, since the coupling names it. Ports, whose terminations are
resistors of the table, and current sources are left out.

Returns `(branches, position)`: the branches in the order they first
appear, each a named tuple of its `type`, the flat `index` of its first
entry, whose name and terminal order its line takes, and its merged
`value`; and the position of each merged branch in `branches` by
`(type, node1, node2)`, the node indices sorted, or `(:K, inductor1,
inductor2)`. Every value must be a real number (see `spicevalue`).
"""
function spicebranches(psc::CompiledCircuit, componentvalues::AbstractVector)
    length(componentvalues) == length(psc.componenttypes) || throw(DimensionMismatch(
        lazy"the circuit has $(length(psc.componenttypes)) components and $(length(componentvalues)) values."))
    coupled = Dict{Int,NTuple{2,Int}}()
    inductorscoupled = Set{Int}()
    for (k, l1, l2) in psc.couplings
        coupled[k] = (l1, l2)
        push!(inductorscoupled, l1, l2)
    end
    branches = @NamedTuple{type::Symbol, index::Int, value::Float64}[]
    position = Dict{Tuple{Symbol,Int,Int},Int}()
    for (i, type) in enumerate(psc.componenttypes)
        (type === :P || type === :I) && continue
        value = spicevalue(componentvalues[i], psc.componentnames[i])
        if type === :L && i in inductorscoupled
            push!(branches, (type = type, index = i, value = value))
            continue
        end
        key = if type === :K
            (type, coupled[i]...)
        else
            (type, minmax(psc.nodeindices[1, i], psc.nodeindices[2, i])...)
        end
        j = get(position, key, 0)
        if j == 0
            push!(branches, (type = type, index = i, value = value))
            position[key] = length(branches)
        else
            b = branches[j]
            branches[j] = (type = b.type, index = b.index,
                value = sumvalues(type, b.value, value))
        end
    end
    return branches, position
end

# the limits of the WRSPICE jj model which the junctions of a netlist
# share: the range of their critical currents relative to the model's,
# their mean, and the largest ratio of its capacitance to its critical
# current, in farads per ampere
const JJMODELICRATIOS = (0.02, 50.0)
const JJMODELMAXCJOIC = 0.99e-6

"""
    calcCjIcmean(Ic::AbstractVector, C::AbstractVector)

The critical current and the capacitance of the one WRSPICE `jj` model the
junctions of a netlist share, from the critical current `Ic` and the
shunt capacitance `C` of each junction branch: the mean `Icmean` of the
critical currents, and `Cj = CjoIc*Icmean` with `CjoIc` the smallest ratio
of a branch's capacitance to its critical current, clamped to the WRSPICE
maximum of `0.99e-6`. Every junction then takes its critical current and
the capacitance the model gives it, and the rest of its shunt capacitance
is a separate capacitor, which is why the smallest ratio is the one. A
junction without shunt capacitance, and critical currents below 0.02 or
above 50 times the mean, which the model cannot span, are refused. Returns
`(Cj, Icmean)`.

# Examples
```jldoctest
julia> Ic = JosephsonCircuits.LjtoIc.([1.0e-9, 1.1e-9]);

julia> JosephsonCircuits.calcCjIcmean(Ic, [1.0e-12, 1.2e-12])
(3.1100514965930345e-13, 3.141466158174782e-7)

julia> Ic = JosephsonCircuits.LjtoIc.([2.0e-9, 1.1e-9]);

julia> JosephsonCircuits.calcCjIcmean(Ic, [1.0e-12, 1.2e-12])
(2.2955141998662873e-13, 2.3187012119861487e-7)
```
"""
function calcCjIcmean(Ic::AbstractVector, C::AbstractVector)
    length(Ic) == length(C) || throw(DimensionMismatch(
        lazy"$(length(Ic)) critical currents and $(length(C)) capacitances."))
    Icmean = 0.0
    Icmax = 0.0
    Icmin = 0.0
    CjoIc = 0.0
    for (n, (ic, c)) in enumerate(zip(Ic, C))
        Icmean = Icmean + (ic - Icmean)/n
        ratio = c/ic
        if n == 1
            CjoIc = ratio
            Icmin = ic
        end
        Icmin = min(Icmin, ic)
        Icmax = max(Icmax, ic)
        if ratio == 0.0
            error(lazy"Cj cannot be zero in the WRSPICE JJ model.")
        end
        CjoIc = min(CjoIc, ratio)
    end

    # the range of junction sizes the jj model allows
    if Icmin/Icmean < first(JJMODELICRATIOS)
        error(lazy"Minimum junction too much smaller than average for WRSPICE.")
    end
    if Icmax/Icmean > last(JJMODELICRATIOS)
        error(lazy"Maximum junction too much larger than average for WRSPICE.")
    end

    # the largest ratio of Cj / Ic WRSPICE allows
    CjoIc = min(CjoIc, JJMODELMAXCJOIC)

    return CjoIc*Icmean, Icmean
end

# SPICE takes an element's type from the first character of its name and does
# not accept "/" in one. A legacy netlist satisfies both already, so its
# output is unchanged; a hierarchical instance path from a typed circuit
# satisfies neither, and is written with the prefix and "_" in place of "/".
function spicename(name::AbstractString, prefix::Char)
    s = replace(String(name), '/' => '_')
    if prefix == 'B' && length(s) > 2 && uppercase(s[1:2]) == "LJ"
        # a legacy junction is named Lj1 and has always been written as B1
        return string(prefix, s[3:end])
    end
    return (isempty(s) || uppercase(first(s)) != prefix) ? string(prefix, s) : s
end

# The first number free to name the phase node of a jj instance: past the
# count of the nets and past every net named by an integer, so a phase node
# never coincides with a net of the circuit.
function firstphasenode(nodenames::AbstractVector{<:AbstractString})
    n = length(nodenames) - 1
    for name in nodenames
        v = tryparse(Int, name)
        isnothing(v) || (n = max(n, v))
    end
    return n + 1
end

"""
    exportnetlist(circuit, circuitdefs::Dict; port::Int = 1, jj::Bool = true,
        vm::Real = 9.9)
    exportnetlist(psc::CompiledCircuit, circuitdefs::Dict; port::Int = 1,
        jj::Bool = true, vm::Real = 9.9)
    exportnetlist(psc::CompiledCircuit, componentvalues::AbstractVector;
        port::Int = 1, jj::Bool = true, vm::Real = 9.9)
    exportnetlist(circuit; port::Int = 1, jj::Bool = true, vm::Real = 9.9)

Export a circuit as a WRSPICE netlist. Returns a named tuple with the
netlist as a string in `netlist`, the port number in `port`, the node
count in `Nnodes`, and in `junctions` one entry per `jj` model instance
written, in the order of the netlist, with the flat component index of
the junction the instance realizes in `index` and the name of its phase
node, whose voltage WRSPICE reports as the junction phase in radians, in
`phasenode`; `portnodes` and `portcurrent` are placeholders fixed at
`1`, since the source nodes and amplitude are given directly to
[`wrspice_input_transient`](@ref) or [`wrspice_input_ac`](@ref). A fully
numeric circuit needs no `circuitdefs`; a compiled circuit whose values
are already numbers can be given those values directly as a vector in
compiled component order. `port` is the port number the sources are
applied to, recorded in the output.

The elements of one kind between the same two nodes are written as one
line, named after the first of them, whatever their order in the circuit:
capacitances add, and inductances, junctions and resistances combine in
parallel (see [`spicebranches`](@ref)). An inductor which a mutual
inductor couples keeps a line of its own, which the coupling names. A
resistor of infinite resistance is an open and writes no line. A port
writes the resistor of its termination, and a current source writes no
line: the drives are given to [`wrspice_input_transient`](@ref) or
[`wrspice_input_ac`](@ref), and [`WRspice`](@ref) writes the sources of a
transient problem itself. An ideal [`TransmissionLine`](@ref) is written as
the SPICE lossless line element with its impedance and delay; any other
scattering block has no SPICE element and is refused, as is a value which
is complex or not a number.

Component values are resolved with `circuitdefs`. With `jj = true` each
Josephson junction is written as an instance of one WRSPICE `jj` model
whose critical current is the mean over the junctions and whose
capacitance to critical current ratio is the smallest such ratio over
them, clamped to the WRSPICE maximum of `0.99e-6` (see
[`calcCjIcmean`](@ref)); the part of each junction's shunt capacitance
above what the model provides is written as a separate capacitor. The
model needs a shunt capacitance on every junction and junctions of
comparable size, and its relation is the sinusoidal one, so a
[`NonlinearInductor`](@ref) with another relation is refused. The phase
nodes of the instances are numbered past the circuit's nets, so none
coincides with a net. `vm`, in volts, is the model's product of the
critical current and the subgap resistance, which sets its subgap loss:
WRSPICE's own default, 16.5 mV, is very lossy and its range 8 to 100 mV,
which the model's `force = 1` lets `vm` exceed, so that the default 9.9 V
makes the loss small, though not zero. With `jj = false` each junction
is written as its linear inductance, which none of the model's
conditions apply to.

# Examples
```jldoctest
circuit = Circuit(
    [:P1 => Port(1; Z0 = :R),
     :C1 => Capacitor(:Cc),
     :Lj1 => JosephsonJunction(:Lj),
     :C2 => Capacitor(:Cj),
     :gnd => Ground()],
    [Net("1", [(:P1, 1), (:C1, 1)]),
     Net("2", [(:C1, 2), (:Lj1, 1), (:C2, 1)]),
     Net("0", [(:P1, 2), (:Lj1, 2), (:C2, 2), (:gnd, 1)])])

circuitdefs = Dict(
    :Lj =>1000.0e-12,
    :Cc => 100.0e-15,
    :Cj => 1000.0e-15,
    :R => 50.0)

println(JosephsonCircuits.exportnetlist(circuit, circuitdefs;port = 1, jj = true).netlist)
println("")
println(JosephsonCircuits.exportnetlist(circuit, circuitdefs;port = 1, jj = false).netlist)

# output
* SPICE Simulation
RP1_termination 1 0 50.0
C1 1 2 100.0f
B1 2 0 3 jjk ics=0.32910597847545336u
C2 2 0 674.1850813093012f
.model jjk jj(rtype=0,cct=1,icrit=0.32910597847545336u,cap=325.8149186906988f,force=1,vm=9.9)

* SPICE Simulation
RP1_termination 1 0 50.0
C1 1 2 100.0f
Lj1 2 0 1000.0000000000001p
C2 2 0 1000.0f
```
```jldoctest
circuit = Circuit(
    [:P1 => Port(1; Z0 = :R),
     :C1 => Capacitor(:Cc),
     :L1 => Inductor(:L1),
     :L2 => Inductor(:L2),
     :C2 => Capacitor(:Cj1),
     :C3 => Capacitor(:Cj2),
     :I1 => CurrentSource(:I1),
     :gnd => Ground()],
    [Net("1", [(:P1, 1), (:C1, 1)]),
     Net("2", [(:C1, 2), (:L1, 1), (:L2, 1), (:C2, 1), (:C3, 1), (:I1, 1)]),
     Net("0", [(:P1, 2), (:L1, 2), (:L2, 2), (:C2, 2), (:C3, 2), (:I1, 2),
      (:gnd, 1)])])

circuitdefs = Dict(
    :L1 =>2000.0e-12,
    :L2 =>2000.0e-12,
    :Cc => 100.0e-15,
    :Cj1 => 500.0e-15,
    :Cj2 => 500.0e-15,
    :R => 50.0,
    :I1 =>0.1)

println(JosephsonCircuits.exportnetlist(circuit, circuitdefs;port = 1, jj = true).netlist)
println("")
println(JosephsonCircuits.exportnetlist(circuit, circuitdefs;port = 1, jj = false).netlist)

# output
* SPICE Simulation
RP1_termination 1 0 50.0
C1 1 2 100.0f
L1 2 0 1000.0000000000001p
C2 2 0 1000.0f

* SPICE Simulation
RP1_termination 1 0 50.0
C1 1 2 100.0f
L1 2 0 1000.0000000000001p
C2 2 0 1000.0f
```
```jldoctest
circuit = Circuit(
    [:P1 => Port(1; Z0 = :Rleft),
     :L1 => Inductor(:L1),
     :Lj1 => JosephsonJunction(:Lj1),
     :L2 => Inductor(:L2),
     :K1 => MutualInductor(:K1, :L1, :L2),
     :C2 => Capacitor(:C2),
     :C3 => Capacitor(:C3),
     :gnd => Ground()],
    [Net("1", [(:P1, 1), (:L1, 1)]),
     Net("2", [(:Lj1, 1), (:L2, 1), (:C2, 1), (:C3, 1)]),
     Net("0", [(:P1, 2), (:L1, 2), (:Lj1, 2), (:L2, 2), (:C2, 2), (:C3, 2),
      (:gnd, 1)])])
circuitdefs = Dict(
    :Rleft => 50.0,
    :L1 => 1000.0e-12,
    :Lj1 => 1000.0e-12,
    :K1 => 0.1,
    :L2 => 1000.0e-12,
    :C2 => 1000.0e-15,
    :C3 => 1000.0e-15)

println(JosephsonCircuits.exportnetlist(circuit, circuitdefs;port = 1, jj = true).netlist)
println("")
println(JosephsonCircuits.exportnetlist(circuit, circuitdefs;port = 1, jj = false).netlist)

# output
* SPICE Simulation
RP1_termination 1 0 50.0
L1 1 0 1000.0000000000001p
B1 2 0 3 jjk ics=0.32910597847545336u
C2 2 0 1674.185081309301f
L2 2 0 1000.0000000000001p
K1 L1 L2 0.1
.model jjk jj(rtype=0,cct=1,icrit=0.32910597847545336u,cap=325.8149186906988f,force=1,vm=9.9)

* SPICE Simulation
RP1_termination 1 0 50.0
L1 1 0 1000.0000000000001p
Lj1 2 0 1000.0000000000001p
L2 2 0 1000.0000000000001p
K1 L1 L2 0.1
C2 2 0 2000.0f
```
```jldoctest
circuit = Circuit(
    [:P1 => Port(1; Z0 = :Rleft),
     :L1 => Inductor(:L1),
     :Lj1 => JosephsonJunction(:Lj1),
     :L2 => Inductor(:L2),
     :K1 => MutualInductor(:K1, :L2, :L1),
     :C2 => Capacitor(:C2),
     :C3 => Capacitor(:C3),
     :gnd => Ground()],
    [Net("1", [(:P1, 1), (:L1, 1)]),
     Net("2", [(:Lj1, 1), (:L2, 1), (:C2, 1), (:C3, 1)]),
     Net("0", [(:P1, 2), (:L1, 2), (:Lj1, 2), (:L2, 2), (:C2, 2), (:C3, 2),
      (:gnd, 1)])])
circuitdefs = Dict(
    :Rleft => 50.0,
    :L1 => 1000.0e-12,
    :Lj1 => 1000.0e-12,
    :K1 => 0.1,
    :L2 => 1000.0e-12,
    :C2 => 1000.0e-15,
    :C3 => 1000.0e-15)

println(JosephsonCircuits.exportnetlist(circuit, circuitdefs;port = 1, jj = true).netlist)
println("")
println(JosephsonCircuits.exportnetlist(circuit, circuitdefs;port = 1, jj = false).netlist)

# output
* SPICE Simulation
RP1_termination 1 0 50.0
L1 1 0 1000.0000000000001p
B1 2 0 3 jjk ics=0.32910597847545336u
C2 2 0 1674.185081309301f
L2 2 0 1000.0000000000001p
K1 L2 L1 0.1
.model jjk jj(rtype=0,cct=1,icrit=0.32910597847545336u,cap=325.8149186906988f,force=1,vm=9.9)

* SPICE Simulation
RP1_termination 1 0 50.0
L1 1 0 1000.0000000000001p
Lj1 2 0 1000.0000000000001p
L2 2 0 1000.0000000000001p
K1 L2 L1 0.1
C2 2 0 2000.0f
```
"""
function exportnetlist(circuit::CompilableCircuit, circuitdefs::Dict;
        port::Int = 1, jj::Bool = true, vm::Real = 9.9)
    return exportnetlist(compile(circuit), circuitdefs; port = port, jj = jj,
        vm = vm)
end

function exportnetlist(psc::CompiledCircuit,circuitdefs::Dict;
        port::Int = 1, jj::Bool = true, vm::Real = 9.9)
    return exportnetlist(psc,
        componentvaluestonumber(psc.componentvalues,circuitdefs);
        port = port, jj = jj, vm = vm)
end

function exportnetlist(psc::CompiledCircuit,componentvalues::AbstractVector;
        port::Int = 1, jj::Bool = true, vm::Real = 9.9)

    # an ideal lossless transmission line is the one scattering block
    # with a SPICE element; any other block has none, and exporting the
    # circuit without it would simulate a different circuit
    tlineblocks = [b for b in psc.scatteringblocks
        if b.definition.provider isa TransmissionLineProvider]
    if length(tlineblocks) != length(psc.scatteringblocks)
        others = [b.path for b in psc.scatteringblocks
            if !(b.definition.provider isa TransmissionLineProvider)]
        blocknames = join(others, ", ")
        throw(ComponentNotSupportedError(
            lazy"the circuit has $(length(others)) scattering block(s) ($(blocknames)), which the netlist export cannot express; export lumped elements and transmission lines only."))
    end

    # the jj model is the sinusoidal junction
    if jj && !isempty(psc.junctioncprs)
        path = psc.componentnames[minimum(keys(psc.junctioncprs))]
        throw(ComponentNotSupportedError(
            lazy"the NonlinearInductor at $(path) has a current-phase relation other than the sinusoidal one, which the WRSPICE jj model cannot express; export with jj = false to write its linear inductance."))
    end

    # placeholders; only a single port is handled
    portnodes = 1
    portcurrent = 1

    Nnodes = length(psc.nodenames)
    componentnames = psc.componentnames
    nodeindexarray = psc.nodeindices
    uniquenodevector = psc.nodenames
    branches, position = spicebranches(psc, componentvalues)
    couplingof = Dict(k => (l1, l2) for (k, l1, l2) in psc.couplings)

    # the capacitance on the branch of each junction, which the jj model
    # takes its share of, wherever the capacitors sit in the table
    function shunt(b)
        n1, n2 = minmax(nodeindexarray[1, b.index], nodeindexarray[2, b.index])
        j = get(position, (:C, n1, n2), 0)
        return j == 0 ? (0.0, 0) : (branches[j].value, branches[j].index)
    end

    # the jj model shared by the junctions
    CjoIc = 0.0
    Icmean = 0.0
    if jj && any(b -> b.type === :Lj, branches)
        junctionbranches = [b for b in branches if b.type === :Lj]
        Cj, Icmean = calcCjIcmean(
            [LjtoIc(b.value) for b in junctionbranches],
            [first(shunt(b)) for b in junctionbranches])
        CjoIc = Cj/Icmean
    end

    # multiply by these scale factors for the prefixes
    femto = 1e15
    pico = 1e12
    micro = 1e6

    # `vm`, (reference icrit)*rsub, exceeds WRSPICE's range of 8e-3 to
    # 100e-3 only with `force=1` in the model's arguments; `force=0` does
    # not turn the flag off, leaving it out does.
    # http://www.wrcad.com/ftp/pub/jj.va

    # define an array of strings for the netlist
    netlist =  ["* SPICE Simulation"]

    # one entry per jj model instance: the flat component index of the
    # junction it realizes and the name of its phase node
    junctions = @NamedTuple{index::Int, phasenode::String}[]
    phasenode = firstphasenode(uniquenodevector)

    for b in branches
        i = b.index
        value = b.value
        if b.type == :K
            # the coupled inductors by the names their own lines carry
            l1, l2 = couplingof[i]
            push!(netlist,"$(spicename(componentnames[i],'K')) $(spicename(componentnames[l1],'L')) $(spicename(componentnames[l2],'L')) $(value)")
            continue
        end
        node1 = uniquenodevector[nodeindexarray[1, i]]
        node2 = uniquenodevector[nodeindexarray[2, i]]
        if b.type == :Lj && jj
            Ictmp = LjtoIc(value)
            push!(netlist,"$(spicename(componentnames[i],'B')) $(node1) $(node2) $(phasenode) jjk ics=$(LjtoIc(value)*micro)u")
            push!(junctions,(index = i, phasenode = string(phasenode)))
            phasenode += 1

            # the shunt capacitance beyond the model's
            capvalue, capindex = shunt(b)
            if capvalue > Ictmp*CjoIc
                push!(netlist,"$(spicename(componentnames[capindex],'C')) $(node1) $(node2) $(femto*(capvalue-Ictmp*CjoIc))f")
            end
        elseif b.type == :Lj || b.type == :L
            push!(netlist,"$(spicename(componentnames[i],'L')) $(node1) $(node2) $(value*pico)p")
        elseif b.type == :C
            # a junction's shunt capacitance is written with the junction
            n1, n2 = minmax(nodeindexarray[1, i], nodeindexarray[2, i])
            jj && haskey(position, (:Lj, n1, n2)) && continue
            push!(netlist,"$(spicename(componentnames[i],'C')) $(node1) $(node2) $(value*femto)f")
        elseif b.type == :R
            # every resistor, the environment a port owns included. The port
            # itself writes no line, so its source impedance reaches the
            # exported circuit only through this one; dropping it as a
            # lowering artifact would export a different circuit. An
            # infinite resistance is an open, which is no element at all.
            isfinite(value) || continue
            push!(netlist,"$(spicename(componentnames[i],'R')) $(node1) $(node2) $(value)")
        end
    end

    # each transmission line as the lossless line element, its two ports
    # between their signal and reference terminals
    for b in tlineblocks
        provider = b.definition.provider
        (isfinite(provider.Z0) && provider.Z0 > 0 && isfinite(provider.delay) && provider.delay > 0) || throw(ArgumentError(
            lazy"the transmission line at $(b.path) needs a finite positive impedance and a positive delay to be written as a SPICE element."))
        node = i -> uniquenodevector[i]
        push!(netlist,"$(spicename(b.path,'T')) $(node(b.signalnodes[1])) $(node(b.refnodes[1])) $(node(b.signalnodes[2])) $(node(b.refnodes[2])) z0=$(provider.Z0) td=$(provider.delay)")
    end

    if !isempty(junctions)
        push!(netlist,".model jjk jj(rtype=0,cct=1,icrit=$(micro*Icmean)u,cap=$(femto*Icmean*CjoIc)f,force=1,vm=$(vm))")
    end

    return  (netlist=join(netlist,"\n"),portnodes=portnodes,port=port,portcurrent=portcurrent,Nnodes = Nnodes,junctions=junctions)
end

# A fully numeric circuit (such as the output of a circuit builder) needs no
# component definitions, so `circuitdefs` is optional.
function exportnetlist(circuit; kwargs...)
    return exportnetlist(circuit, Dict(); kwargs...)
end
