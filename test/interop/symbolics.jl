# The Symbolics extension: component values written as `Num`, the
# deprecated symbolic frequency variable spelled as one, the symbolic
# matrices, and the design sensitivities of such values.
using Symbolics

@testset verbose=true "Symbolics" begin

@variables Rleft Lj Cc Cj wsym

sym_ws = 2*pi*(4.5:0.1:5.0)*1e9
sym_wp = (2*pi*4.75001e9,)
sym_src = [(mode = (1,), port = 1, current = 0.00565e-6)]

symjpa(z, c1, lj, c2) = Circuit([
    (:P1, 1, 0, Port(1; Z0 = z)), (:C1, 1, 2, Capacitor(c1)),
    (:Lj1, 2, 0, JosephsonJunction(lj)), (:C2, 2, 0, Capacitor(c2))])

@testset "the extension loads" begin
    @test !isnothing(Base.get_extension(JosephsonCircuits,
        :JosephsonCircuitsSymbolicsExt))
end

@testset "a value written as a Num" begin
    defs = Dict(Rleft => 50.0, Lj => 1000.0e-12, Cc => 100.0e-15,
        Cj => 1000.0e-15)
    Ssym = hbsolve(sym_ws, sym_wp, sym_src, (2,), (8,),
        symjpa(Rleft, Cc, Lj, Cj), defs; keyedarrays = false).linearized.S
    Snum = hbsolve(sym_ws, sym_wp, sym_src, (2,), (8,),
        symjpa(50.0, 100.0e-15, 1000.0e-12, 1000.0e-15);
        keyedarrays = false).linearized.S
    @test isapprox(Array(Ssym), Array(Snum), rtol = 1e-12)
end

@testset "the symbolic matrices carry the expressions" begin
    # the matrices of a circuit whose values are symbols, substituted
    # afterwards, are the matrices of the circuit at those values
    m = symbolicmatrices(symjpa(Rleft, Cc, Lj, Cj))
    defs = Dict(Rleft => 50.0, Lj => 1000.0e-12, Cc => 100.0e-15,
        Cj => 1000.0e-15)
    n = JosephsonCircuits.numericmatrices(symjpa(Rleft, Cc, Lj, Cj), defs)
    at(v) = Float64(Symbolics.value(Symbolics.substitute(v, defs)))
    @test m.Cnm.colptr == n.Cnm.colptr && m.Cnm.rowval == n.Cnm.rowval
    @test at.(m.Cnm.nzval) ≈ n.Cnm.nzval
    @test at.(m.invLnm.nzval) ≈ n.invLnm.nzval
    @test at(m.Lmean) ≈ n.Lmean
end

@testset "an undefined Num parameter is named" begin
    @variables Lundef
    err = try
        hbsolve(sym_ws, sym_wp, sym_src, (2,), (8,),
            symjpa(50.0, 100.0e-15, Lundef, 1000.0e-15), Dict())
        nothing
    catch e; e; end
    @test err isa ArgumentError
    @test occursin("Lj1", sprint(showerror, err))
    @test occursin("Lundef", sprint(showerror, err))
end

@testset "the deprecated symbolic frequency variable" begin
    law(x) = 50.0*(1 + (x/1e11)^2)
    Ssym = (@test_logs (:warn,) match_mode = :any hbsolve(sym_ws, sym_wp,
        sym_src, (2,), (8,), symjpa(law(wsym), 100.0e-15, 1000.0e-12,
        1000.0e-15), Dict(); symfreqvar = wsym,
        keyedarrays = false)).linearized.S
    Sfun = hbsolve(sym_ws, sym_wp, sym_src, (2,), (8,),
        symjpa(FrequencyDependent(law), 100.0e-15, 1000.0e-12, 1000.0e-15);
        keyedarrays = false).linearized.S
    @test isapprox(Array(Ssym), Array(Sfun), rtol = 1e-12)
    # and a value named by a symbol whose definition is written in it
    Snamed = (@test_logs (:warn,) match_mode = :any hbsolve(sym_ws, sym_wp,
        sym_src, (2,), (8,), symjpa(:Z, 100.0e-15, 1000.0e-12, 1000.0e-15),
        Dict(:Z => law(wsym)); symfreqvar = wsym,
        keyedarrays = false)).linearized.S
    @test isapprox(Array(Snamed), Array(Sfun), rtol = 1e-12)
    # every entry point which takes it warns
    c = symjpa(law(wsym), 100.0e-15, 1000.0e-12, 1000.0e-15)
    @test_logs (:warn,) match_mode = :any hbnlsolve(sym_wp, (2,), sym_src,
        c, Dict(); symfreqvar = wsym)
    @test_logs (:warn,) match_mode = :any hblinsolve(sym_ws, c, Dict();
        symfreqvar = wsym)
    # the rewrite leaves a value depending on something else alone, so the
    # undefined variable is still reported with its component
    @variables Rundef
    c = Circuit([("P1", "1", "0", Port(1; Z0 = Rundef/(1 + wsym*1e-12))),
        ("C1", "1", "0", Capacitor(100.0e-15)),
        ("L1", "1", "0", Inductor(1.0e-9))])
    err = @test_logs (:warn,) match_mode = :any try
        hbnlsolve((2*pi*5e9,), (1,), [(mode=(1,),port=1,current=1e-8)], c,
            Dict(); symfreqvar = wsym)
        nothing
    catch e; e; end
    @test err isa ArgumentError
    @test occursin("Rundef", sprint(showerror, err))
    @test occursin("P1", sprint(showerror, err))
    # a table of values already resolved does not carry the definitions
    # the keyword needs, and is refused
    psc = compile(symjpa(50.0, 100.0e-15, 1000.0e-12, 1000.0e-15))
    sf = JosephsonCircuits.truncfreqs(JosephsonCircuits.calcfreqsdft((2,));
        dc = true, odd = false, even = true, maxintermodorder = Inf)
    vvn = JosephsonCircuits.componentvaluestonumber(psc.componentvalues,
        Dict{Any,Any}())
    @test_throws ArgumentError hblinsolve(sym_ws, psc, vvn, sf;
        symfreqvar = wsym)
end

@testset "design sensitivities of Num values" begin
    # a component written as a Num variable is written as the name of its
    # symbol, which the values function gives
    c = symjpa(50.0, Cc, Lj, 1000.0e-15)
    p = (Lj = 1000.0e-12, Cc = 100.0e-15)
    # the values are the parameters, whose derivatives are stated
    r = designsensitivities(c, identity, p, sym_ws, sym_wp, sym_src, (2,), (8,);
        jacobian = p -> [1.0 0.0; 0.0 1.0])
    @test Set(collect(JosephsonCircuits.AxisKeys.axiskeys(r.dSdp,
        :parameter))) == Set([:Lj, :Cc])
    # against centered finite differences of the whole solve
    for q in (:Lj, :Cc)
        h = 1e-6*p[q]
        at(v) = Dict{Any,Any}(pairs(merge(p, NamedTuple{(q,)}((v,)))))
        Sp = hbsolve(sym_ws, sym_wp, sym_src, (2,), (8,), c,
            at(p[q] + h)).linearized.S
        Sm = hbsolve(sym_ws, sym_wp, sym_src, (2,), (8,), c,
            at(p[q] - h)).linearized.S
        fd = vec((Array(Sp) .- Array(Sm))./(2*h))
        mine = vec(Array(r.dSdp(parameter = q)))
        @test norm(mine .- fd)/norm(fd) < 1e-4
    end
end

@testset "an unwrapped symbolic value is the wrapped one" begin
    # a `Num` is a wrapper; the unwrapped `BasicSymbolic` reaches the
    # package whenever a caller takes a value apart, and must be read the
    # same way
    unwrap(v) = Symbolics.value(v)
    defs = Dict(Lj => 1000.0e-12, Cc => 100.0e-15)
    wrapped = symjpa(50.0, Cc, Lj, 1000.0e-15)
    raw = symjpa(50.0, unwrap(Cc), unwrap(Lj), 1000.0e-15)
    @test isapprox(Array(hblinsolve(sym_ws, raw, defs; keyedarrays = false).S),
        Array(hblinsolve(sym_ws, wrapped, defs; keyedarrays = false).S),
        rtol = 1e-12)
    p = (Lj = 1000.0e-12, Cc = 100.0e-15)
    J(p) = [1.0 0.0; 0.0 1.0]
    @test JosephsonCircuits.designjacobian(raw, identity, p; jacobian = J) ==
        JosephsonCircuits.designjacobian(wrapped, identity, p; jacobian = J)
end

@testset "a value keeps its undefined parameters" begin
    # the definitions are substituted as they are given, so a value which
    # only some of them define stays symbolic and resolves at the rest
    partial = JosephsonCircuits.valuetonumber(Lj + Cc,
        Dict{Any,Any}(Lj => 1000.0e-12))
    @test JosephsonCircuits.checkissymbolic(partial)
    @test JosephsonCircuits.valuetonumber(partial,
        Dict{Any,Any}(Cc => 100.0e-15)) ≈ 1000.0e-12 + 100.0e-15
end

@testset "a cache over definitions keyed by Num" begin
    defs = Dict(Lj => 1000.0e-12, Cc => 100.0e-15)
    c = symjpa(50.0, Cc, Lj, 1000.0e-15)
    cache = hbcache(sym_wp, (4,), sym_src, c, defs; atol = 1e-12)
    moved = hbsolve!(cache, (Lj = 900.0e-12,))
    @test cache.converged
    fresh = JosephsonCircuits.hbnlsolve(sym_wp, (4,), sym_src, c,
        Dict(Lj => 900.0e-12, Cc => 100.0e-15); atol = 1e-12,
        keyedarrays = false)
    @test isapprox(vec(collect(moved.nodeflux)), vec(collect(fresh.nodeflux));
        rtol = 1e-8)
    # a parameter defined under its Num, its symbol and its string moves
    # under every one of them
    aliases = Dict{Any,Any}(Lj => 1000.0e-12, :Lj => 1000.0e-12,
        "Lj" => 1000.0e-12)
    at = JosephsonCircuits.definitionsat(aliases,
        JosephsonCircuits.definitionkeys(aliases), (Lj = 900.0e-12,))
    @test at[Lj] == at[:Lj] == at["Lj"] == 900.0e-12
    @test aliases[Lj] == 1000.0e-12
end

@testset "a Num value defined under any of its names" begin
    for key in (Lj, :Lj, "Lj")
        @test JosephsonCircuits.valuetonumber(2*Lj,
            Dict{Any,Any}(key => 1000.0e-12)) ≈ 2000.0e-12
    end
end

@testset "a Num variable defined by a closure" begin
    # a bare variable is a name, and defined as a frequency dependent value
    # takes that value, as a symbol does
    ws = 2*pi*[4.5e9, 5.0e9]
    S(c, d) = hblinsolve(ws, c, d; keyedarrays = false).S
    @test S(symjpa(50.0, Cc, 1000.0e-12, 1000.0e-15),
        Dict(Cc => FrequencyDependent(w -> 100.0e-15))) ≈
        S(symjpa(50.0, 100.0e-15, 1000.0e-12, 1000.0e-15), Dict())
end

@testset "the symbolic matrices of a circuit written with names" begin
    # a name is the variable of that name, so the matrices are those of the
    # circuit written with the variables
    named = symbolicmatrices(symjpa(:Rleft, :Cc, :Lj, :Cj))
    vars = symbolicmatrices(symjpa(Rleft, Cc, Lj, Cj))
    for f in (:Cnm, :Gnm, :invLnm, :Ljb)
        @test isequal(getfield(named, f), getfield(vars, f))
    end
end

@testset "the examples of symbolicmatrices" begin
    # the examples of its docstring and of the circuits page
    # (docs/src/circuits.md) as the docs write them, which the docs build
    # cannot run without Symbolics, and the values they show
    @variables R
    circuit = Circuit(
        [:p1 => Port(1; Z0 = :Rleft),
         :i1 => CurrentSource(:Ipump),
         :cc => Capacitor(:Cc),
         :jj => JosephsonJunction(:Lj),
         :cj => Capacitor(:Cj),
         :gnd => Ground()],
        [[(:p1, 1), (:i1, 1), (:cc, 1)],
         [(:cc, 2), (:jj, 1), (:cj, 1)],
         [(:p1, 2), (:i1, 2), (:jj, 2), (:cj, 2), (:gnd, 1)]])
    m = symbolicmatrices(circuit)
    @test isequal((m.Cnm[1, 2], m.Cnm[2, 2]), (-Cc, Cc + Cj))
    @test isequal((m.Gnm[1, 1], m.Ljb[2]), (1 / Rleft, Lj))
    jpa = Circuit([
        (:p1, 1, 0, Port(1; Z0 = :R)),
        (:cc, 1, 2, Capacitor(:Cc)),
        (:jj, 2, 0, JosephsonJunction(:Lj)),
        (:cj, 2, 0, Capacitor(:Cj)),
    ])
    m = symbolicmatrices(jpa)
    @test isequal(Matrix(m.Cnm), [Cc -Cc; -Cc Cc + Cj])
    @test isequal((m.Gnm[1, 1], m.Ljb[2]), (1 / R, Lj))
end

@testset "a netlist file holds a Num expression" begin
    # a value written as an expression is read back by Symbolics
    filename = tempname()
    netlist = [("P","1","0",1), ("R","1","0",2*Rleft)]
    @test_logs (:warn,) match_mode = :any JosephsonCircuits.export_netlist(
        filename, netlist)
    back = @test_logs (:warn,) match_mode = :any JosephsonCircuits.import_netlist(
        filename)
    rm(filename)
    @test isequal(back, netlist)
end

@testset "the static stiffness of a symbolic inductance" begin
    # a symbolic value is taken to provide finite, nonzero stiffness
    @test isnothing(JosephsonCircuits.checkstaticstiffnessvalues(
        [:P,:R,:L,:L], Any[1, 50.0, 1e-9, Lj]))
    @test JosephsonCircuits.calcstaticfluxcomponents(
        [:P,:R,:L],[2 2 2;1 1 1],Any[1,50.0,Lj],2) == Vector{Int64}[]
end

@testset "the value handling methods" begin
    defs = Dict{Any,Any}(Lj => 1000.0e-12, Cc => 100.0e-15)
    @test JosephsonCircuits.valuetonumber(2*Lj, defs) ≈ 2000.0e-12
    @test JosephsonCircuits.checkissymbolic(Lj)
    @test !JosephsonCircuits.checkissymbolic(Symbolics.value(Num(2.0)))
    @test Set(Symbolics.tosymbol.(JosephsonCircuits.circuitvariables(Lj*Cc);
        escape = false)) == Set([:Lj, :Cc])
    @test JosephsonCircuits.definitionname(Lj) === :Lj
    # a symbolic value states no frequency of its own
    @test JosephsonCircuits.substitutefreq(Num(3.0), 1e9) == 3.0
    # the deprecated netlist file writes a value resolved as a solve does
    @test JosephsonCircuits.substitutedefs(2*Lj, Dict(:Lj => 1e-9)) ≈ 2e-9
end

@testset "a symbolic resistor across a port which owns its environment" begin
    # The port owns a matched environment and a resistor of its own sits
    # across the same terminals: numerically that is the duplicate load the
    # compiler warns about, but with the resistor symbolic the comparison is
    # a symbolic value and not a Bool, and the circuit must compile, carry
    # the expression into its matrices, and solve at a definition.
    @variables R
    circuit(v) = Circuit([(:p, 1, 0, Port(1)), (:r, 1, 0, Resistor(v)),
        (:c, 1, 0, Capacitor(1e-12))])
    cc = @test_logs compile(circuit(R))
    @test cc.componentvalues[cc.componentnamedict["r"]] === R
    sm = @test_logs symbolicmatrices(circuit(R))
    @test any(isequal(R), Symbolics.get_variables(sum(sm.Gnm)))
    sol = hblinsolve([2*pi*1e9], circuit(R), Dict(R => 75.0))
    @test size(sol.S) == size(hblinsolve([2*pi*1e9], circuit(75.0)).S)
    @test Array(sol.S) ≈ Array(hblinsolve([2*pi*1e9], circuit(75.0)).S)
    # and the numeric circuit it stands for still reports the double load
    @test_logs (:warn,) match_mode = :any compile(circuit(50.0))
end

@testset "a Num port number of a deprecated tuple netlist" begin
    c = @test_logs (:warn,) match_mode = :any Circuit([
        ("P1", "1", "0", Num(1)), ("R1", "1", "0", 50.0),
        ("C1", "1", "0", 1.0e-12)])
    @test compile(c).ports[1].number == 1
end

end
