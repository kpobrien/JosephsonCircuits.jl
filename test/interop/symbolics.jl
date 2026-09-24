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
end

@testset "design sensitivities of Num values" begin
    c = symjpa(50.0, Cc, Lj, 1000.0e-15)
    defs = Dict(Lj => 1000.0e-12, Cc => 100.0e-15)
    r = designsensitivities(c, defs, sym_ws, sym_wp, sym_src, (2,), (8,))
    @test Set(collect(JosephsonCircuits.AxisKeys.axiskeys(r.dSdp,
        :parameter))) == Set([:Lj, :Cc])
    # against centered finite differences of the whole solve
    for q in (Lj, Cc)
        h = 1e-6*defs[q]
        at(v) = merge(Dict{Any,Any}(defs), Dict{Any,Any}(q => v))
        Sp = hbsolve(sym_ws, sym_wp, sym_src, (2,), (8,), c,
            at(defs[q] + h)).linearized.S
        Sm = hbsolve(sym_ws, sym_wp, sym_src, (2,), (8,), c,
            at(defs[q] - h)).linearized.S
        fd = vec((Array(Sp) .- Array(Sm))./(2*h))
        mine = vec(Array(r.dSdp(
            parameter = Symbolics.tosymbol(q; escape = false))))
        @test norm(mine .- fd)/norm(fd) < 1e-4
    end
    # a parameter the value does not depend on moves nothing
    @test JosephsonCircuits.designderivative(Symbolics.value(Cc), :Lj,
        Dict{Any,Any}(Cc => 100.0e-15)) == 0
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
    @test JosephsonCircuits.designjacobian(raw, defs) ==
        JosephsonCircuits.designjacobian(wrapped, defs)
    @test JosephsonCircuits.designderivative(unwrap(Lj*Cc + Cc^2), :Lj, defs) ≈
        100.0e-15
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
    for key in (Lj, :Lj, "Lj", JosephsonCircuits.CircuitValues.Parameter(:Lj))
        @test JosephsonCircuits.valuetonumber(2*Lj,
            Dict{Any,Any}(key => 1000.0e-12)) ≈ 2000.0e-12
    end
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
