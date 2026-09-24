using JosephsonCircuits
using Test

@testset verbose=true "the circuit matrices" begin

    @testset "a coupling must name two inductors" begin
        # a junction as either member of a coupling is refused when the
        # netlist is compiled, naming the coupling and the junction
        for (first, second) in (("Lj1", "L2"), ("L1", "Lj2"))
            inductor(name, L) = startswith(name, "Lj") ? JosephsonJunction(L) : Inductor(L)
            junction = startswith(first, "Lj") ? first : second
            err = try
                circuit = Circuit([("P1", "1", "0", Port(1; Z0 = 50.0)), (first, "1", "0", inductor(first, 1.0e-9)),
                    (second, "2", "0", inductor(second, 4.0e-9)), ("K1", first, second, MutualInductor(0.1)),
                    ("C1", "2", "0", Capacitor(2.0e-12))])
                numericmatrices(circuit, Dict{Symbol,Any}()); nothing
            catch e; e end
            @test err isa ArgumentError && occursin("K1 couples $(junction)", err.msg)
        end
    end

    @testset "definitions in any dictionary" begin
        c = Circuit([(:P1, 1, 0, Port(1; Z0 = 50.0)), (:C1, 1, 0, Capacitor(:C)),
            (:L1, 1, 0, Inductor(:L))])
        ref = numericmatrices(c, Dict(:C => 1e-12, :L => 1e-9))
        other = numericmatrices(c, IdDict{Any,Any}(:C => 1e-12, :L => 1e-9))
        @test other.Cnm == ref.Cnm && other.invLnm == ref.invLnm
        # and the value table of the transient and the solver cache
        cc = compile(c)
        @test JosephsonCircuits.numericvalues(cc,
            IdDict{Any,Any}(:C => 1e-12, :L => 1e-9)) ==
            JosephsonCircuits.numericvalues(cc, Dict(:C => 1e-12, :L => 1e-9))
    end

    @testset "the sign of a mutual inductance does not follow the node names" begin
        # two inductors in series sharing a node and coupled to each other:
        # walking the chain the series inductance is L1 + L2 + 2M, whatever
        # the nodes are called.  The incidence matrix orients a branch by
        # the spanning tree, so the netlist's own order has to be carried
        # over to it.
        L = 1e-9
        K = 0.5
        Z0 = 50.0
        w = 2pi*1e9
        chain(a, b, c) = Circuit([
            (:p1, a, 0, Port(1; Z0 = Z0)),
            (:p2, c, 0, Port(2; Z0 = Z0)),
            (:l1, a, b, Inductor(L)),
            (:l2, b, c, Inductor(L)),
            (:k, :l1, :l2, MutualInductor(K))])
        function seriesinductance(nodes)
            S21 = hblinsolve([w], chain(nodes...)).S((0,), 2, (0,), 1, 1)
            return 2*Z0*imag(1/S21 - 1)/w
        end
        for nodes in ((1, 2, 3), (1, 3, 2), (2, 1, 3), (3, 2, 1), (3, 1, 2))
            @test isapprox(seriesinductance(nodes), 2*L*(1 + K); rtol = 1e-10)
        end
        # and reversing the terminals of one of them opposes the currents
        opposed = Circuit([
            (:p1, 1, 0, Port(1; Z0 = Z0)),
            (:p2, 3, 0, Port(2; Z0 = Z0)),
            (:l1, 1, 2, Inductor(L)),
            (:l2, 3, 2, Inductor(L)),
            (:k, :l1, :l2, MutualInductor(K))])
        S21 = hblinsolve([w], opposed).S((0,), 2, (0,), 1, 1)
        @test isapprox(2*Z0*imag(1/S21 - 1)/w, 2*L*(1 - K); rtol = 1e-10)
    end

    @testset "a bias current threading a loop through two couplings" begin
        # a superconducting loop from a node to ground, an inductor on one
        # arm and a junction in series with an inductor on the other, each
        # arm inductor coupled to one of two bias inductors in series across
        # a port driven with a direct current, as the flux biased SQUID of
        # the snake amplifier example. The bias current adds flux to each
        # arm inductor in the direction its terminals are declared in, so
        # with both declared toward ground the loop, which runs down one arm
        # and up the other, is threaded only through couplings of opposite
        # sign: the junction phase then solves
        # phi0*phi + (La + Lb)*Ic*sin(phi) = 2*M*I, and equal couplings put
        # the same flux on both arms and leave the junction at zero phase.
        JC = JosephsonCircuits
        La = 100e-12; Lb = 100e-12; Lj = 400e-12; Lbias = 200e-12
        K = 0.5; Idc = 3.3e-6
        Ic = JC.LjtoIc(Lj)
        M = K*sqrt(La*Lbias)
        loop(Ka, Kb, lbnodes) = Circuit([
            (:la, 1, 0, Inductor(La)),
            (:jj, 1, 2, JosephsonJunction(Lj)),
            (:lb, lbnodes..., Inductor(Lb)),
            (:p1, 3, 0, Port(1; Z0 = 1000.0)),
            (:lb1, 3, 4, Inductor(Lbias)), (:lb2, 4, 0, Inductor(Lbias)),
            (:ka, :la, :lb1, MutualInductor(Ka)),
            (:kb, :lb, :lb2, MutualInductor(Kb))])
        # the node fluxes are in units of the reduced flux quantum, so the
        # junction phase is the difference of its two
        function junctionphase(c)
            sol = hbnlsolve((2pi*1e9,), (1,),
                [(mode = (0,), port = 1, current = Idc)], c; dc = true)
            @test sol.solverinfo.converged
            return real(sol.nodeflux(outputmode = (0,), node = "1") -
                sol.nodeflux(outputmode = (0,), node = "2"))
        end
        # the phase of the loop equation, whose left side is monotone at a
        # loop inductance below the junction's
        phi = 0.0
        for _ in 1:50
            f = JC.phi0*phi + (La + Lb)*Ic*sin(phi) - 2*M*Idc
            phi -= f/(JC.phi0 + (La + Lb)*Ic*cos(phi))
        end
        @test 0.5 < phi < 2
        opposite = junctionphase(loop(K, -K, (2, 0)))
        @test isapprox(abs(opposite), phi; rtol = 1e-8)
        # reversing the terminals of an arm inductor in place of the sign
        @test isapprox(junctionphase(loop(K, K, (0, 2))), opposite;
            rtol = 1e-8)
        # reversing both signs reverses the flux
        @test isapprox(junctionphase(loop(-K, K, (2, 0))), -opposite;
            rtol = 1e-8)
        # equal couplings thread no flux through the loop
        @test isapprox(junctionphase(loop(K, K, (2, 0))), 0.0; atol = 1e-8)
    end

    @testset "the endpoints of each branch" begin
        # one walk of the incidence matrix names each branch's endpoints the
        # way the matrix itself does
        JC = JosephsonCircuits
        cc = JC.compile(Circuit([(:p1, 1, 0, Port(1)), (:p2, 3, 0, Port(2)),
            (:l1, 1, 2, Inductor(1e-9)), (:l2, 3, 2, Inductor(2e-9)),
            (:k, :l1, :l2, MutualInductor(0.3))]))
        t = cc.topology
        from, to = JC.branchendpoints(t.Rbn, t.Nbranches)
        for bi in 1:t.Nbranches
            from[bi] > 1 && @test t.Rbn[bi, from[bi]-1] == -1
            to[bi] > 1 && @test t.Rbn[bi, to[bi]-1] == 1
            @test from[bi] != to[bi]
        end
    end
end
