using JosephsonCircuits
using LinearAlgebra
using Test

# Direct current through resistors. The harmonic balance state is periodic
# node flux, so a voltage is its time derivative and vanishes at zero
# frequency: a resistor was an open circuit at DC and a current source
# driving one could not develop I*R. The missing coordinate is the average
# voltage, which is constant on each static flux component because an
# inductor or a zero-voltage junction is a short there.

@testset verbose=true "direct current through resistors" begin

    ws = (2*pi*5e9,)
    dcsolve(c, sources) = hbnlsolve(ws, (1,), sources, c, Dict{Any,Any}();
        keyedarrays = false, dc = true, odd = true)

    @testset "a resistor obeys Ohm's law" begin
        # no inductor or junction on the driven node: either is a short at
        # DC and would hold it at zero, which is correct physics but would
        # hide the behavior under test. A capacitor is an open circuit at DC
        # and does not disturb the resistive path.
        R = 50.0; Idc = 1.0e-6
        c = Circuit([:p1 => Port(1; Z0 = R), :c1 => Capacitor(1.0e-12)],
            [[(:p1,1),(:c1,1)], [(:p1,2),(:c1,2), Ground]])
        for sgn in (1, -1)
            sol = dcsolve(c, [(mode=(0,), port=1, current=sgn*Idc)])
            @test sol.solverinfo.converged
            # ground first and identically zero
            @test sol.dcnodevoltage[1] == 0
            @test isapprox(sol.dcnodevoltage[2], sgn*Idc*R; rtol = 1e-9)
        end
    end

    @testset "a resistive divider divides" begin
        # the port environment in parallel with the series pair
        Rp, Rs, Rl, Idc = 50.0, 50.0, 150.0, 1.0e-6
        c = Circuit(
            [:p1 => Port(1; Z0 = Rp), :rs => Resistor(Rs), :rl => Resistor(Rl),
             :c1 => Capacitor(1.0e-12)],
            [[(:p1,1),(:rs,1),(:c1,1)], [(:rs,2),(:rl,1)],
             [(:p1,2),(:rl,2),(:c1,2), Ground]])
        sol = dcsolve(c, [(mode=(0,), port=1, current=Idc)])
        @test sol.solverinfo.converged
        Vtop = Idc/(1/Rp + 1/(Rs+Rl))
        v = sort(sol.dcnodevoltage)
        @test isapprox(v[3], Vtop; rtol = 1e-9)
        @test isapprox(v[2], Vtop*Rl/(Rs+Rl); rtol = 1e-9)
    end

    @testset "an inductor shorts the resistor across it" begin
        # at DC the inductor holds both of its nodes at the same average
        # voltage, and reaching ground through one fixes it at zero, so the
        # resistor carries no direct current however hard it is driven
        c = Circuit(
            [:p1 => Port(1; Z0 = 50.0), :l1 => Inductor(1.0e-9),
             :c1 => Capacitor(1.0e-12)],
            [[(:p1,1),(:l1,1),(:c1,1)], [(:p1,2),(:l1,2),(:c1,2), Ground]])
        sol = dcsolve(c, [(mode=(0,), port=1, current=1.0e-6)])
        @test sol.solverinfo.converged
        @test isnothing(sol.dcnodevoltage) || all(iszero, sol.dcnodevoltage)
    end

    @testset "resistors carry current between floating components" begin
        # two floating nodes, each a capacitor to ground, joined by a
        # resistor and driven in and out. The bridge carries it and develops
        # I*R; the port environments are large enough to divert a part in a
        # ten million, which the tolerance allows for.
        Rb, Rbig, Idc = 100.0, 1.0e9, 1.0e-6
        c = Circuit(
            [:p1 => Port(1; Z0 = Rbig), :p2 => Port(2; Z0 = Rbig),
             :rb => Resistor(Rb), :c1 => Capacitor(1e-12),
             :c2 => Capacitor(1e-12)],
            [[(:p1,1),(:rb,1),(:c1,1)], [(:rb,2),(:p2,1),(:c2,1)],
             [(:p1,2),(:p2,2),(:c1,2),(:c2,2), Ground]])
        sol = dcsolve(c, [(mode=(0,), port=1, current=Idc),
                          (mode=(0,), port=2, current=-Idc)])
        @test sol.solverinfo.converged
        v = sol.dcnodevoltage
        @test isapprox(maximum(v) - minimum(v), Idc*Rb; rtol = 1e-5)
    end

    @testset "no direct current, no voltage" begin
        # the common case: the elimination costs one projection and reports
        # nothing, and the answer is what it always was
        c = Circuit([:p1 => Port(1), :cc => Capacitor(100e-15),
                     :jj => JosephsonJunction(1000e-12),
                     :cj => Capacitor(1000e-15)],
            [[(:p1,1),(:cc,1)], [(:cc,2),(:jj,1),(:cj,1)],
             [(:p1,2),(:jj,2),(:cj,2), Ground]])
        sol = hbnlsolve((2*pi*4.75001e9,), (8,),
            [(mode=(1,), port=1, current=0.00565e-6)], c, Dict{Any,Any}();
            keyedarrays = false)
        @test sol.solverinfo.converged
        @test isnothing(sol.dcnodevoltage)
    end

    @testset "a floating island fixes only its differences" begin
        # A differential port has both terminals off ground, so the
        # environment it owns sits between them and the pair has no
        # conductance to ground at all. Such an island has no absolute
        # voltage: the solve pins one component and reports differences,
        # which are what is physical. The earlier cases all reach ground
        # through a port environment and so never take this path.
        Z0, Idc = 200.0, 1.0e-6
        c = Circuit(
            [:p1 => Port(1; Z0 = Z0), :ca => Capacitor(1e-12),
             :cb => Capacitor(1e-12)],
            [[(:p1,1),(:ca,1)], [(:p1,2),(:cb,1)],
             [(:ca,2),(:cb,2), Ground]])
        sol = dcsolve(c, [(mode=(0,), port=1, current=Idc)])
        @test sol.solverinfo.converged
        v = sol.dcnodevoltage
        @test isapprox(maximum(v) - minimum(v), Idc*Z0; rtol = 1e-9)
        # the reference is held at zero, so one of them is
        @test count(iszero, v) >= 2   # ground, and the pinned component
    end

    @testset "a conductance with no real DC limit is refused" begin
        # A frequency dependent resistance is evaluated at zero frequency
        # like everywhere else. If what comes back is not a finite real
        # conductance then the component has no direct current behavior to
        # use, and taking a part of it would be inventing one.
        Idc = 1.0e-6
        complexatdc = Circuit(
            [:p1 => Port(1; Z0 = 50.0),
             :r2 => Resistor(FrequencyDependent(w -> 50.0 + 10.0im)),
             :c1 => Capacitor(1.0e-12)],
            [[(:p1,1),(:r2,1),(:c1,1)], [(:p1,2),(:r2,2),(:c1,2), Ground]])
        @test_throws ArgumentError dcsolve(complexatdc,
            [(mode=(0,), port=1, current=Idc)])

        # and one whose limit is real is fine, and carries its share
        realatdc = Circuit(
            [:p1 => Port(1; Z0 = 50.0),
             :r2 => Resistor(FrequencyDependent(w -> 100.0*(1 + (w/1e11)^2))),
             :c1 => Capacitor(1.0e-12)],
            [[(:p1,1),(:r2,1),(:c1,1)], [(:p1,2),(:r2,2),(:c1,2), Ground]])
        sol = dcsolve(realatdc, [(mode=(0,), port=1, current=Idc)])
        @test sol.solverinfo.converged
        @test isapprox(sol.dcnodevoltage[2], Idc/(1/50.0 + 1/100.0);
            rtol = 1e-9)
    end

    # The rows themselves, apart from a solve. They are the component sum of
    # the zero frequency Kirchhoff equations, so solving them gives the
    # average voltages directly, and the coupling gives the current those
    # voltages drive into the nodes.
    @testset "the transport rows give the average voltages" begin
        JC = JosephsonCircuits
        parts(c, sources) = JC.hbnlsolve(ws, (1,), sources, c,
            Dict{Any,Any}(); keyedarrays = false, dc = true, odd = true,
            returnsystem = true)

        function agrees(c, sources)
            d = parts(c, sources)
            plan = d.dcplan
            isnothing(plan) && return missing
            t = JC.transportrows(plan, d.bnmsource, d.Nmodes)

            # the pinned rows make it nonsingular even where `Y` is not
            @test isfinite(cond(t.Ytr))
            v = t.Ytr \ t.jtr

            # the residual of the rows vanishes at that solution
            Fv = similar(v)
            JC.transportresidual!(Fv, t, v)
            @test all(iszero, Fv)

            # and the coupling maps the solved voltages back to the
            # currents the unpinned transport rows asked for: `Ytr` is
            # `P' G0 P` away from a pinned row, and `transportcurrent!` is
            # `G0 P`, so projecting it recovers the injected current there
            d2 = zeros(size(t.coupling, 1))
            JC.transportcurrent!(d2, t, v)
            proj = transpose(Matrix(plan.lift)) * d2
            for k in eachindex(t.jtr)
                (k in t.pinned) && continue
                @test isapprox(proj[k], t.jtr[k];
                    atol = 1e-8*max(1, maximum(abs, t.jtr)))
            end
            return t
        end

        # a grounded island of two components joined by a bridge resistor
        Rb, Rbig, Idc = 100.0, 1.0e9, 1.0e-6
        cg = Circuit(
            [:p1 => Port(1; Z0 = Rbig), :p2 => Port(2; Z0 = Rbig),
             :rb => Resistor(Rb), :c1 => Capacitor(1e-12),
             :c2 => Capacitor(1e-12)],
            [[(:p1,1),(:rb,1),(:c1,1)], [(:rb,2),(:p2,1),(:c2,1)],
             [(:p1,2),(:p2,2),(:c1,2),(:c2,2), Ground]])
        t = agrees(cg, [(mode=(0,), port=1, current=Idc),
                        (mode=(0,), port=2, current=-Idc)])
        @test t.pinned == Int[]              # it reaches ground, so none

        # a floating island, where only differences are physical and the
        # elimination holds the lowest numbered component at zero
        Z0 = 200.0
        cf = Circuit(
            [:p1 => Port(1; Z0 = Z0), :ca => Capacitor(1e-12),
             :cb => Capacitor(1e-12)],
            [[(:p1,1),(:ca,1)], [(:p1,2),(:cb,1)],
             [(:ca,2),(:cb,2), Ground]])
        t = agrees(cf, [(mode=(0,), port=1, current=Idc)])
        @test length(t.pinned) == 1
        @test (t.Ytr \ t.jtr)[only(t.pinned)] == 0
    end

    # The explicit block end to end: the average voltages carried as
    # unknowns, their transport rows solved inside Newton, and their
    # resistor current reaching the nodes through the coupling. The answers
    # are analytic, so they are asserted directly rather than against a
    # second implementation.
    @testset "the explicit block solves the direct current" begin
        solve(c, sources) = hbnlsolve(ws, (1,), sources, c, Dict{Any,Any}();
            keyedarrays = false, dc = true, odd = true)

        R, Idc = 50.0, 1.0e-6
        # a port environment carrying the whole injected current
        a = solve(Circuit([:p1 => Port(1; Z0 = R), :c1 => Capacitor(1.0e-12)],
                [[(:p1,1),(:c1,1)], [(:p1,2),(:c1,2), Ground]]),
            [(mode=(0,), port=1, current=Idc)])
        @test a.solverinfo.converged
        @test isapprox(maximum(abs, a.dcnodevoltage), Idc*R; rtol = 1e-9)

        # a bridge between two floating components develops I*R across
        # itself, the port environments diverting a part in ten million
        Rb, Rbig = 100.0, 1.0e9
        b = solve(Circuit(
                [:p1 => Port(1; Z0 = Rbig), :p2 => Port(2; Z0 = Rbig),
                 :rb => Resistor(Rb), :c1 => Capacitor(1e-12),
                 :c2 => Capacitor(1e-12)],
                [[(:p1,1),(:rb,1),(:c1,1)], [(:rb,2),(:p2,1),(:c2,1)],
                 [(:p1,2),(:p2,2),(:c1,2),(:c2,2), Ground]]),
            [(mode=(0,), port=1, current=Idc),
             (mode=(0,), port=2, current=-Idc)])
        @test b.solverinfo.converged
        v = b.dcnodevoltage
        @test isapprox(maximum(v) - minimum(v), Idc*Rb; rtol = 1e-5)

        # a floating island fixes only its differences, the pinned component
        # deciding the reference
        c = solve(Circuit(
                [:p1 => Port(1; Z0 = 200.0), :ca => Capacitor(1e-12),
                 :cb => Capacitor(1e-12)],
                [[(:p1,1),(:ca,1)], [(:p1,2),(:cb,1)],
                 [(:ca,2),(:cb,2), Ground]]),
            [(mode=(0,), port=1, current=Idc)])
        @test c.solverinfo.converged
        vc = c.dcnodevoltage
        @test isapprox(maximum(vc) - minimum(vc), Idc*200.0; rtol = 1e-9)
        @test count(iszero, vc) >= 2       # ground, and the pinned component

        # a drive with no direct current in it does not build the block at
        # all, and reports no voltages
        d = solve(Circuit([:p1 => Port(1; Z0 = R), :c1 => Capacitor(1.0e-12)],
                [[(:p1,1),(:c1,1)], [(:p1,2),(:c1,2), Ground]]),
            [(mode=(1,), port=1, current=Idc)])
        @test d.solverinfo.converged
        @test isnothing(d.dcnodevoltage)
    end

    # The capability the explicit block exists for. A scattering block used
    # to be an open circuit at direct current, because its zero frequency
    # row was `i = 0`; with an average voltage to respond to it obeys its
    # own relation instead, and a block which is a resistor carries what
    # that resistor would.
    @testset "a scattering block carries direct current" begin
        JC = JosephsonCircuits
        Rb, Idc, Zbig = 100.0, 1.0e-6, 1.0e9

        # the same circuit twice: once with the resistor, once with a
        # scattering block that is the resistor
        asres = Circuit(
            [:p1 => Port(1; Z0 = Zbig), :rb => Resistor(Rb),
             :c1 => Capacitor(1e-12)],
            [[(:p1,1),(:rb,1),(:c1,1)], [(:rb,2), Ground],
             [(:p1,2),(:c1,2), Ground]])
        blk = ScatteringParameters(
            w -> JC.ABCDtoS(JC.ABCD_seriesZ(Rb + 0im));
            nports = 2, grounded = true, noise = Lossless())
        asblk = Circuit(
            [:p1 => Port(1; Z0 = Zbig), :b => blk, :c1 => Capacitor(1e-12)],
            [[(:p1,1),(:b,1),(:c1,1)], [(:b,2), Ground],
             [(:p1,2),(:c1,2), Ground]])
        srcs = [(mode=(0,), port=1, current=Idc)]
        # The explicit path carries the applied direct current source, so
        # its initial residual is the size of that source in the solver's
        # units, about 1e8 here, and no step drives a residual below about
        # `eps` times its own size. It reaches 2.4e-16 of the initial
        # residual at every drive level, which is machine precision, so the
        # relative test is the one that means anything here.
        go(c) = hbnlsolve(ws, (1,), srcs, c, Dict{Any,Any}();
            keyedarrays = false, dc = true, odd = true, rtol = 1e-12)

        # the block carries the current and develops I*R, which is what it
        # would do if it were the resistor it describes
        got = go(asblk)
        @test got.solverinfo.converged
        @test isapprox(maximum(abs, got.dcnodevoltage), Idc*Rb; rtol = 1e-4)

        # and it agrees with that resistor
        r = go(asres)
        @test r.solverinfo.converged
        @test isapprox(maximum(abs, got.dcnodevoltage),
            maximum(abs, r.dcnodevoltage); rtol = 1e-6)

        # the number to beat: were the block an open, as it was before the
        # explicit rows, the whole current would go through the port
        # environment instead and the voltage would be Idc*Zbig, seven
        # orders larger
        @test Idc*Zbig / (Idc*Rb) > 1e6

        # the relative test is drive independent, which the absolute one is
        # not: the residual it stops at moves with the source while the
        # accuracy does not
        for I in (1e-9, 1e-3)
            big = hbnlsolve(ws, (1,), [(mode=(0,), port=1, current=I)],
                asblk, Dict{Any,Any}(); keyedarrays = false, dc = true,
                odd = true, rtol = 1e-12)
            @test big.solverinfo.converged
            @test isapprox(maximum(abs, big.dcnodevoltage), I*Rb;
                rtol = 1e-4)
        end
    end

    # A block whose own row does not determine its current: `C(0)` is
    # singular, so the block constrains its voltage and leaves a current
    # direction free. Whether that is solvable is not a property of the
    # block -- it is whether anything else sees the free direction, which is
    # the rank of the direct current subsystem and not a taxonomy of block
    # types.
    @testset "a free port current is determined, or named" begin
        JC = JosephsonCircuits
        R, Idc = 100.0, 1.0e-6
        short() = ScatteringParameters(
            JC.S_short!(ones(Complex{Float64}, 1, 1));
            nports = 1, grounded = true, noise = Lossless())
        go(c) = hbnlsolve(ws, (1,), [(mode=(0,), port=1, current=Idc)],
            c, Dict{Any,Any}(); keyedarrays = false, dc = true, odd = true,
            rtol = 1e-12)

        # A short to ground carries the injected current away. Its own row
        # is `B0 V = 0`, which pins the voltage and says nothing about the
        # current; the transport row of the component determines it, because
        # the current crosses the component's boundary on its way to ground.
        withshort = Circuit(
            [:p1 => Port(1; Z0 = 1.0e9), :r => Resistor(R), :s => short(),
             :c1 => Capacitor(1e-12), :c2 => Capacitor(1e-12)],
            [[(:p1,1),(:r,1),(:c1,1)], [(:r,2),(:s,1),(:c2,1)],
             [(:p1,2),(:c1,2),(:c2,2), Ground]])

        # a path to ground through a scattering block is not in the
        # conductance graph, which is why this circuit was refused before
        # the block's own row existed. It solves now: the driven node sits
        # at I*R and the shorted one at zero.
        s = go(withshort)
        @test s.solverinfo.converged
        v = s.dcnodevoltage
        @test isapprox(maximum(v), Idc*R; rtol = 1e-4)
        @test count(iszero, v) >= 2          # ground, and the shorted node

        # An ideal through in parallel with an inductor is two ideal shorts
        # between the same pair of nodes. Every node voltage is determined
        # and the division of current between the two is not: the free
        # direction stays inside one static flux component and cancels in
        # its transport row. That is a gauge, not an error, so it is pinned
        # -- one of the equivalent divisions is chosen -- and the voltages,
        # which are what is observable, come out the same either way.
        through() = ScatteringParameters(
            w -> JC.ABCDtoS(JC.ABCD_seriesZ(0.0 + 0im));
            nports = 2, grounded = true, noise = Lossless())
        undetermined = Circuit(
            [:p1 => Port(1; Z0 = R), :l => Inductor(1e-9), :t => through(),
             :c1 => Capacitor(1e-12), :c2 => Capacitor(1e-12)],
            [[(:p1,1),(:l,1),(:t,1),(:c1,1)], [(:l,2),(:t,2),(:c2,1)],
             [(:p1,2),(:c1,2),(:c2,2), Ground]])
        b = go(undetermined)
        @test b.solverinfo.converged
        # the two nodes are tied by two ideal shorts, so their average
        # voltages are equal and sit at the drive across the environment
        v = b.dcnodevoltage
        @test isapprox(maximum(v), Idc*R; rtol = 1e-6)
        @test isapprox(v[2], v[3]; rtol = 1e-9)

        # the pinning is what makes that possible: without it the direct
        # current subsystem is singular
        @test JC.freecurrents(JC.dcblockrows(
            JC.scatteringstampsystem(
                JC.compile(undetermined).scatteringblocks, 1;
                auxoffset = 0, Ntotal = 1000, scale = 1.0).blocks,
            zeros(Int, 8), 1, 1, 0, 0, 1.0)) > 0
    end

    # Singular and solvable, against singular and not. The first is pinned
    # and the second refused, and the difference is whether the constant
    # side lies in the range of the matrix.
    #
    # This is exercised on a doctored subsystem rather than a circuit,
    # deliberately. Every direct current source is a port drive and every
    # port owns its termination, which is a resistive path between that
    # port's own terminals, so the injection into a resistive island is
    # balanced by construction and the inconsistent case cannot be reached
    # through the public interface today. It is a guard against an input
    # which does not exist yet, and it is tested as one.
    @testset "an unsolvable direct current block is refused" begin
        JC = JosephsonCircuits

        # a real plan, so the structure around it is genuine
        c = Circuit(
            [:p1 => Port(1; Z0 = 200.0), :ca => Capacitor(1e-12),
             :cb => Capacitor(1e-12)],
            [[(:p1,1),(:ca,1)], [(:p1,2),(:cb,1)],
             [(:ca,2),(:cb,2), Ground]])
        d = hbnlsolve(ws, (1,), [(mode=(0,), port=1, current=1e-6)], c,
            Dict{Any,Any}(); keyedarrays = false, dc = true, odd = true,
            returnsystem = true)
        plan = d.dcplan
        @test !isnothing(plan)
        real = JC.transportrows(plan, d.bnmsource, d.Nmodes)
        nc = length(real.jtr)
        @test nc == 2

        # a two component block whose rows say only that the difference is
        # fixed: the sum is a direction no equation sees
        singular = [1.0 -1.0; -1.0 1.0]
        L = JC.compositelayout(JC.ModeLayout([true], 4), [(0,)]; nvdc = nc)
        function work(Y, j)
            t = JC.TransportRows(plan, Y, j, real.coupling, Int[])
            return JC.CanonicalWork(L, zeros(L.rdim); transport = t,
                nnodaldc = L.ndc)
        end

        # a constant side along the difference is in the range: solvable,
        # and the free direction is pinned rather than refused
        w = work(copy(singular), [1.0, -1.0])
        @test !isnothing(w.pinning)
        @test size(w.pinning.N, 2) == 1
        # the pinned direction is the one nothing determines, the sum
        nvec = vec(w.pinning.N)
        @test isapprox(abs(nvec[1]), abs(nvec[2]); rtol = 1e-8)
        @test sign(nvec[1]) == sign(nvec[2])

        # a constant side along the sum is not in the range: no solution
        @test_throws ArgumentError work(copy(singular), [1.0, 1.0])
    end
end
