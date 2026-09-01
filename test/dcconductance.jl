using JosephsonCircuits
using LinearAlgebra
using Random
using Test

# Direct current through resistors. The harmonic balance state is periodic
# node flux, so a voltage is its time derivative and vanishes at zero
# frequency: a resistor was an open circuit at DC and a current source
# driving one could not develop I*R. The missing coordinate is the average
# voltage, which is constant on each static flux component because an
# inductor or a zero-voltage junction is a short there.

# an inner preconditioner which does nothing, so what is asserted below is
# the direct current step and not the mode coupling solve around it
struct Passthrough <: JosephsonCircuits.AbstractPreconditioner end
JosephsonCircuits.applypreconditioner!(z, ::Passthrough, r) = copyto!(z, r)
JosephsonCircuits.updatepreconditioner!(pc::Passthrough, x) = pc

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
            # ground is excluded, as it is from the node flux, so the
            # first entry is the first real node
            @test length(sol.dcnodevoltage) == 1
            @test isapprox(only(sol.dcnodevoltage), sgn*Idc*R; rtol = 1e-9)
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
        @test isapprox(v[2], Vtop; rtol = 1e-9)
        @test isapprox(v[1], Vtop*Rl/(Rs+Rl); rtol = 1e-9)
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

    @testset "the tolerance follows the scale of the source" begin
        # The scale which nondimensionalizes the system is read off the port
        # reference impedances, so a circuit whose interior sits far from
        # them -- a hundred ohm bridge between ports made nearly open at a
        # teraohm -- is left with a scaled source many orders above one. The
        # residual cannot be pushed below the rounding error of a sum of
        # terms that size, and an absolute tolerance under that floor makes
        # convergence a matter of which way the last rounding fell rather
        # than of whether the iteration found the answer. What is asserted
        # is that it is reported converged and that the answer is the exact
        # one, on both of the methods which carry the block.
        Rb, Rbig, Idc = 100.0, 1.0e12, 1.0e-6
        c = Circuit(
            [:p1 => Port(1; Z0 = Rbig), :p2 => Port(2; Z0 = Rbig),
             :rb => Resistor(Rb), :c1 => Capacitor(1e-12),
             :c2 => Capacitor(1e-12)],
            [[(:p1,1),(:rb,1),(:c1,1)], [(:rb,2),(:p2,1),(:c2,1)],
             [(:p1,2),(:p2,2),(:c1,2),(:c2,2), Ground]])
        for m in (:newton, :newtonkrylov)
            sol = hbnlsolve(ws, (1,), [(mode=(0,), port=1, current=Idc),
                    (mode=(0,), port=2, current=-Idc)], c, Dict{Any,Any}();
                keyedarrays = false, dc = true, odd = true, method = m)
            @test sol.solverinfo.converged
            # and it stops where the arithmetic runs out, not earlier
            @test sol.solverinfo.finalresidual <=
                1e-14*sol.solverinfo.initialresidual
            v = sol.dcnodevoltage
            @test isapprox(maximum(v) - minimum(v), Idc*Rb; rtol = 1e-8)
        end
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
        @test count(iszero, v) >= 1   # the component held as the reference
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
        @test isapprox(sol.dcnodevoltage[1], Idc/(1/50.0 + 1/100.0);
            rtol = 1e-9)
    end

    # The rows themselves, apart from a solve. They are the component sum of
    # the zero frequency Kirchhoff equations, so they give the average
    # voltages directly, and the coupling gives the current those voltages
    # drive into the nodes. They carry no reference: which of them is
    # redundant depends on what else is in the descriptor, and that is
    # decided later, by `dcpinning`.
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

            # a least squares solution, which exists whether or not the rows
            # fix an absolute potential
            v = pinv(t.Y) * t.j

            # the residual of the rows vanishes at that solution
            Fv = similar(v)
            JC.transportresidual!(Fv, t, v)
            @test all(x -> abs(x) <= 1e-8*max(1, maximum(abs, t.j)), Fv)

            # and the coupling maps the solved voltages back to the currents
            # the rows asked for: `Y` is `P' G0 P` and `transportcurrent!`
            # is `G0 P`, so projecting it recovers the injected current
            d2 = zeros(size(t.coupling, 1))
            JC.transportcurrent!(d2, t, v)
            proj = transpose(Matrix(plan.lift)) * d2
            for k in eachindex(t.j)
                @test isapprox(proj[k], t.j[k];
                    atol = 1e-8*max(1, maximum(abs, t.j)))
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
        # it reaches ground, so the rows already fix an absolute potential
        @test isfinite(cond(t.Y))

        # a floating island, where only differences are physical: the rows
        # say exactly that, and are singular by one direction
        Z0 = 200.0
        cf = Circuit(
            [:p1 => Port(1; Z0 = Z0), :ca => Capacitor(1e-12),
             :cb => Capacitor(1e-12)],
            [[(:p1,1),(:ca,1)], [(:p1,2),(:cb,1)],
             [(:ca,2),(:cb,2), Ground]])
        t = agrees(cf, [(mode=(0,), port=1, current=Idc)])
        @test rank(t.Y) == size(t.Y, 1) - 1
        @test isapprox(t.Y * ones(size(t.Y, 2)), zeros(size(t.Y, 1));
            atol = 1e-12*maximum(abs, t.Y))
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
        @test count(iszero, vc) >= 1       # the reference component

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
        @test count(iszero, v) >= 1          # the shorted node

        # An ideal through in parallel with an inductor leaves the current
        # divided between them undetermined, and that division is not a
        # gauge. It cancels in the transport row, because both terminals lie
        # in one static flux component, but it does not cancel at the nodes:
        # moving it injects +d at one and -d at the other, which moves the
        # inductor current and the static flux across it. So the subsystem
        # is singular in a direction the rest of the circuit can see, and
        # the circuit is refused rather than given one of infinitely many
        # answers. What is pinned instead is a direction the nodes cannot
        # see, which is a floating island's common voltage; the difference
        # is `H N`, not the block type.
        through() = ScatteringParameters(
            w -> JC.ABCDtoS(JC.ABCD_seriesZ(0.0 + 0im));
            nports = 2, grounded = true, noise = Lossless())
        undetermined = Circuit(
            [:p1 => Port(1; Z0 = R), :l => Inductor(1e-9), :t => through(),
             :c1 => Capacitor(1e-12), :c2 => Capacitor(1e-12)],
            [[(:p1,1),(:l,1),(:t,1),(:c1,1)], [(:l,2),(:t,2),(:c2,1)],
             [(:p1,2),(:c1,2),(:c2,2), Ground]])
        @test_throws ArgumentError go(undetermined)

        # giving the block a finite series impedance determines the division
        # and the same circuit solves
        finite() = ScatteringParameters(
            w -> JC.ABCDtoS(JC.ABCD_seriesZ(10.0 + 0im));
            nports = 2, grounded = true, noise = Lossless())
        determined = Circuit(
            [:p1 => Port(1; Z0 = R), :l => Inductor(1e-9), :t => finite(),
             :c1 => Capacitor(1e-12), :c2 => Capacitor(1e-12)],
            [[(:p1,1),(:l,1),(:t,1),(:c1,1)], [(:l,2),(:t,2),(:c2,1)],
             [(:p1,2),(:c1,2),(:c2,2), Ground]])
        b = go(determined)
        @test b.solverinfo.converged
        # the inductor is a short at zero frequency, so the two nodes still
        # sit together at the drive across the environment
        v = b.dcnodevoltage
        @test isapprox(maximum(v), Idc*R; rtol = 1e-6)
        @test isapprox(v[1], v[2]; rtol = 1e-9)

    end

    # The reference is chosen after the whole descriptor exists, and not
    # from the resistors alone. These two circuits are the cases where the
    # difference shows: in the first a block makes a row the resistors
    # called redundant necessary, and in the second there are no resistors
    # to ask.
    @testset "the descriptor is complete before a reference is chosen" begin
        JC = JosephsonCircuits
        R, Idc = 100.0, 1.0e-6
        go(c) = hbnlsolve(ws, (1,), [(mode=(0,), port=1, current=Idc)], c,
            Dict{Any,Any}(); keyedarrays = false, dc = true, odd = true,
            rtol = 1e-12)

        # An ideal through from a driven node to one which has no direct
        # current path of its own. In the conductance graph alone the second
        # node is a floating island and its transport row looks redundant,
        # so choosing a reference there would replace that row by `v = 0`
        # and drop the through current which lands in it. The physical
        # answer is that no current crosses -- the second node has nowhere
        # to send it -- and both nodes sit at `I*R`.
        through = ScatteringParameters(
            w -> JC.ABCDtoS(JC.ABCD_seriesZ(0.0 + 0im));
            nports = 2, grounded = true, noise = Lossless())
        bridged = Circuit(
            [:p1 => Port(1; Z0 = R), :t => through,
             :c1 => Capacitor(1e-12), :c2 => Capacitor(1e-12)],
            [[(:p1,1),(:t,1),(:c1,1)], [(:t,2),(:c2,1)],
             [(:p1,2),(:c1,2),(:c2,2), Ground]])
        a = go(bridged)
        @test a.solverinfo.converged
        v = a.dcnodevoltage
        @test isapprox(maximum(v), Idc*R; rtol = 1e-6)
        @test isapprox(v[1], v[2]; rtol = 1e-9)   # the through ties them

        # A circuit with no conductance at all: the port declares its
        # reference impedance and loads nothing, so `G0` is empty and the
        # only direct current path is the block. Gating the plan on `G0`
        # left the artificial `i = 0` rows of the stamp in place and the
        # injected current with nowhere to go.
        shorted = Circuit(
            [:p1 => Port(1; Z0 = R, termination = nothing),
             :s => ScatteringParameters(JC.S_short!(ones(Complex{Float64},1,1));
                 nports = 1, grounded = true, noise = Lossless()),
             :c1 => Capacitor(1e-12)],
            [[(:p1,1),(:s,1),(:c1,1)], [(:p1,2),(:c1,2), Ground]])
        b = go(shorted)
        @test b.solverinfo.converged
        # the short holds the driven node at ground, and that is an answer:
        # a solution which is zero is not the same as no solution
        @test !isnothing(b.dcnodevoltage)
        @test all(iszero, b.dcnodevoltage)
    end

    # The rank decisions above are structural facts about the circuit, so
    # they must not move with the units it is written in. The subsystem
    # mixes volts and amperes and its entries carry whatever impedance scale
    # the circuit happens to use, which is why the rows and columns are
    # equilibrated before anything is decided.
    @testset "the classification does not depend on the impedance scale" begin
        JC = JosephsonCircuits
        through = ScatteringParameters(
            w -> JC.ABCDtoS(JC.ABCD_seriesZ(0.0 + 0im));
            nports = 2, grounded = true, noise = Lossless())
        for k in (1e-3, 1.0, 1e3)
            R, L, C, I = 100.0*k, 1e-9*k, 1e-12/k, 1e-6/k
            go(c) = hbnlsolve(ws, (1,), [(mode=(0,), port=1, current=I)], c,
                Dict{Any,Any}(); keyedarrays = false, dc = true, odd = true,
                rtol = 1e-12)

            # undetermined at every scale
            @test_throws ArgumentError go(Circuit(
                [:p1 => Port(1; Z0 = R), :l => Inductor(L), :t => through,
                 :c1 => Capacitor(C), :c2 => Capacitor(C)],
                [[(:p1,1),(:l,1),(:t,1),(:c1,1)], [(:l,2),(:t,2),(:c2,1)],
                 [(:p1,2),(:c1,2),(:c2,2), Ground]]))

            # and a gauge at every scale, giving the same answer in the
            # units of the problem
            sol = go(Circuit(
                [:p1 => Port(1; Z0 = R), :ca => Capacitor(C),
                 :cb => Capacitor(C)],
                [[(:p1,1),(:ca,1)], [(:p1,2),(:cb,1)],
                 [(:ca,2),(:cb,2), Ground]]))
            @test sol.solverinfo.converged
            v = sol.dcnodevoltage
            @test isapprox(maximum(v) - minimum(v), I*R; rtol = 1e-9)
        end
    end

    # The waves are in units of sqrt(photons/second), whose normalization
    # has no limit at zero frequency, so there is no direct current wave and
    # the zero frequency entries of the pump scattering matrix are
    # identically zero. That is a convention rather than a computation, and
    # it is what keeps the direct current out: the voltage the wave
    # extractor reconstructs is `im*w*phi`, which is zero at zero frequency
    # however large the average voltage is, so a nonzero entry there would
    # be an artifact rather than a measurement.
    @testset "there is no direct current wave" begin
        c = Circuit(
            [:p1 => Port(1; Z0 = 50.0), :cc => Capacitor(100e-15),
             :jj => JosephsonJunction(1000e-12), :cj => Capacitor(1000e-15)],
            [[(:p1,1),(:cc,1)], [(:cc,2),(:jj,1),(:cj,1)],
             [(:p1,2),(:jj,2),(:cj,2), Ground]])
        sol = hbnlsolve((2*pi*4.75e9,), (4,),
            [(mode=(1,), port=1, current=1.2e-6),
             (mode=(0,), port=1, current=1.0e-7)],
            c, Dict{Any,Any}(); dc = true, odd = true, even = true,
            keyedarrays = false, rtol = 1e-12)
        @test sol.solverinfo.converged
        dc = findfirst(m -> all(iszero, m), sol.frequencies.modes)
        @test !isnothing(dc)
        # exactly zero, not a small number and not a NaN from dividing one
        # vanishing wave by another
        @test all(iszero, sol.S[dc, :])
        @test all(iszero, sol.S[:, dc])
        @test !any(isnan, sol.S)
        # and the operating point it stands for is not zero, so the zeros
        # above are the convention and not the absence of direct current
        @test isapprox(maximum(sol.dcnodevoltage), 1.0e-7*50.0; rtol = 1e-6)
    end

    # The static flux partition treats a junction as a short at zero
    # frequency, which assumes it is in the zero voltage state. A junction
    # asked for more direct current than its critical current has no such
    # state, and the solver cannot say so on its own: the branch current is
    # Ic*sin(phi), which is bounded, so it converges to the nearest periodic
    # thing rather than reporting that none exists. What can be reported is
    # the approach to that edge.
    @testset "a junction near its critical current is reported" begin
        Lj = 1000e-12
        Ic = real(JosephsonCircuits.LjtoIc(Lj))
        c = Circuit(
            [:p1 => Port(1; Z0 = 50.0), :jj => JosephsonJunction(Lj),
             :cj => Capacitor(1000e-15)],
            [[(:p1,1),(:jj,1),(:cj,1)], [(:p1,2),(:jj,2),(:cj,2), Ground]])
        go(I) = hbnlsolve(ws, (4,), [(mode=(0,), port=1, current=I)], c,
            Dict{Any,Any}(); keyedarrays = false, dc = true, odd = true,
            even = true)

        # well inside the zero voltage state: nothing to say
        a = @test_logs go(0.3*Ic)
        @test a.solverinfo.converged

        # at the edge of it: named, with the fraction it reached
        b = (@test_logs (:warn,) go(0.996*Ic))
        @test b.solverinfo.converged

        # and past it there is no periodic solution to find
        d = (@test_logs (:warn,) match_mode = :any go(1.03*Ic))
        @test !d.solverinfo.converged
    end

    # A block's direct current behavior is its scattering matrix at zero,
    # and the default is to ask the block for it. Two kinds of block cannot
    # answer: measured data which starts at gigahertz has no entry there,
    # and a closed form whose limit exists may not be evaluable there. So
    # the limit can be stated instead.
    @testset "a block may state its zero frequency model" begin
        JC = JosephsonCircuits
        R, Idc = 100.0, 1.0e-6
        # a series capacitance in closed form: an open circuit at direct
        # current, and infinite at zero
        blk(dc) = ScatteringParameters(
            w -> JC.ABCDtoS(JC.ABCD_seriesZ(1/(im*w*1e-12)));
            nports = 2, grounded = true, noise = Lossless(), dcmodel = dc)
        mk(b) = Circuit(
            [:p1 => Port(1; Z0 = R), :x => b, :r2 => Resistor(R),
             :c1 => Capacitor(1e-12)],
            [[(:p1,1),(:x,1),(:c1,1)], [(:x,2),(:r2,1)],
             [(:p1,2),(:r2,2),(:c1,2), Ground]])
        go(b) = hbnlsolve(ws, (1,), [(mode=(0,), port=1, current=Idc)], mk(b),
            Dict{Any,Any}(); keyedarrays = false, dc = true, odd = true,
            rtol = 1e-12)

        # asked, it cannot answer, and says what to do about it
        @test_throws ArgumentError go(blk(JC.ScatteringLimit()))

        # an open circuit sends the whole current through the port's own
        # environment and leaves the far node at ground through `r2`
        a = go(blk(OpenDC()))
        @test a.solverinfo.converged
        @test isapprox(a.dcnodevoltage[1], Idc*R; rtol = 1e-9)
        @test isapprox(a.dcnodevoltage[2], 0; atol = 1e-12)

        # a through puts the two resistors in parallel, and both nodes sit
        # at half of that
        b = go(blk(ThroughDC()))
        @test b.solverinfo.converged
        @test isapprox(b.dcnodevoltage[1], Idc*R/2; rtol = 1e-9)
        @test isapprox(b.dcnodevoltage[2], b.dcnodevoltage[1]; rtol = 1e-9)

        # a short holds both of its ports at ground
        c = go(blk(ShortDC()))
        @test c.solverinfo.converged
        @test all(x -> isapprox(x, 0; atol = 1e-12), c.dcnodevoltage)

        # and the named models are shorthand for a stated matrix
        d = go(blk(ScatteringDC([0 1.0; 1.0 0])))
        @test isapprox(d.dcnodevoltage, b.dcnodevoltage; rtol = 1e-12)

        # a stated matrix is data like the block's own, and is checked the
        # same way at construction
        @test_throws ArgumentError ScatteringParameters([0 1;1 0];
            dcmodel = ScatteringDC([0 2.0; 2.0 0]))          # active
        @test_throws DimensionMismatch ScatteringParameters([0 1;1 0];
            dcmodel = ScatteringDC(fill(0.0, 3, 3)))         # wrong size
        @test_throws ArgumentError ScatteringDC([0 im; im 0]) # complex
        @test_throws ArgumentError ScatteringParameters(
            reshape([0.0+0im],1,1); dcmodel = ThroughDC())    # not two port

        # an active zero frequency model is allowed where an active block is
        @test ScatteringParameters([0 1;1 0];
            noise = NoiseCovariance([1.0 0.0; 0.0 1.0]),
            dcmodel = ScatteringDC([0 2.0; 2.0 0])) isa ScatteringParameters
    end

    # The direct current block is held in Float64 whatever precision the
    # periodic solve runs in, because it is small, exactly solved, and the
    # worst conditioned part of the problem. A single precision solve is
    # then as accurate as single precision allows, which is what it should
    # be, and needs a tolerance it can meet: the default is absolute and
    # sized for double precision.
    @testset "a single precision solve carries the block" begin
        c = Circuit(
            [:p1 => Port(1; Z0 = 50.0), :c1 => Capacitor(1.0e-12)],
            [[(:p1,1),(:c1,1)], [(:p1,2),(:c1,2), Ground]])
        src = [(mode=(0,), port=1, current=1e-6),
               (mode=(1,), port=1, current=1e-6)]
        go(P; kw...) = hbnlsolve(ws, (4,), src, c, Dict{Any,Any}();
            keyedarrays = false, dc = true, odd = true, even = true,
            precision = P, kw...)

        a = go(Float64)
        b = go(Float32; rtol = 1e-6)
        @test a.solverinfo.converged
        @test b.solverinfo.converged
        # the same answer, to what single precision can hold
        @test isapprox(only(a.dcnodevoltage), only(b.dcnodevoltage);
            rtol = 1e-6)
        # and it stops where single precision runs out rather than earlier
        r = b.solverinfo.finalresidual/b.solverinfo.initialresidual
        @test r <= eps(Float32)
    end

    # The preconditioner solves the direct current subsystem exactly and
    # writes the answer over whatever the inner preconditioner guessed at
    # those coordinates. On a device the same factors and the same
    # substitutions run there rather than the window crossing the bus, so
    # the two paths are the same arithmetic; here the host path is checked
    # against the factorization it is meant to reproduce.
    @testset "the preconditioner solves the block exactly" begin
        JC = JosephsonCircuits
        Rb, Rbig, Idc = 100.0, 1.0e9, 1.0e-6
        c = Circuit(
            [:p1 => Port(1; Z0 = Rbig), :p2 => Port(2; Z0 = Rbig),
             :rb => Resistor(Rb), :c1 => Capacitor(1e-12),
             :c2 => Capacitor(1e-12)],
            [[(:p1,1),(:rb,1),(:c1,1)], [(:rb,2),(:p2,1),(:c2,1)],
             [(:p1,2),(:p2,2),(:c1,2),(:c2,2), Ground]])
        d = hbnlsolve(ws, (1,), [(mode=(0,), port=1, current=Idc),
                                 (mode=(0,), port=2, current=-Idc)], c,
            Dict{Any,Any}(); keyedarrays = false, dc = true, odd = true,
            returnsystem = true)
        w = d.canonicalwork
        L = w.layout
        pc = JC.CanonicalPreconditioner(Passthrough(), w)

        n = JC.canonicaldim(L)
        Random.seed!(9)
        r = randn(n)
        z = zeros(n)
        JC.applypreconditioner!(z, pc, r)

        # the flux coordinates carry the inner preconditioner's answer
        idx = JC.dcsubsystemindices(w)
        rest = setdiff(1:n, idx)
        @test z[rest] == r[rest]
        # and the subsystem's own are its exact solution
        A = JC.dcsubsystem(w)
        @test z[idx] == A \ r[idx]
        # which is a real solve and not the identity: this subsystem is ill
        # conditioned enough that an inverse would not do
        @test cond(A) > 1e6
        @test z[idx] != r[idx]
    end

    # The operating point of a circuit with a direct current block, and the
    # sensitivities taken there. The implicit function theorem applies to
    # the canonical system, not to the harmonic part of it, so both the
    # residual derivative and the solve have to be in those coordinates: a
    # resistor which carries direct current moves the average voltages, and
    # nothing in the harmonic rows knows that.
    @testset "sensitivities carry the direct current block" begin
        JC = JosephsonCircuits
        R1, R2, Idc, Iac = 50.0, 200.0, 1.0e-6, 1.0e-6
        wp = (2*pi*4.75e9,)
        mk(r2) = Circuit(
            [:p1 => Port(1; Z0 = R1), :r2 => Resistor(r2),
             :c1 => Capacitor(1e-12), :jj => JosephsonJunction(1000e-12),
             :cj => Capacitor(1000e-15)],
            [[(:p1,1),(:r2,1),(:c1,1)], [(:c1,2),(:jj,1),(:cj,1)],
             [(:p1,2),(:r2,2),(:jj,2),(:cj,2), Ground]])
        srcs = [(mode=(0,), port=1, current=Idc),
                (mode=(1,), port=1, current=Iac)]
        go(r2; kw...) = hbnlsolve(wp, (4,), srcs, mk(r2), Dict{Any,Any}();
            dc = true, odd = true, even = true, keyedarrays = false,
            rtol = 1e-13, method = :newton, kw...)

        sol = go(R2; returnoperatingpoint = true)
        op = sol.operatingpoint
        @test !isnothing(op.dc)
        L = op.dc.work.layout
        @test length(op.dc.u) == JC.canonicaldim(L)

        psc = JC.compile(mk(R2))
        cg = calccircuitgraph(psc)
        nm = numericmatrices(psc, cg, Dict{Any,Any}(); Nmodes = op.Nmodes)
        ir2 = findfirst(==("r2"), psc.componentnames)
        dFr = JC.calcresidualsensitivity(op, psc, cg, nm, [ir2])
        # the residual derivative is in the canonical coordinates, and it
        # reaches the transport rows, which the harmonic one cannot
        @test size(dFr, 1) == JC.canonicaldim(L)
        @test !iszero(dFr[end, 1])

        dx = JC.calcnodefluxsensitivity(op, dFr)
        dv = JC.dcvoltagesensitivity(op, dFr)

        # against a central difference of a re-solve, in the same relative
        # parameter the sensitivity is taken in
        h = 1e-6
        sp = go(R2*(1+h)); sm = go(R2*(1-h))
        fdv = (sp.dcnodevoltage .- sm.dcnodevoltage)./(2h)
        fdflux = (sp.nodeflux .- sm.nodeflux)./(2h)
        @test isapprox(vec(dv), fdv; rtol = 1e-5, atol = 1e-12)
        @test isapprox(dx, fdflux; rtol = 1e-4,
            atol = 1e-8*maximum(abs, fdflux))

        # the analytic value: two resistors in parallel carry the direct
        # current, so d(R1||R2)/dlog(R2) is R2*R1^2/(R1+R2)^2
        @test isapprox(maximum(abs, dv), Idc*R2*R1^2/(R1+R2)^2; rtol = 1e-6)

        # a circuit with no block still answers, and has no voltages
        plain = hbnlsolve(wp, (4,),
            [(mode=(1,), port=1, current=Iac)], mk(R2), Dict{Any,Any}();
            dc = true, odd = true, even = true, keyedarrays = false,
            method = :newton, returnoperatingpoint = true)
        @test isnothing(plain.operatingpoint.dc)
        @test isnothing(JC.dcvoltagesensitivity(plain.operatingpoint,
            zeros(1,1)))
    end

    # The whole analysis, not just the nonlinear solve: `hbsolve` runs the
    # linearized solve against the operating point, and its scattering
    # parameter sensitivities contract through the Jacobian of the system
    # which was solved -- the canonical one when a block is active.
    @testset "the linearized analysis carries the direct current block" begin
        JC = JosephsonCircuits
        R1, R2, Idc, Iac = 50.0, 200.0, 1.0e-6, 1.2e-6
        wp = (2*pi*4.75e9,)
        wsig = 2*pi*[4.6e9, 5.0e9]
        mk(r2) = Circuit(
            [:p1 => Port(1; Z0 = R1), :r2 => Resistor(r2),
             :c1 => Capacitor(100e-15), :jj => JosephsonJunction(1000e-12),
             :cj => Capacitor(1000e-15)],
            [[(:p1,1),(:r2,1),(:c1,1)], [(:c1,2),(:jj,1),(:cj,1)],
             [(:p1,2),(:r2,2),(:jj,2),(:cj,2), Ground]])
        srcs = [(mode=(0,), port=1, current=Idc),
                (mode=(1,), port=1, current=Iac)]
        go(r2; kw...) = hbsolve(wsig, wp, srcs, (2,), (4,), mk(r2),
            Dict{Any,Any}(); dc = true, keyedarrays = false, kw...)

        # the drive frequency survives the solve. It did not: the direct
        # current work was assigned to a local named `w`, which is this
        # function's drive frequency, so every circuit which injected direct
        # current reported the work object as its drive and the linearized
        # solve could not compute its mode frequencies from it. No circuit
        # with a direct current drive reached `hbsolve` at all.
        nl = hbnlsolve(wp, (4,), srcs, mk(R2), Dict{Any,Any}();
            dc = true, odd = true, even = true, keyedarrays = false)
        @test nl.w isa Tuple
        @test only(nl.w) ≈ only(wp)

        base = go(R2)
        @test all(isfinite, Array(base.linearized.S))

        # the two contraction orders are the same computation
        fw = go(R2; sensitivitynames = ["r2"], returnSsensitivity = true,
            sensitivityoperatingpoint = true, sensitivitymode = :forward)
        rv = go(R2; sensitivitynames = ["r2"], returnSsensitivity = true,
            sensitivityoperatingpoint = true, sensitivitymode = :reverse)
        dSf = Array(fw.linearized.Ssensitivity)
        dSr = Array(rv.linearized.Ssensitivity)
        @test maximum(abs, dSf .- dSr) < 1e-12*maximum(abs, dSr)

        # and both are the derivative of the whole analysis, against a
        # central difference of a re-solve in the same relative parameter
        h = 1e-5
        Sp = Array(go(R2*(1+h)).linearized.S)
        Sm = Array(go(R2*(1-h)).linearized.S)
        fd = (Sp .- Sm)./(2h)
        @test isapprox(dSr[:,:,1,:], fd; rtol = 1e-6,
            atol = 1e-8*maximum(abs, fd))
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
        nc = length(real.j)
        @test nc == 2

        # a two component block whose rows say only that the difference is
        # fixed: the sum is a direction no equation sees
        singular = [1.0 -1.0; -1.0 1.0]
        L = JC.compositelayout(JC.ModeLayout([true], 4), [(0,)]; nvdc = nc)
        function work(Y, j)
            t = JC.TransportRows(plan, Y, j, real.coupling)
            return JC.CanonicalWork(L, zeros(L.rdim); transport = t,
                nnodaldc = L.ndc)
        end

        # a constant side along the difference is in the range: solvable,
        # and the free direction is pinned rather than refused
        w = work(copy(singular), [1.0, -1.0])
        @test !isnothing(w.pinning)
        @test length(w.pinning.rows) == 1
        # one redundant equation is given up, and one voltage coordinate is
        # held at zero in its place. The undetermined direction is the sum,
        # which the coupling of this floating island cannot see, so it is a
        # gauge and either coordinate is a legitimate reference
        @test only(w.pinning.cols) in 1:nc
        @test only(w.pinning.rows) in 1:nc

        # a constant side along the sum is not in the range: no solution
        @test_throws ArgumentError work(copy(singular), [1.0, 1.0])
    end

    # The average voltages and the blocks' zero frequency rows live outside
    # the `HBSystem`, so an interface which hands that object out, or stores
    # it to be differentiated later, is describing a different problem from
    # the one solved. `returnsystem` and `debugJacobian` hand back the
    # canonical work beside it; `returnoperatingpoint` and sensitivities
    # cannot yet, and say so.
    @testset "interfaces do not hand out the unaugmented system" begin
        c = Circuit(
            [:p1 => Port(1; Z0 = 50.0), :c1 => Capacitor(1.0e-12)],
            [[(:p1,1),(:c1,1)], [(:p1,2),(:c1,2), Ground]])
        dcsrc = [(mode=(0,), port=1, current=1.0e-6)]
        acsrc = [(mode=(1,), port=1, current=1.0e-6)]

        # with direct current, the parts interface carries the direct
        # current block too
        d = hbnlsolve(ws, (1,), dcsrc, c, Dict{Any,Any}();
            keyedarrays = false, dc = true, odd = true, returnsystem = true)
        @test d.dcexplicit
        @test !isnothing(d.canonicalwork)
        @test d.canonicalwork.layout.nvdc > 0

        # without it there is no block and nothing to carry
        a = hbnlsolve(ws, (1,), acsrc, c, Dict{Any,Any}();
            keyedarrays = false, dc = true, odd = true, returnsystem = true)
        @test !a.dcexplicit
        @test isnothing(a.canonicalwork)

        # the operating point carries the block rather than defining a
        # second one without it
        op = hbnlsolve(ws, (1,), dcsrc, c, Dict{Any,Any}();
            keyedarrays = false, dc = true, odd = true,
            returnoperatingpoint = true).operatingpoint
        @test !isnothing(op.dc)
        @test size(op.dc.jacobian, 1) ==
            JosephsonCircuits.canonicaldim(op.dc.work.layout)

        # and sensitivities are taken through it rather than refused
        sens = hbnlsolve(ws, (1,), dcsrc, c, Dict{Any,Any}();
            keyedarrays = false, dc = true, odd = true,
            sensitivitynames = ["c1"])
        @test sens.solverinfo.converged

        # and both are available when no direct current is injected
        @test hbnlsolve(ws, (1,), acsrc, c, Dict{Any,Any}();
            keyedarrays = false, dc = true, odd = true,
            returnoperatingpoint = true).solverinfo.converged
    end
end
