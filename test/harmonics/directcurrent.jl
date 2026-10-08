using JosephsonCircuits
using LinearAlgebra
using Random
using SparseArrays
using Test

# Direct current through resistors. The harmonic balance state is periodic
# node flux, so a voltage is its time derivative and vanishes at zero
# frequency: on the flux alone a resistor is an open circuit at DC and a
# current source driving one cannot develop I*R. The missing coordinate is
# the average voltage, which is constant on each static flux component
# because an inductor or a zero-voltage junction is a short there.

# an inner preconditioner which does nothing, so what is asserted below is
# the direct current step and not the mode coupling solve around it
struct Passthrough <: JosephsonCircuits.AbstractPreconditioner end
JosephsonCircuits.applypreconditioner!(z, ::Passthrough, r) = copyto!(z, r)
JosephsonCircuits.updatepreconditioner!(pc::Passthrough, x) = pc

# a window whose rows marked `unread` throw when read, on the host backend
struct UnreadRows <: AbstractVector{Float64}
    x::Vector{Float64}
    unread::BitVector
end
Base.size(v::UnreadRows) = size(v.x)
Base.getindex(v::UnreadRows, i::Int) = v.unread[i] ? error("row $i was read") : v.x[i]
Base.setindex!(v::UnreadRows, a, i::Int) = (v.x[i] = a)
JosephsonCircuits.KernelAbstractions.get_backend(::UnreadRows) =
    JosephsonCircuits.KernelAbstractions.CPU()

# a grounded two port block which is a series impedance `Z` at every
# frequency, an ideal through when `Z` is zero; the tests build their
# series blocks here, so that the blocks share one type and their
# evaluation is compiled once
seriesblock(Z) = ScatteringParameters(
    w -> JosephsonCircuits.ABCDtoS(JosephsonCircuits.ABCD_seriesZ(Z + 0im));
    nports = 2, grounded = true, noise = Lossless())

@testset verbose=true "direct current through resistors" begin

    ws = (2*pi*5e9,)
    dcsolve(c, sources) = hbnlsolve(ws, (1,), sources, c, Dict{Any,Any}();
        keyedarrays = false, dc = true, odd = true)

    @testset "a constant current source of the netlist" begin
        # a junction biased by a current source component: the phase of
        # the junction is the analytic one, the same source at a port
        # gives the same static flux, and the transient settles to it
        JC = JosephsonCircuits
        Lj, I = 1e-9, 1e-7
        phase = asin(I*Lj/JC.phi0)
        probe = [(mode = (1,), port = 1, current = 1e-13)]
        dcmode(out) = findfirst(m -> all(iszero, m), out.modes)
        biased = Circuit([(:p1, 1, 0, Port(1; termination = nothing)),
            (:jj, 1, 0, JosephsonJunction(Lj)), (:c, 1, 0, Capacitor(1e-12)),
            (:i, 0, 1, CurrentSource(I))])
        out = dcsolve(biased, probe)
        @test out.solverinfo.converged
        @test out.nodeflux[dcmode(out)] ≈ phase rtol=1e-8
        ported = Circuit([(:p1, 1, 0, Port(1; termination = nothing)),
            (:jj, 1, 0, JosephsonJunction(Lj)), (:c, 1, 0, Capacitor(1e-12))])
        atport = dcsolve(ported, vcat(probe, [(mode = (0,), port = 1, current = I)]))
        @test out.nodeflux ≈ atport.nodeflux rtol=1e-10
        # the current leaves the source's first terminal: with the terminals
        # the other way round it leaves the node, and the phase reverses
        reversed = Circuit([(:p1, 1, 0, Port(1; termination = nothing)),
            (:jj, 1, 0, JosephsonJunction(Lj)), (:c, 1, 0, Capacitor(1e-12)),
            (:i, 1, 0, CurrentSource(I))])
        rev = dcsolve(reversed, probe)
        @test rev.nodeflux[dcmode(rev)] ≈ -phase rtol=1e-8
        # the transient reads the same source, and settles to the same flux
        damped = Circuit([(:p1, 1, 0, Port(1)), (:l1, 1, 2, Inductor(1e-9)),
            (:jj, 2, 0, JosephsonJunction(Lj)), (:c2, 2, 0, Capacitor(1e-12)),
            (:i, 0, 2, CurrentSource(I))])
        hb = dcsolve(damped, probe)
        nm = length(hb.modes)
        sol = transientsolve(transientproblem(damped), (0.0, 5e-9); dt = 5e-12)
        @test sol.finalflux[2] ≈ hb.nodeflux[(2 - 1)*nm + dcmode(hb)] rtol=1e-6
        @test sol.finalflux[2] ≈ phase rtol=1e-6
        # without the zero frequency mode a nonzero source is an error, a
        # zero one is not, and a symbolic value is read from the definitions
        @test_throws ArgumentError hbnlsolve(ws, (1,), probe, biased; keyedarrays = false)
        @test hbnlsolve(ws, (1,), probe, Circuit([(:p1, 1, 0, Port(1)), (:jj, 1, 0, JosephsonJunction(Lj)),
            (:c, 1, 0, Capacitor(1e-12)), (:i, 0, 1, CurrentSource(0.0))]); keyedarrays = false).solverinfo.converged
        symbolic = Circuit([(:p1, 1, 0, Port(1; termination = nothing)),
            (:jj, 1, 0, JosephsonJunction(Lj)), (:c, 1, 0, Capacitor(1e-12)),
            (:i, 0, 1, CurrentSource(:Ibias))])
        sym = hbnlsolve(ws, (1,), probe, symbolic, Dict(:Ibias => I); keyedarrays = false, dc = true, odd = true)
        @test sym.nodeflux ≈ out.nodeflux rtol=1e-10
    end

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
        for m in (Newton(), NewtonKrylov())
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
        # the common case: the block is classified and not carried, the
        # answer is that of the flux alone, and the average voltage it
        # reports is the zero every node sits at
        c = Circuit([:p1 => Port(1), :cc => Capacitor(100e-15),
                     :jj => JosephsonJunction(1000e-12),
                     :cj => Capacitor(1000e-15)],
            [[(:p1,1),(:cc,1)], [(:cc,2),(:jj,1),(:cj,1)],
             [(:p1,2),(:jj,2),(:cj,2), Ground]])
        sol = hbnlsolve((2*pi*4.75001e9,), (8,),
            [(mode=(1,), port=1, current=0.00565e-6)], c, Dict{Any,Any}();
            keyedarrays = false, dc = true)
        @test sol.solverinfo.converged
        @test sol.dcnodevoltage == zeros(2)
        # and without a zero frequency mode there is no average voltage
        sol = hbnlsolve((2*pi*4.75001e9,), (8,),
            [(mode=(1,), port=1, current=0.00565e-6)], c, Dict{Any,Any}();
            keyedarrays = false, dc = false)
        @test sol.solverinfo.converged
        @test isnothing(sol.dcnodevoltage)
    end

    @testset "the network is classified whether or not it is driven" begin
        # Whether a direct current is determined is a property of the
        # network and not of the drive: a short or an ideal through in
        # parallel with an inductive path leaves the division undetermined
        # with a drive and without one, and the answer the stamp's `i = 0`
        # row would give is one of infinitely many either way. So the
        # descriptor is assembled and classified whenever there is a zero
        # frequency mode, and the block is carried through the solve only
        # when something is injected.
        JC = JosephsonCircuits
        R = 100.0
        through() = seriesblock(0.0)
        finite() = seriesblock(10.0)
        pump = [(mode=(1,), port=1, current=1e-6)]
        direct = [(mode=(0,), port=1, current=1e-6)]
        go(c, srcs) = hbnlsolve(ws, (1,), srcs, c, Dict{Any,Any}();
            keyedarrays = false, dc = true, odd = true)

        # the loop which is refused with a direct current drive is the same
        # loop without one
        loop(t) = Circuit(
            [:p1 => Port(1; Z0 = R), :l => Inductor(1e-9), :t => t,
             :c1 => Capacitor(1e-12), :c2 => Capacitor(1e-12)],
            [[(:p1,1),(:l,1),(:t,1),(:c1,1)], [(:l,2),(:t,2),(:c2,1)],
             [(:p1,2),(:c1,2),(:c2,2), Ground]])
        @test_throws ArgumentError go(loop(through()), pump)
        @test go(loop(finite()), pump).solverinfo.converged

        # every node held at ground by an inductor: there is no voltage to
        # find, and there is still a current to classify. The classification
        # runs with no floating component in the circuit, or this loop would
        # solve with all of its direct current through the inductors.
        grounded(t) = Circuit(
            [:p1 => Port(1; Z0 = R), :lg => Inductor(1e-9),
             :l => Inductor(1e-9), :t => t,
             :c1 => Capacitor(1e-12), :c2 => Capacitor(1e-12)],
            [[(:p1,1),(:lg,1),(:l,1),(:t,1),(:c1,1)],
             [(:l,2),(:t,2),(:c2,1)],
             [(:p1,2),(:lg,2),(:c1,2),(:c2,2), Ground]])
        @test_throws ArgumentError go(grounded(through()), pump)
        @test_throws ArgumentError go(grounded(through()), direct)
        # a finite block between two nodes held at ground carries nothing,
        # and the whole drive goes to ground through the inductor
        s = go(grounded(finite()), direct)
        @test s.solverinfo.converged
        @test s.dcnodevoltage == zeros(2)

        # a block whose data does not reach zero frequency is open there
        # when nothing asks it to carry direct current, and has to state its
        # limit when something does
        f = 2*pi*(1e9:1e9:10e9)
        S = zeros(Complex{Float64}, 2, 2, length(f))
        for i in eachindex(f)
            S[:,:,i] = JC.ABCDtoS(JC.ABCD_seriesZ(10.0 + 0im))
        end
        tab = ScatteringParameters((collect(f), S); nports = 2,
            grounded = true)
        c = Circuit(
            [:p1 => Port(1; Z0 = R), :x => tab, :r2 => Resistor(R),
             :c1 => Capacitor(1e-12)],
            [[(:p1,1),(:x,1),(:c1,1)], [(:x,2),(:r2,1)],
             [(:p1,2),(:r2,2),(:c1,2), Ground]])
        @test go(c, pump).solverinfo.converged
        @test_throws ArgumentError go(c, direct)

        # Skipping the block when nothing is injected is exact and not an
        # approximation: a nonlinear inductor whose current-phase relation
        # has an even term rectifies a single tone, and a floating island
        # joined to the rest by resistors could in principle develop an
        # average voltage from it. It does not, because a rectified current
        # circulates inside the island's inductive paths and the transport
        # row sums it away. The same circuit driven by a direct current too
        # small to matter carries the block explicitly, and lands on the
        # same point.
        j = Circuit(
            [:p1 => Port(1), :cc => Capacitor(100e-15), :rb => Resistor(200.0),
             :jj => NonlinearInductor(500e-12, PolynomialCPR([1.0, 0.3, -1/6])),
             :cj => Capacitor(500e-15),
             :l3 => Inductor(2e-9), :c3 => Capacitor(1e-12),
             :r3 => Resistor(300.0), :c4 => Capacitor(1e-12)],
            [[(:p1,1),(:cc,1)],
             [(:cc,2),(:rb,1),(:jj,1),(:cj,1),(:l3,1)],
             [(:rb,2),(:jj,2),(:cj,2),(:c3,1),(:l3,2)],
             [(:c3,2),(:r3,1)], [(:r3,2),(:c4,1)],
             [(:p1,2),(:c4,2), Ground]])
        ws1 = (2*pi*4.75e9,)
        tone = [(mode=(1,), port=1, current=1.2e-6)]
        kw = (; keyedarrays = false, dc = true, odd = true, even = true,
            atol = 1e-11)
        a = hbnlsolve(ws1, (8,), tone, j, Dict{Any,Any}(); kw...)
        @test a.solverinfo.converged
        @test a.dcnodevoltage == zeros(5)
        # the tone does rectify: the island carries a zero frequency flux
        k0 = findfirst(m -> all(iszero, m), a.frequencies.modes)
        dcflux = reshape(a.nodeflux, length(a.frequencies.modes), :)[k0, :]
        @test maximum(abs, dcflux) > 1e-3
        b = hbnlsolve(ws1, (8,),
            vcat(tone, [(mode=(0,), port=1, current=1e-30)]), j,
            Dict{Any,Any}(); kw...)
        @test b.solverinfo.converged
        @test maximum(abs, b.dcnodevoltage) < 1e-20
        @test isapprox(a.nodeflux, b.nodeflux;
            atol = 1e-12*maximum(abs, a.nodeflux))
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
            v = pinv(Matrix(t.Y)) * t.j

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
        @test isfinite(cond(Matrix(t.Y)))

        # a floating island, where only differences are physical: the rows
        # say exactly that, and are singular by one direction
        Z0 = 200.0
        cf = Circuit(
            [:p1 => Port(1; Z0 = Z0), :ca => Capacitor(1e-12),
             :cb => Capacitor(1e-12)],
            [[(:p1,1),(:ca,1)], [(:p1,2),(:cb,1)],
             [(:ca,2),(:cb,2), Ground]])
        t = agrees(cf, [(mode=(0,), port=1, current=Idc)])
        @test rank(Matrix(t.Y)) == size(t.Y, 1) - 1
        @test isapprox(t.Y * ones(size(t.Y, 2)), zeros(size(t.Y, 1));
            atol = 1e-12*maximum(abs, t.Y))
    end

    # The capability the explicit block exists for. With an average voltage
    # to respond to, a scattering block obeys its own relation at direct
    # current rather than the open circuit row `i = 0`, and a block which
    # is a resistor carries what that resistor would.
    @testset "a scattering block carries direct current" begin
        Rb, Idc, Zbig = 100.0, 1.0e-6, 1.0e9

        # the same circuit twice: once with the resistor, once with a
        # scattering block that is the resistor
        asres = Circuit(
            [:p1 => Port(1; Z0 = Zbig), :rb => Resistor(Rb),
             :c1 => Capacitor(1e-12)],
            [[(:p1,1),(:rb,1),(:c1,1)], [(:rb,2), Ground],
             [(:p1,2),(:c1,2), Ground]])
        blk = seriesblock(Rb)
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

        # were the block an open, the whole current would go through the
        # port environment instead and the voltage would be Idc*Zbig, seven
        # orders larger, so the agreement above is not an accident

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
        # conductance graph, so it is the block's own row which makes this
        # circuit solvable: the driven node sits at I*R and the shorted one
        # at zero.
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
        through() = seriesblock(0.0)
        undetermined = Circuit(
            [:p1 => Port(1; Z0 = R), :l => Inductor(1e-9), :t => through(),
             :c1 => Capacitor(1e-12), :c2 => Capacitor(1e-12)],
            [[(:p1,1),(:l,1),(:t,1),(:c1,1)], [(:l,2),(:t,2),(:c2,1)],
             [(:p1,2),(:c1,2),(:c2,2), Ground]])
        @test_throws ArgumentError go(undetermined)

        # giving the block a finite series impedance determines the division
        # and the same circuit solves
        finite() = seriesblock(10.0)
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
        through = seriesblock(0.0)
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
        through = seriesblock(0.0)
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
        go(I; kw...) = hbnlsolve(ws, (4,), [(mode=(0,), port=1, current=I)], c,
            Dict{Any,Any}(); keyedarrays = false, dc = true, odd = true,
            even = true, kw...)

        # well inside the zero voltage state: nothing to say
        a = @test_logs go(0.3*Ic)
        @test a.solverinfo.converged

        # at the edge of it: named, with the fraction it reached
        b = (@test_logs (:warn,) go(0.996*Ic))
        @test b.solverinfo.converged

        # and past it there is no periodic solution to find
        # (the search cannot succeed, so its budget is what it costs)
        d = (@test_logs (:warn,) match_mode = :any go(1.03*Ic; iterations = 40))
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
        # to the same tolerance: a dissipation `I - S S'` of -0.014,
        # below -atol, though the singular value is within 1 + atol
        @test_throws ArgumentError ScatteringParameters([0 1;1 0]; atol = 1e-2,
            dcmodel = ScatteringDC([0 1.007; 1.007 0]))
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

    # The preconditioner solves the direct current subsystem exactly and
    # writes the answer over whatever the inner preconditioner guessed at
    # those coordinates. On a device a subsystem this small is solved there,
    # through dense factors of the same matrix, rather than the window
    # crossing the bus; here the host path is checked against the sparse
    # factorization it is meant to reproduce.
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
        r = randn(n)
        z = zeros(n)
        JC.applypreconditioner!(z, pc, r)

        # the flux coordinates carry the inner preconditioner's answer
        idx = JC.dcsubsystemindices(w)
        rest = setdiff(1:n, idx)
        @test z[rest] == r[rest]
        # and the subsystem's own are its exact solution
        A = JC.dcsubsystem(w)
        @test z[idx] == JC.kluordered(A) \ r[idx]
        @test A*z[idx] ≈ r[idx]
        # which is a real solve and not the identity: this subsystem is ill
        # conditioned enough that an inverse would not do
        @test cond(Matrix(A)) > 1e6
        @test z[idx] != r[idx]

        # The block in its matrix form, which a device applies, agrees with
        # the scalar form where the rows it writes over hold a NaN: the
        # internal residual does not write them, and on a device they hold
        # whatever the memory did.
        up = JC.dcupdate(w)
        nw = length(up.keep)
        uw, Fw = randn(nw), randn(nw)
        Fw[iszero.(up.keep)] .= NaN
        @test any(iszero, up.keep)
        @test JC.applydcupdate!(copy(Fw), uw, up) ≈ JC.addtransportwindow!(copy(Fw), uw, w)
        # and does not read them at all
        Fu = UnreadRows(copy(Fw), iszero.(up.keep))
        JC.applydcupdate!(Fu, uw, up)
        @test Fu.x ≈ JC.addtransportwindow!(copy(Fw), uw, w)
    end

    @testset "the held point is recognized through a view of the canonical state" begin
        # the solve hands the system views of its canonical state: the
        # residual sets the point through one, and the preconditioner's
        # update sets it again through another. That point is the one the
        # system holds, so what was evaluated there stays, and a point which
        # moved is set.
        JC = JosephsonCircuits
        c = Circuit([(:p1, 1, 0, Port(1)), (:r, 1, 2, Resistor(20.0)),
            (:jj, 2, 0, JosephsonJunction(1e-9)),
            (:c2, 2, 0, Capacitor(1e-12))])
        d = hbnlsolve(ws, (1,), [(mode = (1,), port = 1, current = 2e-6),
            (mode = (0,), port = 1, current = 1e-7)], c, Dict{Any,Any}();
            keyedarrays = false, dc = true, odd = true, returnsystem = true)
        sys, L = d.sys, d.canonicalwork.layout
        u = 0.1 .* sin.(1:JC.canonicaldim(L))
        JC.residual!(zeros(L.rdim), JC.setpoint!(sys, JC.internalpart(u, L)))
        JC.cosphimatrix(sys)
        JC.setpoint!(sys, JC.internalpart(copy(u), L))
        @test sys.sincurrent[] && sys.cosfdcurrent[]
        u[1] += 0.01
        JC.setpoint!(sys, JC.internalpart(u, L))
        @test !sys.sincurrent[] && sys.xr == JC.internalpart(u, L)
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
            rtol = 1e-13, method = Newton(), kw...)

        sol = go(R2; returnoperatingpoint = true)
        op = sol.operatingpoint
        @test !isnothing(op.dc)
        L = op.dc.work.layout
        @test length(op.dc.u) == JC.canonicaldim(L)

        psc = JC.compile(mk(R2))
        nm = numericmatrices(psc, Dict{Any,Any}(); Nmodes = op.Nmodes)
        ir2 = findfirst(==("r2"), psc.componentnames)
        dFr = JC.calcresidualsensitivity(op, psc, nm, [ir2])
        # the residual derivative is in the canonical coordinates, and it
        # reaches the transport rows, which the harmonic one cannot
        @test size(dFr, 1) == JC.canonicaldim(L)
        @test !iszero(dFr[end, 1])

        dx = JC.calcnodefluxsensitivity(op, dFr)

        # against a central difference of a re-solve, in the same relative
        # parameter the sensitivity is taken in
        h = 1e-6
        sp = go(R2*(1+h)); sm = go(R2*(1-h))
        fdflux = (sp.nodeflux .- sm.nodeflux)./(2h)
        @test isapprox(dx, fdflux; rtol = 1e-4,
            atol = 1e-8*maximum(abs, fdflux))
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
    # A current source into a resistive island which reaches the rest of
    # the circuit only through capacitors injects a net direct current the
    # island cannot carry away, and the solve refuses it. Both cases are
    # then taken on a doctored subsystem, where the constant side can be
    # put along either direction of a singular block.
    @testset "an unsolvable direct current block is refused" begin
        JC = JosephsonCircuits

        island = Circuit([(:p1, 1, 0, Port(1)), (:cp, 1, 2, Capacitor(1e-12)),
            (:c1, 2, 0, Capacitor(1e-12)), (:r, 2, 3, Resistor(100.0)),
            (:c2, 3, 0, Capacitor(1e-12)), (:i, 0, 2, CurrentSource(1e-6))])
        @test_throws ArgumentError hbnlsolve(ws, (1,),
            [(mode = (1,), port = 1, current = 1e-13)], island,
            Dict{Any,Any}(); keyedarrays = false, dc = true, odd = true)

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
        singular = sparse([1.0 -1.0; -1.0 1.0])
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
    # it to be differentiated later, would describe a different problem from
    # the one solved. `returnsystem` and `debugJacobian` hand back the
    # canonical work beside it, and the operating point carries the block.
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


        # and both are available when no direct current is injected
        @test hbnlsolve(ws, (1,), acsrc, c, Dict{Any,Any}();
            keyedarrays = false, dc = true, odd = true,
            returnoperatingpoint = true).solverinfo.converged
    end

    # The direct current subsystem has one unknown per floating island and
    # one per block port current, so a chain of capacitively coupled
    # islands, or a cascade of blocks, makes it as long as the circuit. It
    # is assembled, classified and factorized sparse, block by block of its
    # pattern, so a fourfold step allocates four times as much, not the
    # sixteen or more of a dense subsystem's square or cube. The block's
    # matrix form, which the canonical Jacobian and a device read off the
    # residual, takes as many passes over it for the long chains as for the
    # short ones.
    @testset "a long chain of islands or of blocks" begin
        islandchain(n) = Circuit(vcat(
            Any[(:p1, 1, 0, Port(1; Z0 = 50.0)), (:r2, 2n, 0, Resistor(50.0))],
            [(Symbol(:j, k), 2k - 1, 2k, JosephsonJunction(100e-12)) for k in 1:n],
            [(Symbol(:ca, k), 2k - 1, 0, Capacitor(40e-15)) for k in 1:n],
            [(Symbol(:cb, k), 2k, 0, Capacitor(40e-15)) for k in 1:n],
            [(Symbol(:cc, k), 2k, 2k + 1, Capacitor(1e-12)) for k in 1:n-1]))
        blockchain(n) = Circuit(vcat(
            Any[(:p1, 1, 0, Port(1; Z0 = 50.0)), (:r2, 2n + 1, 0, Resistor(50.0))],
            [(Symbol(:t, k), 2k - 1, 2k, seriesblock(10.0)) for k in 1:n],
            [(Symbol(:j, k), 2k, 2k + 1, JosephsonJunction(100e-12)) for k in 1:n],
            [(Symbol(:c, k), 2k + 1, 0, Capacitor(40e-15)) for k in 1:n]))
        # the bytes of the setup after a warm call, and the passes the
        # block's matrix form takes
        function measure(make, n)
            c = JosephsonCircuits.compile(make(n))
            run() = hbnlsolve(ws, (1,), [(mode = (0,), port = 1, current = 1e-7)],
                c, Dict{Any,Any}(); keyedarrays = false, dc = true, odd = true,
                returnsystem = true)
            passes = length(first(JosephsonCircuits.dcprobes(run().canonicalwork)))
            return @allocated(run()), passes
        end
        for make in (islandchain, blockchain)
            (b64, p64), (b256, p256) = measure(make, 64), measure(make, 256)
            @test b256 < 6*b64
            @test p256 == p64
        end
    end

    # A floating chain of resistors longer than the blocks the subsystem
    # decomposes by their singular values, so it is classified by sparse QR.
    # Its common voltage is a gauge and is pinned: driven across a port of
    # its own, the chain carries the drop of its resistance in parallel with
    # the port's environment. Driven at one end alone, the current has
    # nowhere to go, and it is refused.
    @testset "a long floating chain is classified by sparse QR" begin
        m, R, Z0, Idc = 80, 10.0, 50.0, 1.0e-6
        chain(extra...) = Circuit(vcat(Any[(:p1, 1, m, Port(1; Z0 = Z0))],
            [(Symbol(:r, k), k, k + 1, Resistor(R)) for k in 1:m-1],
            [(Symbol(:c, k), k, 0, Capacitor(1e-12)) for k in 1:m], Any[extra...]))
        sol = dcsolve(chain(), [(mode = (0,), port = 1, current = Idc)])
        @test sol.solverinfo.converged
        v = sol.dcnodevoltage
        Rc = (m - 1)*R
        @test isapprox(maximum(v) - minimum(v), Idc*Z0*Rc/(Z0 + Rc); rtol = 1e-9)
        @test_throws ArgumentError dcsolve(chain((:i, 0, 1, CurrentSource(Idc))),
            [(mode = (1,), port = 1, current = 1e-13)])
    end

    # The rank of the direct current subsystem is decided in the circuit's
    # units. In the solver's an average voltage is V/phi0 and a block current
    # I Lscale/phi0, and equilibrating the matrix as it stands leaves a row
    # which mixes the two far out of balance, so that a well posed subsystem
    # can be taken for a singular one. A resistor between a port and a
    # resistive one port block: its subsystem, put in the circuit's units by
    # the plan's voltage scale and equilibrated, has the singular values of
    # the same equations written from the element values in volts and in
    # amperes times the reference impedance. And a circuit of bias resistors
    # from 114 ohm to 260 Mohm around two near-through sections, nonsingular,
    # which was refused: a current driven through it develops its I*R.
    @testset "the rank is decided in the circuit's units" begin
        JC = JosephsonCircuits
        Zp, R, Rb, z = 50.0, 30.0, 80.0, 50.0
        load = ScatteringParameters(fill((Rb - z)/(Rb + z) + 0im, 1, 1); zref = z)
        c = Circuit([(:p1, 1, 0, Port(1; Z0 = Zp)), (:c1, 1, 0, Capacitor(1e-12)),
            (:c2, 2, 0, Capacitor(1e-12)), (:r, 1, 2, Resistor(R)), (:b, 2, load)])
        w = hbnlsolve(ws, (1,), [(mode = (0,), port = 1, current = 1e-9)], c,
            Dict{Any,Any}(); keyedarrays = false, dc = true, odd = true,
            returnsystem = true).canonicalwork
        B, _, _ = JC.equilibrate(JC.dcsubsystem(w), JC.nvoltages(w.transport),
            w.transport.plan.voltagescale)
        # Kirchhoff's law at the two nodes and the block's relation
        # (1 - S) V/sqrt(z) - (1 + S) sqrt(z) I = 0, in V1, V2 and Z0*I, with
        # Z0 the reference impedance, here the port's
        S, Z0 = (Rb - z)/(Rb + z), Zp
        M = [Z0/Zp + Z0/R  -Z0/R  0.0;
             -Z0/R  Z0/R  1.0;
             0.0  (1 - S)/sqrt(z)  -(1 + S)*sqrt(z)/Z0]
        Mb, _, _ = JC.equilibrate(sparse(M), 0, 1.0)
        @test svdvals(Matrix(B)) ≈ svdvals(Matrix(Mb)) rtol = 1e-10

        I = 1e-9
        bias = Circuit(vcat(Any[(:p1, 1, 0, Port(1; Z0 = 50.0))],
            [(Symbol(:c, k), k, 0, Capacitor(1e-12)) for k in 1:6],
            Any[(:r7, 2, 3, Resistor(34185.0)), (:r8, 2, 6, Resistor(800733.0)),
            (:t9, 4, 2, seriesblock(0.019)), (:t10, 2, 5, seriesblock(0.038)),
            (:r11, 6, 1, Resistor(113.74)), (:r12, 4, 2, Resistor(2.5752e8)),
            (:i5, 0, 5, CurrentSource(I))]))
        sol = dcsolve(bias, [(mode = (1,), port = 1, current = 1e-13)])
        @test sol.solverinfo.converged
        # the source's current has one way to ground: through t10, r8, r11
        # and the port
        names = JC.compile(bias).nodenames
        v = Dict(n => sol.dcnodevoltage[k - 1] for (k, n) in enumerate(names) if k > 1)
        @test v["2"] - v["6"] ≈ I*800733.0 rtol = 1e-8
    end
end
