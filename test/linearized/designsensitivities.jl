using JosephsonCircuits
using LinearAlgebra
using SparseArrays
using Test

# The design Jacobian of a typed circuit and the sensitivities built on it,
# against finite differences of the full solve.
@testset verbose=true "design sensitivities" begin

    wp = (2*pi*4.75001*1e9,)
    src = [(mode=(1,), port=1, current=0.00565e-6)]
    ws = 2*pi*(4.5:0.1:5.0)*1e9
    # finite differences of the scattering parameters at definitions moved
    # by a relative step in one parameter
    function fdS(circuit, defs, q; h = 1e-6, kw...)
        at(v) = merge(Dict{Any,Any}(defs), Dict{Any,Any}(q => v))
        x = defs[q]
        Sp = hbsolve(ws, wp, src, (2,), (8,), circuit, at(x*(1 + h)); kw...).linearized.S
        Sm = hbsolve(ws, wp, src, (2,), (8,), circuit, at(x*(1 - h)); kw...).linearized.S
        return vec((Array(Sp) .- Array(Sm))./(2*h*x))
    end

    @testset "the definitions are normalized once for the Jacobian" begin
        # an expression valued component resolves its derivative against
        # the definitions normalized once for the whole Jacobian, so an
        # expression costs what the symbol spelling costs
        CV = JosephsonCircuits.CircuitValues
        n = 64
        pars = [CV.Parameter(Symbol(:c, i)) for i in 1:n]
        expr = compile(Circuit(vcat(Any[(:P1, 1, 0, Port(1))],
            Any[(Symbol(:C, i), 1, 0, Capacitor(pars[i])) for i in 1:n])))
        syms = compile(Circuit(vcat(Any[(:P1, 1, 0, Port(1))],
            Any[(Symbol(:C, i), 1, 0, Capacitor(Symbol(:c, i))) for i in 1:n])))
        edefs = Dict{Any,Any}(pars[i] => 1e-13*i for i in 1:n)
        sdefs = Dict{Any,Any}(Symbol(:c, i) => 1e-13*i for i in 1:n)
        je() = JosephsonCircuits.designjacobian(expr, edefs)
        js() = JosephsonCircuits.designjacobian(syms, sdefs)
        je(); js()
        @test je()[3] == js()[3]
        @test @allocated(je()) < 2*@allocated(js())
    end

    @testset "the derivative of an expression" begin
        JosephsonCircuits.@params a b
        CV = JosephsonCircuits.CircuitValues
        d(e, n, defs) = JosephsonCircuits.designderivative(e, n, defs)
        defs = Dict(a => 2.0, b => 3.0)
        @test d(a, :a, defs) == 1 && d(a, :b, defs) == 0
        @test d(3.0, :a, defs) == 0 && d(:a, :a, defs) == 1 && d("b", :b, defs) == 1
        @test d(a*b, :a, defs) == 3 && d(a*b, :b, defs) == 2
        @test d(a/b, :b, defs) ≈ -2/9
        @test d(a^3, :a, defs) ≈ 12
        @test d(2^a, :a, defs) ≈ 4*log(2)
        @test d(sqrt(a), :a, defs) ≈ 1/(2*sqrt(2))
        @test d(exp(a*b), :b, defs) ≈ 2*exp(6)
        @test d(log(a), :a, defs) ≈ 0.5
        @test d(inv(a), :a, defs) ≈ -0.25
        @test d(1/(1 + im*a), :a, defs) ≈ -im/(1 + 2im)^2
        # a value which does not depend on the parameter collapses to zero
        @test CV.derivative(CV.tocv(b) + 1, :a) == CV.Constant(0)
        # an undefined parameter in the derivative is refused, as is a
        # frequency dependent value
        @test_throws ArgumentError d(a*b, :a, Dict(a => 2.0))
        @test_throws ArgumentError d(FrequencyDependent(w -> 1.0), :a, defs)
    end

    @testset "designjacobian" begin
        JosephsonCircuits.@params L C
        circuit = Circuit(Any[
            ("P1", "1", "0", Port(1; Z0 = 50.0)), ("C1", "1", "2", Capacitor(C)),
            ("Lj1", "2", "0", JosephsonJunction(L)), ("C2", "2", "0", Capacitor(C/4))])
        defs = Dict(L => 1000.0e-12, C => 400.0e-15)
        names, v0, J = JosephsonCircuits.designjacobian(circuit, defs)
        # the port is not among them: its reference impedance appears under
        # the termination generated for it
        @test names == ["P1/termination","C1","Lj1","C2"]
        @test v0 == [50.0, 400e-15, 1000e-12, 100e-15]
        # the parameters default to the defined ones, sorted by name, and
        # the rows are exact: dC1/dC = 1, dC2/dC = 1/4, dLj1/dL = 1
        iC1 = findfirst(==("C1"), names); iC2 = findfirst(==("C2"), names)
        iL = findfirst(==("Lj1"), names); iR = findfirst(==("P1/termination"), names)
        jC, jL = 1, 2
        @test J[iC1, jC] == 1.0 && J[iC2, jC] == 0.25 && J[iL, jL] == 1.0
        @test iszero(J[iR, jL]) && iszero(J[iR, jC]) && iszero(J[iC1, jL])
        # the parameters may be named, in the order given, by name or key
        _, _, JL = JosephsonCircuits.designjacobian(circuit, defs; parameters = (:L,))
        @test size(JL) == (4, 1) && JL[:, 1] == J[:, jL]
        _, _, JK = JosephsonCircuits.designjacobian(circuit, defs; parameters = (L, C))
        @test JK == J[:, [jL, jC]]
        # a parameter left undefined is an error, not a garbage derivative
        @test_throws ArgumentError JosephsonCircuits.designjacobian(circuit, Dict(L => 1e-9))
    end

    @testset "designsensitivities matches full-solve finite differences" begin
        Lj, Cc, Cj = JosephsonCircuits.@params Lj Cc Cj
        jpa = Circuit(Any[
            ("P1", "1", "0", Port(1; Z0 = 50.0)), ("C1", "1", "2", Capacitor(Cc)),
            ("Lj1", "2", "0", JosephsonJunction(Lj)), ("C2", "2", "0", Capacitor(Cj)),
            ("C3", "2", "0", Capacitor(Cj/4))])     # derived value: chain-rule factor 1/4
        defs = Dict(Lj => 1000.0e-12, Cc => 100.0e-15, Cj => 800.0e-15)
        r = designsensitivities(jpa, defs, ws, wp, src, (2,), (8,))
        @test size(r.dSdp, 5) == 3
        for q in (Cj, Lj)
            fd = fdS(jpa, defs, q)
            mine = vec(Array(r.dSdp(parameter = q.name)))
            @test norm(mine .- fd)/norm(fd) < 1e-4
        end
        # a parameter no component depends on is rejected
        @test_throws ArgumentError designsensitivities(jpa,
            merge(defs, Dict(:unused => 1.0)), ws, wp, src, (2,), (8,);
            parameters = (:unused,))
    end

    @testset "three tones, direct current and negative sidebands" begin
        # the chain rule has to carry every retained mode of a mixed grid:
        # three tones with both mixing orders and direct current, whose
        # signal modes include ones at negative physical frequency. Forward
        # and reverse must agree, and both must agree with differences of
        # the whole solve.
        JosephsonCircuits.@params scale loss
        c = Circuit([(:p, 1, 0, Port(1)),
            (:cc, 1, 2, Capacitor(100e-15/(1 + im*loss))),
            (:jj, 2, 0, JosephsonJunction(1e-9/scale)),
            (:cj, 2, 0, Capacitor(1e-12*scale))])
        defs = Dict(scale => 1.0, loss => 0.02)
        pumps = 2*pi .* (4.75001e9, 1.17003e9, 0.63007e9)
        drive = [(mode = (1,0,0), port = 1, current = 1e-8),
            (mode = (0,1,0), port = 1, current = 5e-9),
            (mode = (0,0,1), port = 1, current = 3e-9),
            (mode = (0,0,0), port = 1, current = 1e-9)]
        signals = 2*pi*[4.41e9, 4.59e9]
        opts = (; dc = true, threewavemixing = true, fourwavemixing = true,
            atol = 1e-12)
        forward = designsensitivities(c, defs, signals, pumps, drive,
            (1,1,1), (2,1,1); sensitivitymode = :forward, opts...)
        reverse = designsensitivities(c, defs, signals, pumps, drive,
            (1,1,1), (2,1,1); sensitivitymode = :reverse, opts...)
        @test forward.out.nonlinear.solverinfo.converged
        @test reverse.out.nonlinear.solverinfo.converged
        @test isapprox(Array(forward.dSdp), Array(reverse.dSdp), rtol = 1e-8)
        # a retained signal mode sits below zero frequency
        @test any(signals[1] + sum(m .* pumps) < 0
            for m in forward.out.linearized.modes)
        for parameter in (scale, loss)
            h = 1e-5*defs[parameter]
            at(v) = merge(defs, Dict(parameter => v))
            # the two reference solves are checked before their scattering
            # matrices are differenced, so a difference which disagrees with
            # the sensitivity is a derivative failure and not an unconverged
            # operating point at one of the perturbed points
            up = hbsolve(signals, pumps, drive, (1,1,1), (2,1,1), c,
                at(defs[parameter] + h); opts...)
            dn = hbsolve(signals, pumps, drive, (1,1,1), (2,1,1), c,
                at(defs[parameter] - h); opts...)
            @test up.nonlinear.solverinfo.converged
            @test dn.nonlinear.solverinfo.converged
            fd = (Array(up.linearized.S) .- Array(dn.linearized.S))./(2*h)
            d = Array(forward.dSdp(parameter = parameter.name))
            @test norm(d .- fd)/norm(fd) < 1e-4
        end
    end

    @testset "rotating directions and shared parameters are exact" begin
        # a design parameter which ROTATES a complex component value (a
        # loss tangent) has a direction which is not parallel to the value,
        # so a relative (logarithmic) derivative cannot represent it: the
        # solver must carry the direction dv/dp itself. Also here: one
        # parameter (Ic) touching two components of different kinds (the
        # junction inductance and its capacitance), which exercises the
        # per-parameter stamp merging and the mixed-kind grouping.
        phi0 = 3.29105976e-16
        Ic, adens, t = JosephsonCircuits.@params Ic adens t
        circuit = Circuit(Any[
            ("P1", "1", "0", Port(1; Z0 = 50.0)), ("C1", "1", "2", Capacitor(100e-15/(1 + im*t))),
            ("Lj1", "2", "0", JosephsonJunction(phi0/Ic)),
            ("C2", "2", "0", Capacitor(Ic*adens))])
        defs = Dict(Ic => phi0/1000e-12, adens => 1000e-15/(phi0/1000e-12), t => 2e-3)
        r = designsensitivities(circuit, defs, ws, wp, src, (2,), (8,))
        for q in (Ic, adens, t)
            fd = fdS(circuit, defs, q)
            mine = vec(Array(r.dSdp(parameter = q.name)))
            @test norm(mine .- fd)/norm(fd) < 1e-4
        end
    end

    @testset "reverse contraction agrees with forward" begin
        # the reverse contraction accumulates the operating point shift into
        # the same design parameter slots as the forward one
        Lj, Cc = JosephsonCircuits.@params Lj Cc
        circuit = Circuit(Any[
            ("P1", "1", "0", Port(1; Z0 = 50.0)), ("C1", "1", "2", Capacitor(Cc)),
            ("Lj1", "2", "0", JosephsonJunction(Lj)), ("C2", "2", "0", Capacitor(800e-15)),
            ("C3", "2", "3", Capacitor(100e-15)), ("Lj2", "3", "0", JosephsonJunction(2*Lj)),
            ("C4", "3", "0", Capacitor(800e-15))])
        defs = Dict(Lj => 1000.0e-12, Cc => 100.0e-15)
        rf = designsensitivities(circuit, defs, ws, wp, src, (2,), (8,);
            sensitivitymode = :forward)
        rr = designsensitivities(circuit, defs, ws, wp, src, (2,), (8,);
            sensitivitymode = :reverse)
        a = vec(Array(rf.dSdp)); b = vec(Array(rr.dSdp))
        @test norm(a .- b)/norm(a) < 1e-8
    end

    @testset "port impedance design parameter" begin
        circuit = Circuit(Any[
            ("P1", "1", "0", Port(1; Z0 = :R)), ("C1", "1", "2", Capacitor(:Cc)),
            ("Lj1", "2", "0", JosephsonJunction(1000e-12)), ("C2", "2", "0", Capacitor(1000e-15))])
        defs = Dict(:R => 50.0, :Cc => 100.0e-15)
        r = designsensitivities(circuit, defs, ws, wp, src, (2,), (8,))
        fd = fdS(circuit, defs, :R)
        mine = vec(Array(r.dSdp(parameter = :R)))
        @test norm(mine .- fd)/norm(fd) < 1e-4
    end

    @testset "scattering block design parameters" begin
        # a theta-dependent series-impedance two port feeding a JPA: the
        # block parameter and a lumped parameter share the output axis,
        # and the block moves the pump operating point. The block depends
        # on theta through the derivative it states.
        Z0 = 50.0
        function makeblk(theta; analytic = true)
            seriesS(w) = (z = 1/(im*w*theta*Z0);
                [z/(z+2) 2/(z+2); 2/(z+2) z/(z+2)])
            dS(w) = (z = 1/(im*w*theta*Z0); d = -2*z/(theta*(z+2)^2);
                [d -d; -d d])
            blk = analytic ?
                ScatteringParameters(seriesS; nports = 2, grounded = false,
                    derivatives = (theta = dS,)) :
                ScatteringParameters(seriesS; nports = 2, grounded = false)
            return Circuit(
                [:p1 => Port(1; termination = nothing), :r1 => Resistor(50.0), :cc => blk,
                 :jj => JosephsonJunction(:Lj),
                 :c2 => Capacitor(1000.0e-15)],
                [((:p1, 1), (:r1, 1), (:cc, 1, 1)),
                 ((:cc, 2, 1), (:jj, 1), (:c2, 1)),
                 ((:cc, 1, 2), (:cc, 2, 2), (:jj, 2), (:c2, 2), (:r1, 2),
                  (:p1, 2), Ground)])
        end
        theta = 100.0e-15
        defs = Dict(:Lj => 1000.0e-12)
        circuit = makeblk(theta)
        # the block's parameter is among the defaults, beside the defined one
        r = designsensitivities(circuit, defs, ws, wp, src, (2,), (8,))
        @test collect(JosephsonCircuits.AxisKeys.axiskeys(r.dSdp, :parameter)) == [:Lj, :theta]
        fdL = fdS(circuit, defs, :Lj)
        @test norm(vec(Array(r.dSdp(parameter = :Lj))) .- fdL)/norm(fdL) < 1e-4
        h = 1e-6*theta
        Sp = hbsolve(ws, wp, src, (2,), (8,), makeblk(theta + h), defs).linearized.S
        Sm = hbsolve(ws, wp, src, (2,), (8,), makeblk(theta - h), defs).linearized.S
        fdt = vec((Array(Sp) .- Array(Sm))./(2h))
        @test norm(vec(Array(r.dSdp(parameter = :theta))) .- fdt)/norm(fdt) < 1e-4

        # a block which states no derivative depends on no parameter and
        # costs nothing
        plain = makeblk(theta; analytic = false)
        @test isempty(JosephsonCircuits.designblockjacobian(plain, [:theta, :Lj]))
        rp = designsensitivities(plain, defs, ws, wp, src, (2,), (8,))
        @test collect(JosephsonCircuits.AxisKeys.axiskeys(rp.dSdp, :parameter)) == [:Lj]
        @test all(isfinite, Array(rp.dSdp))

        # the reverse contraction rejects block parameters
        @test_throws ArgumentError designsensitivities(circuit, defs, ws, wp,
            src, (2,), (8,); sensitivitymode = :reverse)
    end
end
