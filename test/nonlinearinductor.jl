using JosephsonCircuits
using LinearAlgebra
using Random
using Test

# The current-phase relation of a nonlinear inductor, given as a
# polynomial, and the solvers which evaluate it. The relation enters the
# solvers only as itself, its first derivative and its second, in place of
# `sin`, `cos` and `-sin`, so the checks here are that each of the three is
# the right one: the residual and the gain against the Josephson junction
# whose expansion the polynomial is, the Jacobian and the Hessian against
# finite differences, and a quadratic term against the three wave mixing it
# exists to produce.
@testset "the polynomial current-phase relation" begin
    JC = JosephsonCircuits
    L0 = 1e-9
    # the Taylor expansion of `sin` to order n, as a relation: the linear
    # coefficient is one, so L0 stays the small signal inductance
    taylor(n) = [k % 2 == 1 ? (-1.0)^((k-1)÷2)/factorial(k) : 0.0 for k in 1:n]
    jpa(jj) = Circuit([("p1","1","0",Port(1)), ("c1","1","2",Capacitor(100e-15)),
        ("lj","2","0",jj), ("c2","2","0",Capacitor(1e-12))])

    @testset "the relation, its derivative and its second derivative" begin
        p = PolynomialCPR(taylor(5))
        d = JC.cprderivative(p)
        d2 = JC.cprderivative(d)
        for phi in (0.0, 0.3, -0.8)
            @test p(phi) ≈ phi - phi^3/6 + phi^5/120
            @test d(phi) ≈ 1 - phi^2/2 + phi^4/24
            @test d2(phi) ≈ -phi + phi^3/6
            # the expansion of the Josephson relation, to its own order
            @test isapprox(p(phi), sin(phi), atol = abs(phi)^7)
            @test isapprox(d(phi), cos(phi), atol = abs(phi)^6)
            @test isapprox(d2(phi), -sin(phi), atol = abs(phi)^5)
        end
        # differentiating past the degree gives the zero polynomial rather
        # than an error, so a caller need not know the degree
        zero5 = foldl((f, _) -> JC.cprderivative(f), 1:6; init = d2)
        @test zero5(0.7) == 0
        # the linear coefficient is the normalization which makes `L0` the
        # small signal inductance
        @test_throws ArgumentError PolynomialCPR([2.0, 0.0, -1/6])
        @test_throws ArgumentError PolynomialCPR(Float64[])
    end

    @testset "a polynomial junction compiles as a junction" begin
        c = jpa(NonlinearInductor(L0, PolynomialCPR(taylor(5))))
        psc = JC.compile(c)
        @test psc.componenttypes == [:P, :R, :C, :Lj, :C]
        # it is a junction like any other, and its relation rides beside
        # the table keyed by the same index
        @test psc.junctions == [4]
        @test haskey(psc.junctioncprs, 4)
        @test psc.junctioncprs[4].a == [0.0; taylor(5)]
        # a sinusoidal one records nothing
        @test isempty(JC.compile(jpa(JosephsonJunction(L0))).junctioncprs)
        # a relation the solvers cannot write down is refused
        @test_throws JC.ComponentNotSupportedError JC.compile(
            jpa(NonlinearInductor(L0, tanh, x -> sech(x)^2)))
        # the transient carries the same relation, in the order of its
        # junction rows
        tp = transientproblem(jpa(NonlinearInductor(L0, PolynomialCPR(taylor(5)))))
        @test tp.relations.value[1, :] == [0.0; taylor(5)]
        @test tp.relations.sinusoidal == [false]
        @test isnothing(transientproblem(jpa(JosephsonJunction(L0))).relations)
    end

    @testset "the expansion of the Josephson relation approaches it" begin
        w, ip = 2*pi*4.75e9, 4*0.00565e-6
        src = [(mode = (1,), port = 1, current = ip)]
        ref = hbnlsolve((w,), (10,), src, jpa(JosephsonJunction(L0)),
            Dict{Symbol,Float64}(); keyedarrays = false)
        # the junction is driven to two thirds of a radian, where the
        # expansion is not trivially exact
        @test 0.6 < 2*maximum(abs, ref.nodeflux) < 0.7
        errors = map((3, 5, 7)) do n
            s = hbnlsolve((w,), (10,), src,
                jpa(NonlinearInductor(L0, PolynomialCPR(taylor(n)))),
                Dict{Symbol,Float64}(); keyedarrays = false)
            maximum(abs, s.nodeflux .- ref.nodeflux)/maximum(abs, ref.nodeflux)
        end
        # each pair of terms buys about two orders of magnitude here
        @test errors[1] < 1e-2
        @test errors[2] < errors[1]/50
        @test errors[3] < errors[2]/50
    end

    @testset "an array of junctions is one element with its relation" begin
        # N junctions in series carry one current and divide the phase, so
        # the array is one element of relation `N*sin(phi/N)` and small
        # signal inductance `N*Lj`. Written as a polynomial, that element
        # replaces the array.
        N, fp, fs, ip = 3, 4.75e9, 4.74e9, 4*0.00565e-6
        Lj = L0/N
        arrayrelation(n) = [k % 2 == 1 ?
            (-1.0)^((k-1)÷2)/(factorial(k)*N^(k-1)) : 0.0 for k in 1:n]
        explicit = Circuit(vcat(
            Any[("p1","1","0",Port(1)), ("c1","1","2",Capacitor(100e-15)),
                ("c2","2","0",Capacitor(1e-12))],
            [("jj$k", k == 1 ? "2" : "n$(k-1)", k == N ? "0" : "n$k",
                JosephsonJunction(Lj)) for k in 1:N]))
        run(c) = hbsolve([2*pi*fs], (2*pi*fp,),
            [(mode = (1,), port = 1, current = ip)], (8,), (16,), c; ftol = 1e-14)
        a = run(explicit)
        gain = abs2(a.linearized.S((0,),1,(0,),1,1))
        @test gain > 2                      # the array is driven into gain
        for (n, rtol) in ((5, 1e-4), (9, 1e-7))
            b = run(jpa(NonlinearInductor(L0, PolynomialCPR(arrayrelation(n)))))
            @test abs2(b.linearized.S((0,),1,(0,),1,1)) ≈ gain rtol=rtol
            @test b.linearized.QE((0,),1,(0,),1,1) ≈
                a.linearized.QE((0,),1,(0,),1,1) rtol=rtol
            # the two circuits do not have the same nodes, so what is
            # compared is what a user measures: the port quantities
        end
    end

    @testset "a quadratic term mixes three waves" begin
        # what a SNAIL is for: biased away from its symmetric point its
        # relation has a quadratic term, and the amplifier is pumped at
        # twice the signal. The Josephson relation is odd and does not mix
        # three waves at all.
        fs, ip = 4.74e9, 0.1e-6
        run(jj) = hbsolve([2*pi*fs], (2*pi*2*4.75e9,),
            [(mode = (1,), port = 1, current = ip)], (8,), (16,), jpa(jj);
            ftol = 1e-14, threewavemixing = true)
        snail = run(NonlinearInductor(L0, PolynomialCPR([1.0, 0.3, -1/6])))
        junction = run(JosephsonJunction(L0))
        @test abs2(snail.linearized.S((0,),1,(0,),1,1)) > 3
        # the odd Josephson relation has no quadratic term, so the pump
        # at twice the signal does not couple the signal to its idler and
        # the circuit is a passive reflection
        @test abs2(junction.linearized.S((0,),1,(0,),1,1)) ≈ 1 rtol=1e-6
        @test 0 < snail.linearized.QE((0,),1,(0,),1,1) < 1
    end

    @testset "a polynomial junction beside a sinusoidal one" begin
        # a circuit holding both takes the mixed path, which writes the
        # sinusoidal columns over the polynomial ones: to a high enough
        # order the two are the same circuit
        w, ip = 2*pi*4.75e9, 4*0.00565e-6
        src = [(mode = (1,), port = 1, current = ip)]
        pair(second) = Circuit([("p1","1","0",Port(1)),
            ("c1","1","2",Capacitor(100e-15)), ("jj","2","3",JosephsonJunction(L0)),
            ("lj","3","0",second), ("c2","2","0",Capacitor(1e-12))])
        both = hbnlsolve((w,), (10,), src, pair(JosephsonJunction(L0)),
            Dict{Symbol,Float64}(); keyedarrays = false)
        mixed = hbnlsolve((w,), (10,), src,
            pair(NonlinearInductor(L0, PolynomialCPR(taylor(11)))),
            Dict{Symbol,Float64}(); keyedarrays = false)
        @test maximum(abs, mixed.nodeflux .- both.nodeflux) <
            1e-9*maximum(abs, both.nodeflux)
    end

    @testset "the second and third derivatives of the problem" begin
        # the problem interface, which an external solver drives, takes
        # directional derivatives to third order: for the Josephson
        # relation those are `-sin`, `cos`, `-sin` and `-cos`, and for a
        # polynomial the differentiated polynomials. Each is checked
        # against a central difference of the one below it.
        prob = JC.hbnonlinearproblem((2*pi*4.75001e9,), (8,),
            [(mode = (1,), port = 1, current = 0.02e-6)],
            jpa(NonlinearInductor(L0, PolynomialCPR([1.0, 0.25, -1/6]))),
            Dict{Symbol,Float64}(); Nevaluationharmonics = (16,))
        n = length(prob)
        u = 0.05 .* randn(n)
        v, w = randn(n), randn(n)
        z = randn(n)./sqrt(n); zh = z./norm(z)
        h = 1e-6*max(norm(u), 1)
        relerr(a, b) = maximum(abs, a .- b)/max(maximum(abs, b), 1e-30)
        # symmetric in its directions
        @test relerr(JC.hbd3F!(zeros(n), prob, u, v, w, zh),
            JC.hbd3F!(zeros(n), prob, u, zh, v, w)) < 1e-12
        # the second derivative differentiates the Jacobian, the third the
        # second
        fd2 = (JC.hbjvp!(zeros(n), prob, u .+ h.*zh, v) .-
               JC.hbjvp!(zeros(n), prob, u .- h.*zh, v))./(2h)
        @test relerr(JC.hbd2F!(zeros(n), prob, u, v, zh), fd2) < 1e-5
        fd3 = (JC.hbd2F!(zeros(n), prob, u .+ h.*zh, v, w) .-
               JC.hbd2F!(zeros(n), prob, u .- h.*zh, v, w))./(2h)
        @test relerr(JC.hbd3F!(zeros(n), prob, u, v, w, zh), fd3) < 1e-4
    end

    @testset "the sensitivity of a polynomial junction's amplifier" begin
        # the reverse contraction of the linearized outputs reads the
        # second derivative of the relation at the pump: the sensitivity
        # of the gain to a component value against a central difference of
        # the gain itself
        fp, fs, ip = 4.75e9, 4.74e9, 4*0.00565e-6
        Cc = 100e-15
        run(c; kw...) = hbsolve([2*pi*fs], (2*pi*fp,),
            [(mode = (1,), port = 1, current = ip)], (8,), (16,),
            Circuit([("p1","1","0",Port(1)), ("c1","1","2",Capacitor(c)),
                ("lj","2","0",NonlinearInductor(L0,
                    PolynomialCPR([1.0, 0.0, -1/6, 0.0, 1/120]))),
                ("c2","2","0",Capacitor(1e-12))]);
            ftol = 1e-14, keyedarrays = false, kw...)
        s = run(Cc; returnSsensitivity = true, sensitivitynames = ["c1"])
        @test size(s.linearized.Ssensitivity, 3) == 1
        dS = s.linearized.Ssensitivity[1, 1, 1]
        # the sensitivity is with respect to the logarithm of the value
        d = 1e-4
        plus, minus = run(Cc*exp(d)), run(Cc*exp(-d))
        @test isapprox(dS, (plus.linearized.S[1,1] - minus.linearized.S[1,1])/(2d);
            rtol = 1e-3)
    end

    @testset "the transient steps the same relation" begin
        # the time domain solver against harmonic balance on the same
        # circuit: the pumped amplitude of the fundamental after the pump
        # has settled. A quadratic term rectifies, so the comparison keeps
        # the direct current and the even harmonics which it generates.
        fp, ip = 4.75e9, 4*0.00565e-6
        ramp(t) = t <= 0 ? 0.0 : t >= 2e-9 ? 1.0 : (1 - cospi(t/2e-9))/2
        pump(t) = 2*ip*ramp(t)*cospi(2*fp*t)
        window(t) = 40e-9 <= t <= 60e-9 ? sinpi((t - 40e-9)/20e-9)^2 : 0.0
        for jj in (JosephsonJunction(L0),
                   NonlinearInductor(L0, PolynomialCPR(taylor(9))),
                   NonlinearInductor(L0, PolynomialCPR([1.0, 0.25, -1/6])))
            c = jpa(jj)
            hb = hbnlsolve((2*pi*fp,), (10,), [(mode = (1,), port = 1, current = ip)],
                c, Dict{Symbol,Float64}(); keyedarrays = false, dc = true,
                even = true, odd = true)
            k = findfirst(==((1,)), hb.frequencies.modes)
            expected = 2im*2*pi*fp*JC.phi0*hb.nodeflux[k]
            sol = transientsolve(transientproblem(c;
                sources = [TransientSource(1, pump)]), (0.0, 60e-9);
                dt = 4e-12, method = GaussLegendre())
            got = transientdemodulate(sol, 1, fp; quantity = :voltage, window)
            @test got ≈ expected rtol=1e-4
        end
    end

    @testset "the transient's tangent, adjoint and batch" begin
        # the derivative of the relation drives the tangent and the
        # adjoint, and a batch steps every condition on one factorization:
        # the tangent against a central difference of the solve, the
        # adjoint against the tangent, and the batch against its members
        L, T, dt = 1e-9, 0.5e-9, 2e-12
        cpr = PolynomialCPR([1.0, 0.25, -1/6])
        drive(t) = t <= 0 ? 0.0 : 0.3e-6*sinpi(2*3e9*t)
        probe(t, p) = 1e-8*sinpi(2*1.1e9*t + p)
        prob = transientproblem(jpa(NonlinearInductor(L, cpr));
            sources = [TransientSource(1, drive)])
        rec = transientsolve(prob, (0.0, T); dt, record = :phases,
            method = GaussLegendre(), rtol = 1e-12)
        currents = [probe(t, 1.0) for _ in 1:1, t in rec.times]
        tg = transienttangent(rec, currents)
        eps = 1e-3
        shifted(s) = transientproblem(prob; sources = [TransientSource(1,
            let s = s; t -> drive(t) + s*eps*probe(t, 1.0); end)])
        sp = transientsolve(shifted(1), (0.0, T); dt, method = GaussLegendre(), rtol = 1e-12)
        sm = transientsolve(shifted(-1), (0.0, T); dt, method = GaussLegendre(), rtol = 1e-12)
        @test tg.outgoing ≈ (sp.outgoing .- sm.outgoing)./(2*eps) rtol=1e-5
        # the adjoint contracts the same responses as the tangent
        weights = [cospi(2*1.7e9*t) for _ in 1:1, t in rec.times]
        ad = transientadjoint(rec, weights)
        @test sum(weights .* tg.outgoing) ≈ sum(ad.currents .* currents) rtol=1e-8
        # a batch of two conditions equals its members
        half = transientproblem(prob; sources = [TransientSource(1, t -> drive(t)/2)])
        b = transientsolve([prob, half], (0.0, T); dt, record = :phases,
            method = GaussLegendre(), rtol = 1e-12)
        @test b[1].outgoing ≈ rec.outgoing rtol=1e-10
        @test transienttangent(b, currents).outgoing[:, :, 1] ≈ tg.outgoing rtol=1e-9
    end

    @testset "the transient's projected constraint and noise" begin
        # a junction alone across an unterminated port is an algebraic
        # constraint, which the Gauss-Legendre rule projects onto: the
        # projection evaluates the relation of its own junctions
        L, w, a = 1e-9, 2*pi*1e9, 0.2
        p = transientproblem(Circuit([("p","1","0",Port(1; termination = nothing)),
            ("lj","1","0",NonlinearInductor(L, PolynomialCPR([1.0, 0.0, -1/6])))]);
            sources = [TransientSource(1, t -> JC.phi0/L*a*(1 - cos(w*t)))])
        @test p.inertialess == [[1]] && p.algebraic == [[1]]
        # the flux solves the constraint `f(phi) = L*I(t)/phi0` at every
        # time, which for this relation is a cubic in the phase
        s = transientsolve(p, (0.0, 1e-9); dt = 1e-9/80, method = GaussLegendre(),
            record = :states, rtol = 1e-12, atol = 1e-13)
        phi = vec(s.flux)
        drive = [a*(1 - cos(w*t)) for t in s.times]
        @test maximum(abs, (phi .- phi.^3 ./ 6) .- drive) < 1e-9
        # the noise linearizes about the same relation: its derivative is
        # the differential inductance the baths see
        rest = transientproblem(jpa(
            NonlinearInductor(L, PolynomialCPR([1.0, 0.0, -1/6, 0.0, 1/120]))))
        nsol = transientsolve(rest, (0.0, 1e-9 - 2e-12); dt = 2e-12,
            record = :phases, method = GaussLegendre())
        plan = transientquantumplan(nsol.times, [4/1e-9])
        noise = transientnoise(nsol, plan; frequencies = [3/1e-9, 4/1e-9],
            weights = [1e9, 1e9], inputs = plan)
        @test noise.diagnostics.passed
        # at rest the amplifier is a passive reflection with vacuum noise
        @test noise.covariance ≈ [0.5 0.0; 0.0 0.5] atol=1e-6
    end

    @testset "the residual, the Jacobian and the Hessian" begin
        # the three evaluations of the relation, against finite differences
        # of each other: the Jacobian is the derivative of the residual and
        # the Hessian the derivative of the Jacobian, which is what pins
        # the second derivative of the polynomial and its sign
        w = 2*pi*4.75e9
        src = [(mode = (1,), port = 1, current = 4*0.00565e-6)]
        d = hbnlsolve((w,), (8,), src,
            jpa(NonlinearInductor(L0, PolynomialCPR([1.0, 0.25, -1/6]))),
            Dict{Symbol,Float64}(); debugJacobian = true)
        sys = d.sys
        nr = length(d.xr)
        xr = 0.4*randn(nr)
        vr, wr = randn(nr), randn(nr)
        h = 1e-7
        Jvr, Fp, Fm = zeros(nr), zeros(nr), zeros(nr)
        JC.setpoint!(sys, xr)
        JC.jacobianvectorproduct!(Jvr, sys, vr)
        JC.setpoint!(sys, xr .+ h.*vr); JC.residual!(Fp, sys)
        JC.setpoint!(sys, xr .- h.*vr); JC.residual!(Fm, sys)
        @test isapprox(Jvr, (Fp .- Fm)./(2*h), atol = 1e-4,
            norm = v -> maximum(abs, v))
        # the assembled Jacobian is the same operator
        JC.setpoint!(sys, xr)
        JC.jacobian!(d.Jr, sys)
        @test isapprox(Jvr, d.Jr*vr, atol = 1e-11, norm = v -> maximum(abs, v))
        # the Hessian, symmetric in its directions and the derivative of
        # the Jacobian
        Hvw, Hwv, Jp, Jm = zeros(nr), zeros(nr), zeros(nr), zeros(nr)
        JC.setpoint!(sys, xr)
        JC.hessianvectorproduct!(Hvw, sys, vr, wr)
        JC.hessianvectorproduct!(Hwv, sys, wr, vr)
        @test isapprox(Hvw, Hwv, atol = 1e-12, norm = v -> maximum(abs, v))
        JC.setpoint!(sys, xr .+ h.*wr); JC.jacobianvectorproduct!(Jp, sys, vr)
        JC.setpoint!(sys, xr .- h.*wr); JC.jacobianvectorproduct!(Jm, sys, vr)
        @test isapprox(Hvw, (Jp .- Jm)./(2*h), atol = 1e-4,
            norm = v -> maximum(abs, v))
    end
end
