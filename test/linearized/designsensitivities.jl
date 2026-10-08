using JosephsonCircuits
using LinearAlgebra
using SparseArrays
using Test

# the bytes `f()` allocates, measured in a function: `@allocated` compiles
# the code around it, which in a testset is the whole testset
bytesallocated(f) = @allocated f()

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
        @test bytesallocated(je) < 2*bytesallocated(js)
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
        # a value which does not depend on the parameter collapses to zero,
        # a frequency dependent leaf among them
        @test CV.derivative(CV.tocv(b) + 1, :a) == CV.Constant(0)
        @test d(FrequencyDependent(w -> 1.0), :a, defs) == 0
        # an undefined parameter in the derivative is refused, as is a
        # direction which depends on the frequency
        @test_throws ArgumentError d(a*b, :a, Dict(a => 2.0))
        @test_throws ArgumentError d(a*FrequencyDependent(w -> 1.0), :a, defs)
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

    @testset "a parameter defined in terms of another" begin
        # one capacitor written in a name defined by an expression in the
        # design parameter, one in the parameter itself: the derivative
        # follows the definition by the chain rule, 2 and 1, whichever
        # spelling names the parameters, and the sensitivity of the whole
        # solve is the finite difference of the definitions
        JosephsonCircuits.@params C0 Cd
        two(c1, c2) = Circuit([(:p, 1, 0, Port(1)), (:c1, 1, 0, Capacitor(c1)),
            (:c2, 1, 0, Capacitor(c2))])
        spellings = (
            (two(:Cd, :C0), Dict(:Cd => 2*C0, :C0 => 1e-12)),
            (two("Cd", "C0"), Dict("Cd" => 2*C0, "C0" => 1e-12)),
            (two(Cd, C0), Dict(Cd => 2*C0, C0 => 1e-12)))
        wl, wpl = [2*pi*2e9], (2*pi*1e9,)
        srcl = [(mode = (1,), port = 1, current = 0.0)]
        opts = (; returnQE = false, returnCM = false)
        Sat(c, defs, x) = Array(hbsolve(wl, wpl, srcl, (1,), (1,), c,
            merge(Dict{Any,Any}(defs), Dict{Any,Any}(:C0 => x)); opts...).linearized.S)
        c, defs = first(spellings)
        h = 1e-6*1e-12
        fd = vec((Sat(c, defs, 1e-12 + h) .- Sat(c, defs, 1e-12 - h))./(2*h))
        for (c, defs) in spellings
            names, _, J = JosephsonCircuits.designjacobian(c, defs)
            @test names == ["p/termination", "c1", "c2"] && J == ComplexF64[0; 2; 1;;]
            # the derived name is no default parameter of its own
            r = designsensitivities(c, defs, wl, wpl, srcl, (1,), (1,); opts...)
            @test collect(JosephsonCircuits.AxisKeys.axiskeys(r.dSdp, :parameter)) == [:C0]
            @test norm(vec(Array(r.dSdp)) .- fd)/norm(fd) < 1e-6
        end
    end

    @testset "a frequency dependent value beside a design parameter" begin
        # a lossy dielectric written as a closure of the frequency, which no
        # parameter moves, beside the coupling capacitance which is the
        # design parameter: it is no row of the Jacobian, and the solve's
        # sensitivity is the finite difference; a frequency dependent value
        # which the parameter moves is refused
        JosephsonCircuits.@params Cc
        lossy = FrequencyDependent(w -> 1000e-15*(1 - 1e-4im*sign(w)))
        jpa(cj) = Circuit([(:p1, 1, 0, Port(1)), (:cc, 1, 2, Capacitor(Cc)),
            (:jj, 2, 0, JosephsonJunction(1000e-12)), (:cj, 2, 0, Capacitor(cj))])
        defs = Dict(Cc => 100e-15)
        @test JosephsonCircuits.designjacobian(jpa(lossy), defs)[1] ==
            ["p1/termination", "cc", "jj"]
        r = designsensitivities(jpa(lossy), defs, ws, wp, src, (2,), (8,))
        fd = fdS(jpa(lossy), defs, Cc)
        @test norm(vec(Array(r.dSdp)) .- fd)/norm(fd) < 1e-4
        @test_throws ArgumentError JosephsonCircuits.designjacobian(
            jpa(:Cj), Dict(Cc => 100e-15, :Cj => Cc*lossy/(100e-15)))
    end

    @testset "plain arrays" begin
        # with keyedarrays = false the derivative is the solve's own plain
        # scattering sensitivity, whose component axis is the parameters
        c = Circuit([(:p1, 1, 0, Port(1)), (:c1, 1, 0, Capacitor(:C1)),
            (:c2, 1, 0, Capacitor(:C2))])
        defs = Dict(:C1 => 1e-12, :C2 => 2e-12)
        args = (2*pi*[3e9, 4e9], (2*pi*5e9,),
            [(mode = (1,), port = 1, current = 0.0)], (0,), (1,))
        keyed = designsensitivities(c, defs, args...)
        plain = designsensitivities(c, defs, args...; keyedarrays = false)
        @test plain.dSdp isa Array{ComplexF64,4} && size(plain.dSdp) == (1, 1, 2, 2)
        @test vec(plain.dSdp) == vec(Array(keyed.dSdp))
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

    @testset "direct current, both mixing orders and negative sidebands" begin
        # the chain rule has to carry every retained mode of a mixed grid:
        # a tone with both mixing orders and direct current, whose signal
        # modes include ones at negative physical frequency. Forward and
        # reverse must agree, and both must agree with differences of the
        # whole solve.
        JosephsonCircuits.@params scale loss
        c = Circuit([(:p, 1, 0, Port(1)),
            (:cc, 1, 2, Capacitor(100e-15/(1 + im*loss))),
            (:jj, 2, 0, JosephsonJunction(1e-9/scale)),
            (:cj, 2, 0, Capacitor(1e-12*scale))])
        defs = Dict(scale => 1.0, loss => 0.02)
        pumps = (2*pi*4.75001e9,)
        drive = [(mode = (1,), port = 1, current = 0.00565e-6),
            (mode = (0,), port = 1, current = 1e-9)]
        signals = 2*pi*[4.41e9, 4.59e9]
        opts = (; dc = true, threewavemixing = true, fourwavemixing = true,
            atol = 1e-12)
        forward = designsensitivities(c, defs, signals, pumps, drive,
            (2,), (8,); sensitivitymode = :forward, opts...)
        reverse = designsensitivities(c, defs, signals, pumps, drive,
            (2,), (8,); sensitivitymode = :reverse, opts...)
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
            up = hbsolve(signals, pumps, drive, (2,), (8,), c,
                at(defs[parameter] + h); opts...)
            dn = hbsolve(signals, pumps, drive, (2,), (8,), c,
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

    @testset "the residual columns of a shared parameter merge sparse" begin
        # the pump of a chain solved once, then the sweep with the pairs of
        # parameters which each move the junctions, or the capacitors, of
        # two cells: the residual derivative columns of a parameter's pairs
        # merge into one and stay sparse, so the bytes of the sweep grow in
        # proportion to the chain
        function sweep(n)
            c = compile(Circuit(vcat(Any[("P1", "1", "0", Port(1))],
                Any[x for i in 1:n for x in (
                    ("Lj$i", "$i", "$(i + 1)", JosephsonJunction(100e-12)),
                    ("C$i", "$i", "0", Capacitor(40e-15)))],
                Any[("R2", "$(n + 1)", "0", Resistor(50.0))])))
            K = cld(n, 2)
            pairs = vcat([("Lj$i", cld(i, 2), 1/100e-12 + 0im) for i in 1:n],
                [("C$i", K + cld(i, 2), 1/40e-15 + 0im) for i in 1:n])
            nl = hbnlsolve((2*pi*7e9,), (8,), [(mode = (1,), port = 1,
                current = 1e-8)], c; returnoperatingpoint = true,
                keyedarrays = false)
            op = nl.operatingpoint
            dF = JosephsonCircuits.calcresidualsensitivity(op, c,
                numericmatrices(c, Dict(); Nmodes = op.Nmodes),
                [JosephsonCircuits.componentindex(c, p[1]) for p in pairs],
                [p[3] for p in pairs])
            return () -> hblinsolve([2*pi*5e9], c; nonlinear = nl,
                nbatches = 1, keyedarrays = false, sensitivitypairs = pairs,
                nsensitivityparameters = 2K, sensitivityresidual = dF,
                returnSsensitivity = true)
        end
        function bytes(n)
            run = sweep(n)
            run()
            return @allocated run()
        end
        @test bytes(512) < 6*bytes(128)
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

    @testset "the reference impedance of a port which owns no termination" begin
        # both ports normalize their waves to the design parameter R and own
        # no termination: the first is loaded by a resistor of R, the second
        # by a fixed one, so R moves the first port's load and both ports'
        # normalizations, which the sensitivities carry as stamps with no
        # entries, one per port
        circuit = Circuit(Any[
            ("P1", "1", "0", Port(1; Z0 = :R, termination = nothing)),
            ("R1", "1", "0", Resistor(:R)), ("C1", "1", "2", Capacitor(100e-15)),
            ("Lj1", "2", "0", JosephsonJunction(1000e-12)),
            ("C2", "2", "0", Capacitor(1000e-15)), ("C3", "2", "3", Capacitor(10e-15)),
            ("P2", "3", "0", Port(2; Z0 = :R, termination = nothing)),
            ("R2", "3", "0", Resistor(50.0))])
        defs = Dict(:R => 50.0)
        r = designsensitivities(circuit, defs, ws, wp, src, (2,), (8,))
        fd = fdS(circuit, defs, :R)
        @test norm(vec(Array(r.dSdp)) .- fd)/norm(fd) < 1e-4
    end

    @testset "scattering block design parameters" begin
        # a theta-dependent series-impedance two port feeding a JPA: the
        # block parameter and a lumped parameter share the output axis,
        # and the block moves the pump operating point. The block depends
        # on theta through the derivative it states.
        Z0 = 50.0
        # the block and its derivative, one closure type for every block
        seriesS(theta) = w -> (z = 1/(im*w*theta*Z0);
            [z/(z+2) 2/(z+2); 2/(z+2) z/(z+2)])
        dS(theta) = w -> (z = 1/(im*w*theta*Z0); d = -2*z/(theta*(z+2)^2);
            [d -d; -d d])
        function makeblk(theta; analytic = true)
            blk = analytic ?
                ScatteringParameters(seriesS(theta); nports = 2,
                    grounded = false, derivatives = (theta = dS(theta),)) :
                ScatteringParameters(seriesS(theta); nports = 2,
                    grounded = false)
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
        r = designsensitivities(circuit, defs, ws, wp, src, (2,), (8,);
            sensitivitymode = :forward)
        @test collect(JosephsonCircuits.AxisKeys.axiskeys(r.dSdp, :parameter)) == [:Lj, :theta]
        fdL = fdS(circuit, defs, :Lj)
        @test norm(vec(Array(r.dSdp(parameter = :Lj))) .- fdL)/norm(fdL) < 1e-4
        h = 1e-6*theta
        Sp = hbsolve(ws, wp, src, (2,), (8,), makeblk(theta + h), defs).linearized.S
        Sm = hbsolve(ws, wp, src, (2,), (8,), makeblk(theta - h), defs).linearized.S
        fdt = vec((Array(Sp) .- Array(Sm))./(2h))
        @test norm(vec(Array(r.dSdp(parameter = :theta))) .- fdt)/norm(fdt) < 1e-4

        # three instances of one definition along a line share the
        # parameter, whose pairs make one stamp
        function makeline(theta)
            blk = ScatteringParameters(seriesS(theta); nports = 2,
                derivatives = (theta = dS(theta),))
            comps = Any[(:p1, 1, 0, Port(1; Z0 = Z0))]
            for k in 1:3
                push!(comps, (Symbol(:b, k), k, k + 1, blk))
                push!(comps, (Symbol(:c, k), k + 1, 0, Capacitor(10e-15)))
            end
            push!(comps, (:cc, 4, 5, Capacitor(100e-15)),
                (:jj, 5, 0, JosephsonJunction(:Lj)),
                (:cj, 5, 0, Capacitor(1000e-15)))
            return Circuit(comps)
        end
        rl = designsensitivities(makeline(theta), defs, ws, wp, src, (2,),
            (8,); parameters = (:theta,), sensitivitymode = :forward)
        Sp = hbsolve(ws, wp, src, (2,), (8,), makeline(theta + h), defs).linearized.S
        Sm = hbsolve(ws, wp, src, (2,), (8,), makeline(theta - h), defs).linearized.S
        fdl = vec((Array(Sp) .- Array(Sm))./(2h))
        @test norm(vec(Array(rl.dSdp(parameter = :theta))) .- fdl)/norm(fdl) < 1e-4

        # a block which states no derivative depends on no parameter and
        # costs nothing
        plain = makeblk(theta; analytic = false)
        @test isempty(JosephsonCircuits.designblockjacobian(plain, [:theta, :Lj]))
        rp = designsensitivities(plain, defs, ws, wp, src, (2,), (8,))
        @test collect(JosephsonCircuits.AxisKeys.axiskeys(rp.dSdp, :parameter)) == [:Lj]
        @test all(isfinite, Array(rp.dSdp))

        # the reverse order takes the residual column of a block parameter
        # as it takes a component's, the three instances' merged into one,
        # and agrees with the forward order
        for (c, q, fw) in ((circuit, nothing, r), (makeline(theta), (:theta,), rl))
            rv = designsensitivities(c, defs, ws, wp, src, (2,), (8,);
                parameters = q, sensitivitymode = :reverse)
            @test isapprox(Array(rv.dSdp), Array(fw.dSdp); rtol = 1e-10)
        end
    end

    @testset "the block pairs of many blocks and parameters" begin
        # a line of instances of one block, which states derivatives for
        # two parameters, among as many parameters as blocks: each block
        # pairs with its two, in the order of the parameters, and its
        # derivatives are read once, so the bytes grow in proportion to
        # the line
        seriesS(theta) = w -> (z = 1/(im*w*theta*50.0);
            [z/(z+2) 2/(z+2); 2/(z+2) z/(z+2)])
        blk = ScatteringParameters(seriesS(1e-12); nports = 2,
            derivatives = (thetb = seriesS(2e-12), theta = seriesS(3e-12)))
        line(n) = compile(Circuit(vcat(Any[(:p1, 1, 0, Port(1))],
            Any[(Symbol(:b, k), k, k + 1, blk) for k in 1:n],
            Any[(:r2, n + 1, 0, Resistor(50.0))])))
        names(n) = [:theta; [Symbol(:c, k) for k in 2:(n - 1)]; :thetb]
        pairs = JosephsonCircuits.designblockjacobian(line(3), names(3))
        @test [(p[1], p[2]) for p in pairs] ==
            [("b$k", j) for k in 1:3 for j in (1, 3)]
        function bytes(n)
            c, q = line(n), names(n)
            run() = JosephsonCircuits.designblockjacobian(c, q)
            run()
            return @allocated run()
        end
        @test bytes(128) < 6*bytes(32)
    end

    @testset "many block parameters take the reverse order" begin
        # a fifty ohm line of instances of one series inductance block
        # ending in a pumped junction resonator, each block its own
        # parameter through a block pair: past eight parameters per output
        # port and mode pair the default order is the reverse one, whose
        # cost does not grow with the parameters, so the bytes grow in
        # proportion to the line
        seriesL(L) = w -> (z = im*w*L/50.0;
            [z/(z+2) 2/(z+2); 2/(z+2) z/(z+2)])
        dseriesL(L) = w -> (z = im*w*L/50.0; d = 2im*w/(50.0*(z+2)^2);
            [d -d; -d d])
        blk = ScatteringParameters(seriesL(1e-9); nports = 2)
        dblk = ScatteringParameters(dseriesL(1e-9); nports = 2)
        function bytes(n)
            c = compile(Circuit(vcat(Any[(:p1, 1, 0, Port(1))],
                Any[x for k in 1:n for x in ((Symbol(:b, k), k, k + 1, blk),
                    (Symbol(:c, k), k + 1, 0, Capacitor(400e-15)))],
                Any[(:cc, n + 1, n + 2, Capacitor(100e-15)),
                    (:jj, n + 2, 0, JosephsonJunction(1000e-12)),
                    (:cj, n + 2, 0, Capacitor(1000e-15))])))
            pairs = [("b$k", k, dblk) for k in 1:n]
            run() = hbsolve(ws[1:1], wp, src, (2,), (8,), c;
                sensitivityblockpairs = pairs, nsensitivityparameters = n,
                returnSsensitivity = true, keyedarrays = false, nbatches = 1)
            run()
            return @allocated run()
        end
        @test bytes(512) < 6*bytes(128)
    end

    @testset "a pumped block among the components" begin
        # a short written as a block of the linearized kind which does not
        # convert and as an ordinary block, behind a coupling capacitor
        # which is the design parameter: the block states no derivative in
        # either form, and the two give one sensitivity
        wpb = (2*pi*4.75e9,)
        wsb = 2*pi*[4.4e9, 4.6e9]
        withblock(b) = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)),
            (:cc, 1, 2, Capacitor(:Cc)), (:b, 2, b), (:c2, 2, 0, Capacitor(500e-15))])
        pumped = LinearizedScattering([fill(-1.0 + 0im, 1, 1)], wpb[1];
            harmonics = [0], nports = 1)
        ordinary = ScatteringParameters(fill(-1.0 + 0im, 1, 1); zref = 50.0,
            noise = Lossless())
        defs = Dict(:Cc => 100e-15)
        dS(b) = Array(designsensitivities(withblock(b), defs, wsb, wpb, [],
            (2,), (4,); parameters = (:Cc,)).dSdp)
        @test isapprox(dS(pumped), dS(ordinary); rtol = 1e-10)
        # a block pair is of an ordinary block, with an ordinary derivative:
        # a pumped block states none
        for (b, d) in ((pumped, ordinary), (pumped, pumped), (ordinary, pumped))
            @test_throws ArgumentError hbsolve(wsb, wpb, [], (2,), (4,),
                withblock(b), defs; sensitivityblockpairs = [("b", 1, d)],
                nsensitivityparameters = 1)
        end
    end

    @testset "a block design parameter with direct current" begin
        # a series resistance block in the path of a direct current, which
        # it reads from its own zero frequency data: its parameter moves the
        # direct current rows of the pump solve as well as the harmonic ones
        Z0 = 50.0
        function makeres(R)
            seriesS(w) = (z = R/Z0; [z/(z+2) 2/(z+2); 2/(z+2) z/(z+2)])
            dS(w) = (z = R/Z0; d = 2/(Z0*(z+2)^2); [d -d; -d d])
            blk = ScatteringParameters(seriesS; nports = 2, grounded = false,
                derivatives = (R = dS,))
            return Circuit(
                [:p1 => Port(1; termination = nothing), :r1 => Resistor(50.0),
                 :cc => blk, :jj => JosephsonJunction(:Lj),
                 :c2 => Capacitor(1000.0e-15)],
                [((:p1, 1), (:r1, 1), (:cc, 1, 1)),
                 ((:cc, 2, 1), (:jj, 1), (:c2, 1)),
                 ((:cc, 1, 2), (:cc, 2, 2), (:jj, 2), (:c2, 2), (:r1, 2),
                  (:p1, 2), Ground)])
        end
        R = 5.0
        srcdc = [(mode = (1,), port = 1, current = 0.00565e-6),
            (mode = (0,), port = 1, current = 1e-8)]
        opts = (; dc = true, threewavemixing = true, fourwavemixing = true)
        defs = Dict(:Lj => 1000.0e-12)
        r = designsensitivities(makeres(R), defs, ws, wp, srcdc, (2,), (8,);
            parameters = (:R,), sensitivitymode = :forward, opts...)
        h = 1e-6*R
        Sp = hbsolve(ws, wp, srcdc, (2,), (8,), makeres(R + h), defs; opts...).linearized.S
        Sm = hbsolve(ws, wp, srcdc, (2,), (8,), makeres(R - h), defs; opts...).linearized.S
        fd = vec((Array(Sp) .- Array(Sm))./(2h))
        @test norm(vec(Array(r.dSdp(parameter = :R))) .- fd)/norm(fd) < 1e-4
        # the reverse order, through the block's column in the canonical
        # coordinates, its zero frequency rows included, agrees
        rv = designsensitivities(makeres(R), defs, ws, wp, srcdc, (2,), (8,);
            parameters = (:R,), sensitivitymode = :reverse, opts...)
        @test isapprox(Array(rv.dSdp), Array(r.dSdp); rtol = 1e-10)
    end
end
