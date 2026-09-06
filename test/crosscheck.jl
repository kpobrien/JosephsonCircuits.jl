using JosephsonCircuits
using LinearAlgebra
using SpecialFunctions
using Test

isdefined(Main, :crosscheckcases) || include("testcircuits.jl")

# The systematic comparisons: every circuit class of `crosscheckcases`
# through every route to the same answer the package offers. A comparison
# is a function of a case; the table says which apply to which.

const CASES = crosscheckcases()

# the two tone form of a single pump source list: the pump modes gain a
# trailing zero, and the weak signal is the trailing tone
function twotonesources(case)
    np = length(case.wp)
    srcs = [(mode = (s.mode..., 0), port = s.port, current = s.current)
        for s in case.sources]
    push!(srcs, (mode = (ntuple(_ -> 0, np)..., 1), port = case.signalport,
        current = case.Is))
    return srcs
end

# the relative deviation of `a` from `b`, with `b`'s scale, and zero when
# both are negligible against `scale`
reldev(a, b, scale) = abs(a - b)/max(abs(b), scale)

@testset verbose = true "cross checks" begin

for case in CASES
@testset "$(case.name)" begin
    np = length(case.wp)
    nlkw = nonlinearkw(case.kw)
    solve(; kw...) = hbsolve(case.ws, case.wp, case.sources, case.Nmod,
        case.Npump, case.circuit, case.defs; case.kw..., ftol = 1e-12, kw...)

    if case.pumped
        # A. the weak second tone of the nonlinear solve against the
        # linearized solve: every signal and idler mode the two share, with
        # the conjugation rule for a mode the nonlinear solve keeps at the
        # negative of the linearized frequency. The nonlinear response
        # differs from the linearized one at first order in the ratio of
        # the weak tone to the drive, and the weak response is resolved to
        # the residual tolerance over that ratio, so the two are compared
        # to ten times the sum of both, on the scale of the largest entry
        @testset "weak tone against the linearized solve" begin
            ftol = 1e-14
            lin = solve(ftol = ftol)
            ratio = case.Is/maximum(s.current for s in case.sources)
            tol = 10*(ratio + ftol/ratio)
            modes = lin.linearized.modes
            ports = lin.linearized.portnumbers
            nzero = ntuple(_ -> 0, np)
            for (i, ws) in enumerate(case.ws)
                # (pump..., signal) order
                nl = hbnlsolve((case.wp..., ws), (case.Npump..., case.Nmod[1]),
                    twotonesources(case), case.circuit, case.defs;
                    nlkw..., ftol = ftol)
                @test nl.solverinfo.converged
                smax = maximum(abs(lin.linearized.S(m, p, n, q, i))
                    for m in modes, p in ports, n in modes, q in ports)
                compared = 0
                for m in modes, pout in ports
                    f = sum(m .* case.wp) + ws
                    key = f > 0 ? (m..., 1) : ((-1 .* m)..., -1)
                    key in nl.modes || continue
                    s1 = lin.linearized.S(m, pout, nzero, case.signalport, i)
                    s2 = nl.S(key, pout, (nzero..., 1), case.signalport)
                    f > 0 || (s2 = conj(s2))
                    @test abs(s2 - s1) <= tol*smax
                    compared += 1
                end
                @test compared >= length(ports)
            end
        end

        # B. the solvers agree: the scattering matrix at the pump
        # frequencies, the direct current voltages, the node fluxes of the
        # driven modes, and the zero frequency node fluxes modulo whole flux
        # quanta. The zero frequency flux of a node is a soft direction of
        # the residual (a junction's stiffness there is its critical current
        # against the drive), so it is resolved to the residual tolerance
        # over that stiffness rather than to machine precision; the voltages
        # and the scattering matrix are what the direct current determines
        @testset "the solvers agree" begin
            methods = Any[("Newton", Newton()), ("QuasiNewton", QuasiNewton()),
                ("NewtonKrylov", NewtonKrylov()),
                ("NewtonKrylov+Floquet", NewtonKrylov(preconditioner = Floquet())),
                ("Staged", Staged())]
            ref = hbnlsolve(case.wp, case.Npump, case.sources, case.circuit,
                case.defs; nlkw..., ftol = 1e-14, keyedarrays = false,
                method = Newton())
            @test ref.solverinfo.converged
            nm = length(ref.modes)
            isdc = [all(iszero, ref.modes[(k - 1) % nm + 1])
                for k in eachindex(ref.nodeflux)]
            acscale = maximum(abs, ref.nodeflux[.!isdc])
            for (name, m) in methods
                @testset "$name" begin
                    sol = hbnlsolve(case.wp, case.Npump, case.sources,
                        case.circuit, case.defs; nlkw..., ftol = 1e-14,
                        keyedarrays = false, method = m)
                    @test sol.solverinfo.converged
                    @test maximum(abs, sol.S .- ref.S) < 1e-8
                    @test maximum(abs, sol.nodeflux[.!isdc] .-
                        ref.nodeflux[.!isdc]) < 1e-7*acscale
                    # the real part only: QuasiNewton solves a complex flux
                    # at zero frequency too and its imaginary part there is
                    # spurious (see its docstring)
                    turns = real.(sol.nodeflux[isdc] .- ref.nodeflux[isdc]) ./
                        (2*pi)
                    @test all(x -> isapprox(x, round(x); atol = 1e-4), turns)
                    if nlkw.dc
                        @test isapprox(sol.dcnodevoltage, ref.dcnodevoltage;
                            rtol = 1e-6, atol = 1e-12)
                    end
                end
            end
        end
    end

    # C. keyed and plain outputs are the same numbers through the index maps
    # the documentation states: the mode index fastest within a port
    @testset "keyed and plain outputs" begin
        keyed = case.pumped ? solve(returnSnoise = true) :
            hblinsolve(case.ws, case.circuit, case.defs; returnSnoise = true)
        plain = case.pumped ? solve(returnSnoise = true, keyedarrays = false) :
            hblinsolve(case.ws, case.circuit, case.defs; returnSnoise = true,
                keyedarrays = false)
        kl = case.pumped ? keyed.linearized : keyed
        pl = case.pumped ? plain.linearized : plain
        modes = kl.modes
        ports = kl.portnumbers
        nm = length(modes)
        # the plain arrays read back through the keyed axes
        slots = [(m, p) for p in ports for m in modes]
        @test [kl.S(m, p, n, q, i) for (m, p) in slots, (n, q) in slots,
            i in eachindex(case.ws)] == pl.S
        @test [kl.QE(m, p, n, q, i) for (m, p) in slots, (n, q) in slots,
            i in eachindex(case.ws)] == pl.QE
        @test [kl.CM(m, p, i) for (m, p) in slots, i in eachindex(case.ws)] == pl.CM
        @test size(pl.Snoise, 3) == length(case.ws)
        @test vec(Array(kl.Snoise)) == vec(pl.Snoise)
    end

    # D. what is true of every scattering matrix: the commutation relations
    # close to the sign of the mode frequency, the quantum efficiency is
    # bounded by the ideal amplifier's, a lossless circuit has no noise
    # channels, a passive lossless one is unitary and reciprocal, and a
    # passive lossy one has Bosma's noise covariance
    @testset "properties of the scattering matrix" begin
        out = case.pumped ? solve(returnSnoise = true, returnCnoise = true,
                keyedarrays = false).linearized :
            hblinsolve(case.ws, case.circuit, case.defs; returnSnoise = true,
                returnCnoise = true, keyedarrays = false)
        modes = out.modes; ports = out.portnumbers; nm = length(modes)
        for i in eachindex(case.ws), (mo, m) in enumerate(modes),
                (po, p) in enumerate(ports)
            f = sum(m .* case.wp) + case.ws[i]
            @test isapprox(out.CM[(po-1)*nm + mo, i], sign(f); atol = 1e-8)
        end
        @test all(q -> -1e-12 <= q <= 1 + 1e-12, out.QE)
        @test all(out.QE .<= out.QEideal .+ 1e-12)
        if case.lossless
            @test isempty(out.Snoise)
            @test all(iszero, out.Cnoise)
        else
            @test !isempty(out.Snoise)
            @test all(isfinite, out.Snoise)
        end
        if case.linear
            for i in eachindex(case.ws)
                S = out.S[:, :, i]
                @test isapprox(S, transpose(S); atol = 1e-12)
                if case.lossless
                    @test isapprox(S'*S, I; atol = 1e-10)
                else
                    @test isapprox(out.Cnoise[:, :, i],
                        JosephsonCircuits.calcCnoise(S); atol = 1e-10)
                end
            end
        end
    end
end
end

# E. the noise against the connection algebra: hblinsolve's scattering
# and added noise of a circuit of blocks against solveS cascading the
# blocks' own scattering matrices and Bosma covariances, first passive,
# then with a pumped amplifier whose multi mode (S, Cnoise) is one of the
# networks
@testset "noise against the cascade" begin
    Z0 = 50.0
    attS(dB) = JosephsonCircuits.ABCDtoS(
        JosephsonCircuits.ABCD_attenuator_T(Z0, dB); portimpedances = Z0)
    bosma(S) = JosephsonCircuits.calcCnoise(S)
    # a two port block over nm modes: the same two port at every mode, in
    # the port major, mode fastest order of both hblinsolve and solveS
    function multimode(S2, nm)
        M = zeros(ComplexF64, 2*nm, 2*nm)
        for p in 1:2, q in 1:2, m in 1:nm
            M[(p-1)*nm + m, (q-1)*nm + m] = S2[p, q]
        end
        return M
    end
    ws = 2*pi*[4.6e9, 5.0e9]

    @testset "two attenuators, passive" begin
        Sa, Sb = attS(3.0), attS(6.0)
        circuit = Circuit([(:p1, 1, 0, Port(1; Z0 = Z0)),
            (:a, 1, 2, ScatteringParameters(Sa; nports = 2, zref = Z0)),
            (:b, 2, 3, ScatteringParameters(Sb; nports = 2, zref = Z0)),
            (:p2, 3, 0, Port(2; Z0 = Z0))])
        out = hblinsolve(ws, circuit; returnSnoise = true, returnCnoise = true,
            keyedarrays = false)
        # the two blocks connected by the algebra, with their Bosma covariances
        cascade = JosephsonCircuits.solveS([("a", Sa), ("b", Sb)],
            [[("a", 2), ("b", 1)]]; noise = true)
        @test cascade.ports == [("a", 1), ("b", 2)]
        Sc, Cc = JosephsonCircuits.interconnectS(Sa, Sb, bosma(Sa), bosma(Sb), 2, 1)
        @test isapprox(Matrix(cascade.S), Sc; atol = 1e-12)
        @test isapprox(Matrix(cascade.C), Cc; atol = 1e-12)
        for i in eachindex(ws)
            @test isapprox(out.S[:, :, i], Sc; atol = 1e-10)
            @test isapprox(out.Cnoise[:, :, i], Cc; atol = 1e-10)
            # which is Bosma's theorem for the cascade itself
            @test isapprox(out.Cnoise[:, :, i], bosma(out.S[:, :, i]); atol = 1e-10)
        end
    end

    @testset "a resistor's noise channel is a port in vacuum" begin
        # the lossy JPA, and the same circuit with its resistor replaced by
        # a port of that impedance: the noise scattering matrix must be
        # that port's column of S, and the added noise covariance that
        # port's vacuum through S
        lossy = CASES[findfirst(c -> c.name == "4WM JPA with loss", CASES)]
        R = 2.0e5
        ported = Circuit([(:p1, 1, 0, Port(1; Z0 = Z0)),
            (:cc, 1, 2, Capacitor(100e-15)),
            (:jj, 2, 0, JosephsonJunction(1000e-12)),
            (:cj, 2, 0, Capacitor(1000e-15)), (:p2, 2, 0, Port(2; Z0 = R))])
        a = hbsolve(ws, lossy.wp, lossy.sources, lossy.Nmod, lossy.Npump,
            lossy.circuit, lossy.defs; lossy.kw..., ftol = 1e-14,
            returnSnoise = true, returnCnoise = true, keyedarrays = false)
        b = hbsolve(ws, lossy.wp, lossy.sources, lossy.Nmod, lossy.Npump,
            ported, lossy.defs; lossy.kw..., ftol = 1e-14,
            keyedarrays = false)
        nm = length(a.linearized.modes)
        @test size(a.linearized.Snoise) == (nm, nm, length(ws))
        for i in eachindex(ws)
            Sp = b.linearized.S[1:nm, nm+1:2*nm, i]
            @test isapprox(transpose(a.linearized.Snoise[:, :, i]), Sp;
                atol = 1e-10)
            @test isapprox(a.linearized.Cnoise[:, :, i], Sp*Sp'; atol = 1e-10)
            @test isapprox(a.linearized.S[:, :, i], b.linearized.S[1:nm, 1:nm, i];
                atol = 1e-10)
        end
    end

    @testset "a pumped JPA behind an attenuator" begin
        # the circuit with the block inside it
        case = CASES[findfirst(c -> c.name == "scattering block JPA", CASES)]
        total = hbsolve(ws, case.wp, case.sources, case.Nmod, case.Npump,
            case.circuit, case.defs; case.kw..., ftol = 1e-14,
            returnSnoise = true, returnCnoise = true, keyedarrays = false)
        Stot, Ctot = total.linearized.S, total.linearized.Cnoise
        nm = length(total.linearized.modes)
        # the amplifier alone, pumped by the wave the matched attenuator
        # passes: the attenuator reflects nothing back, so the incident
        # pump at the amplifier is the transmission times the incident
        # pump at the port, and a current source on a matched port launches
        # a wave proportional to its current
        Satt = attS(3.0)
        tau = Satt[2, 1]
        jpa = CASES[findfirst(c -> c.name == "4WM JPA", CASES)]
        Ip = only(case.sources).current
        alone = hbsolve(ws, case.wp, [(mode = (1,), port = 1, current = tau*Ip)],
            case.Nmod, case.Npump, jpa.circuit, jpa.defs; jpa.kw...,
            ftol = 1e-14, returnCnoise = true, keyedarrays = false)
        @test alone.linearized.modes == total.linearized.modes
        Sm, Cm = multimode(Satt, nm), bosma(multimode(Satt, nm))
        for i in eachindex(ws)
            Sj = alone.linearized.S[:, :, i]
            Cj = zeros(ComplexF64, nm, nm)
            @test all(iszero, alone.linearized.Cnoise[:, :, i])
            cascade = JosephsonCircuits.solveS([("jpa", Sj, Cj), ("att", Sm, Cm)],
                [[("att", 2), ("jpa", 1)]]; noise = true, Nmodes = nm)
            # the external port is the attenuator's input, all of its modes
            @test length(cascade.ports) == nm
            @test all(p -> p[1] == "att", cascade.ports)
            @test isapprox(Matrix(cascade.S), Stot[:, :, i];
                atol = 1e-8*maximum(abs, Stot[:, :, i]))
            @test isapprox(Matrix(cascade.C), Ctot[:, :, i];
                atol = 1e-8*maximum(abs, Ctot[:, :, i]))
        end
    end
end

# F. analytic references: circuits small enough that the harmonic balance
# solution has a closed form to compare with
@testset "analytic references" begin
    phi0 = JosephsonCircuits.phi0
    Z0 = 50.0

    @testset "the harmonics of an imposed flux are Bessel functions" begin
        # a junction shunted by a small inductor and driven through a port
        # sees an imposed flux phi_p cos(w t): the shunt carries the drive,
        # the junction responds with Ic sin(phi), whose harmonics are the
        # Jacobi-Anger coefficients (-1)^k J_{2k+1}(phi_p), and those
        # currents divide between the shunt and the port resistor. The
        # relation is exact to first order in eps = Ls Ic/phi0, the ratio of
        # the junction's back action to the drive, and the error of every
        # harmonic is of that order on the scale of the largest one, so the
        # harmonics are compared to ten times eps times the largest
        Lj = 1.0e-9; Ls = 1.0e-13
        Ic = phi0/Lj
        eps = Ls*Ic/phi0
        circuit = Circuit([(:p1, 1, 0, Port(1; Z0 = Z0)),
            (:ls, 1, 0, Inductor(Ls)), (:jj, 1, 0, JosephsonJunction(Lj))])
        w = 2*pi*5.0e9
        for phip in (0.4, 1.2)
            Ip = phip*phi0/Ls
            sol = hbnlsolve((w,), (7,), [(mode = (1,), port = 1, current = Ip)],
                circuit, Dict{Any,Any}(); ftol = 1e-12, keyedarrays = false)
            @test sol.solverinfo.converged
            modes = sol.modes
            A(n) = sol.nodeflux[findfirst(==((n,)), modes)]
            # the fundamental is the imposed flux, phi(t) = 2 Re(A1 e^{iwt})
            # (a source of amplitude Ip drives a current 2 Re(Ip e^{iwt}), so
            # the imposed flux is 2 phip)
            A1 = A(1)
            @test isapprox(2*abs(A1), 2*phip; rtol = 10*eps)
            theta = angle(A1)
            scale = maximum(abs(A(n)) for n in (3, 5, 7))
            for (k, n) in enumerate((3, 5, 7))
                cn = (-1)^k*besselj(n, 2*abs(A1))*exp(im*n*theta)
                An = -Ic*cn/(phi0*(1/Ls + im*n*w/Z0))
                @test abs(A(n) - An) <= 10*eps*scale
            end
        end
    end

    @testset "the flux biased SQUID linearizes to its bias dependent inductance" begin
        # the flux biased SQUID of the table with the pump turned down to
        # nothing: the direct current phases of its two junctions set their
        # small signal inductances Lj/cos(phi), the loop inductor is loaded
        # by the flux line through the mutual inductance, and the port sees
        # the resulting resonator through the coupling capacitor
        case = CASES[findfirst(c -> c.name == "flux biased SQUID", CASES)]
        Lj = 400e-12; Cj = 0.6e-12; Ll = 40e-12; Ldc = 200e-12; K = 0.9
        Cc = 60e-15; R2 = 1000.0
        M = K*sqrt(Ll*Ldc)
        sources = [(mode = (0,), port = 2, current = 8.0e-6),
            (mode = (1,), port = 1, current = 1e-15)]
        ws = case.ws
        sol = hbsolve(ws, case.wp, sources, case.Nmod, case.Npump, case.circuit,
            case.defs; case.kw..., ftol = 1e-14, keyedarrays = false)
        @test sol.nonlinear.solverinfo.converged
        nm = length(sol.nonlinear.modes)
        dc = findfirst(==((0,)), sol.nonlinear.modes)
        phi2 = real(sol.nonlinear.nodeflux[(2-1)*nm + dc])   # node 2, junction 1
        phi3 = real(sol.nonlinear.nodeflux[(3-1)*nm + dc])   # node 3, junction 2
        @test abs(phi2) > 0.1 && abs(phi3) > 0.1
        nl = length(sol.linearized.modes)
        s0 = findfirst(==((0,)), sol.linearized.modes)
        for (i, w) in enumerate(ws)
            Z3 = 1/(im*w*Cj + cos(phi3)/(im*w*Lj))
            Zloop = im*w*Ll + (w*M)^2/(im*w*Ldc + R2) + Z3
            Y2 = im*w*Cj + cos(phi2)/(im*w*Lj) + 1/Zloop
            Zin = 1/(im*w*Cc) + 1/Y2
            S11 = (Zin - Z0)/(Zin + Z0)
            @test isapprox(sol.linearized.S[s0, s0, i], S11; atol = 1e-8)
        end
    end
end

end
