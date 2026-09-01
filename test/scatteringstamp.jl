using Symbolics
using JosephsonCircuits
using LinearAlgebra
using Test

@testset verbose=true "scattering parameter stamps" begin

    @testset "hybrid coefficient evaluation" begin
        # a shunt capacitor as a one port block: the constitutive equation
        # B*v = C*i must encode i = im*w*C*v, and an ideal short (S = -1,
        # where C = 0 and no admittance exists) must still stamp
        C = 300e-15
        Z0 = 50.0
        f(w) = fill((1 - im*w*C*Z0)/(1 + im*w*C*Z0), 1, 1)
        blk = ScatteringParameters(f; nports = 1, grounded = true)
        ws = [0.0, 2*pi*5e9, -2*pi*5e9]
        B = zeros(Complex{Float64}, 1, 1, 3)
        Cm = zeros(Complex{Float64}, 1, 1, 3)
        JosephsonCircuits.evaluatehybrid!(B, Cm, blk, ws)
        # zero frequency rows are i = 0
        @test B[1,1,1] == 0 && Cm[1,1,1] == 1
        @test isapprox(B[1,1,2]/Cm[1,1,2], im*ws[2]*C; rtol = 1e-12)
        # the negative frequency rule gives the conjugate
        @test isapprox(B[1,1,3]/Cm[1,1,3], conj(B[1,1,2]/Cm[1,1,2]);
            rtol = 1e-12)
        # an ideal short: B = 2/sqrt(Z0) forces v = 0, C = 0
        short = ScatteringParameters(w -> fill(-1.0 + 0.0im, 1, 1);
            nports = 1, grounded = true)
        JosephsonCircuits.evaluatehybrid!(B, Cm, short,
            [2*pi*5e9, 0.0, 2*pi*1e9])
        @test isapprox(B[1,1,1], 2/sqrt(Z0); rtol = 1e-12)
        @test Cm[1,1,1] == 0
        # a lossy attenuator stamps its own coefficients, dissipation and
        # all; what its loss costs in noise is planned elsewhere
        att = ScatteringParameters([0.0 0.5; 0.5 0.0])
        B2 = zeros(Complex{Float64}, 2, 2, 1)
        C2m = zeros(Complex{Float64}, 2, 2, 1)
        JosephsonCircuits.evaluatehybrid!(B2, C2m, att, [2*pi*5e9])
        @test isapprox(B2[1,2,1], -0.5/sqrt(Z0); rtol = 1e-12)
        @test isapprox(C2m[1,2,1], 0.5*sqrt(Z0); rtol = 1e-12)
    end

    # the JPA from the README with the shunt capacitor C2 expressed as a
    # one port grounded scattering block. the stamped system must reproduce
    # the lumped circuit through the full nonlinear + linearized solve.
    function jpacircuit(shunt)
        return Circuit(
            [:p1 => Port(1; termination = nothing), :r1 => Resistor(50.0),
             :cc => Capacitor(100.0e-15),
             :jj => JosephsonJunction(1000.0e-12),
             :c2 => shunt],
            [((:p1, 1), (:r1, 1), (:cc, 1)),
             ((:cc, 2), (:jj, 1), (:c2, 1)),
             ((:jj, 2), (:c2, 2), (:r1, 2), (:p1, 2), Ground)],
        )
    end
    capS(C, Z0) = w -> fill((1 - im*w*C*Z0)/(1 + im*w*C*Z0), 1, 1)

    ws = 2*pi*(4.5:0.01:5.0)*1e9
    wp = (2*pi*4.75001*1e9,)
    sources = [(mode=(1,), port=1, current=0.00565e-6)]

    @testset "one port block reproduces a lumped capacitor in hbsolve" begin
        C2 = 1000.0e-15
        lumped = jpacircuit(Capacitor(C2))
        # a grounded one port block has a single scalar terminal, so wrap
        # the connection endpoints accordingly: rebuild with the block
        blk = ScatteringParameters(capS(C2, 50.0); nports = 1, grounded = true)
        stamped = Circuit(
            [:p1 => Port(1; termination = nothing), :r1 => Resistor(50.0),
             :cc => Capacitor(100.0e-15),
             :jj => JosephsonJunction(1000.0e-12),
             :c2 => blk],
            [((:p1, 1), (:r1, 1), (:cc, 1)),
             ((:cc, 2), (:jj, 1), (:c2, 1)),
             ((:jj, 2), (:r1, 2), (:p1, 2), Ground)],
        )
        sol1 = hbsolve(ws, wp, sources, (8,), (16,), lumped)
        sol2 = hbsolve(ws, wp, sources, (8,), (16,), stamped)
        S1 = sol1.linearized.S((0,), 1, (0,), 1, :)
        S2 = sol2.linearized.S((0,), 1, (0,), 1, :)
        @test isapprox(S1, S2; rtol = 1e-8)
        @test maximum(abs2.(S1)) > 1.0 # it is still an amplifier
        # the pump operating point also agrees
        @test isapprox(sol1.nonlinear.nodeflux[:],
            sol2.nonlinear.nodeflux[:]; rtol = 1e-8)
        # the quantum efficiency and commutation relations require the
        # adjoint solutions, which the scattering path assembles and
        # factorizes independently instead of reusing the transposed
        # forward factorization; they must agree with the lumped circuit
        @test isapprox(sol1.linearized.QE((0,), 1, (0,), 1, :),
            sol2.linearized.QE((0,), 1, (0,), 1, :); rtol = 1e-8)
        @test isapprox(sol1.linearized.CM((0,), 1, :),
            sol2.linearized.CM((0,), 1, :); rtol = 1e-6)
    end

    @testset "two port floating block reproduces a series capacitor" begin
        # the coupling capacitor as a floating two port series element:
        # S for a series impedance z = 1/(im*w*C*Z0):
        # S11 = S22 = z/(z+2), S21 = S12 = 2/(z+2)
        Cc = 100.0e-15
        Z0 = 50.0
        function seriesS(w)
            z = 1/(im*w*Cc*Z0)
            return [z/(z+2) 2/(z+2); 2/(z+2) z/(z+2)]
        end
        blk = ScatteringParameters(seriesS; nports = 2, grounded = false)
        lumped = jpacircuit(Capacitor(1000.0e-15))
        stamped = Circuit(
            [:p1 => Port(1; termination = nothing), :r1 => Resistor(50.0),
             :cc => blk,
             :jj => JosephsonJunction(1000.0e-12),
             :c2 => Capacitor(1000.0e-15)],
            [((:p1, 1), (:r1, 1), (:cc, 1, 1)),
             ((:cc, 2, 1), (:jj, 1), (:c2, 1)),
             ((:cc, 1, 2), (:cc, 2, 2), (:jj, 2), (:c2, 2), (:r1, 2),
              (:p1, 2), Ground)],
        )
        sol1 = hbsolve(ws, wp, sources, (8,), (16,), lumped)
        sol2 = hbsolve(ws, wp, sources, (8,), (16,), stamped)
        @test isapprox(sol1.linearized.S((0,), 1, (0,), 1, :),
            sol2.linearized.S((0,), 1, (0,), 1, :); rtol = 1e-8)
    end

    @testset "tabulated block: interpolation and extrapolation policies" begin
        C2 = 1000.0e-15
        Z0 = 50.0
        lumped = jpacircuit(Capacitor(C2))
        # a dense grid covering signal band and pump harmonics
        fgrid = collect(range(0.05e9, 90e9, 24000))
        wgrid = 2*pi*fgrid
        Sgrid = zeros(Complex{Float64}, 1, 1, length(wgrid))
        for (k, w) in enumerate(wgrid)
            Sgrid[1,1,k] = (1 - im*w*C2*Z0)/(1 + im*w*C2*Z0)
        end
        maketab(extrapolation) = Circuit(
            [:p1 => Port(1; termination = nothing), :r1 => Resistor(50.0),
             :cc => Capacitor(100.0e-15),
             :jj => JosephsonJunction(1000.0e-12),
             :c2 => ScatteringParameters((wgrid, Sgrid); grounded = true,
                 extrapolation = extrapolation)],
            [((:p1, 1), (:r1, 1), (:cc, 1)),
             ((:cc, 2), (:jj, 1), (:c2, 1)),
             ((:jj, 2), (:r1, 2), (:p1, 2), Ground)],
        )
        sol1 = hbsolve(ws, wp, sources, (8,), (16,), lumped)
        sol2 = hbsolve(ws, wp, sources, (8,), (16,), maketab(:error))
        # agreement is limited by the linear interpolation of the
        # tabulated scattering data on the grid
        @test isapprox(sol1.linearized.S((0,), 1, (0,), 1, :),
            sol2.linearized.S((0,), 1, (0,), 1, :); rtol = 1e-3)
        # a short grid: the pump harmonics fall outside and the default
        # extrapolation policy raises an informative error
        wshort = 2*pi*collect(range(4e9, 6e9, 100))
        Sshort = zeros(Complex{Float64}, 1, 1, length(wshort))
        for (k, w) in enumerate(wshort)
            Sshort[1,1,k] = (1 - im*w*C2*Z0)/(1 + im*w*C2*Z0)
        end
        shortcircuit = Circuit(
            [:p1 => Port(1; termination = nothing), :r1 => Resistor(50.0),
             :cc => Capacitor(100.0e-15),
             :jj => JosephsonJunction(1000.0e-12),
             :c2 => ScatteringParameters((wshort, Sshort); grounded = true)],
            [((:p1, 1), (:r1, 1), (:cc, 1)),
             ((:cc, 2), (:jj, 1), (:c2, 1)),
             ((:jj, 2), (:r1, 2), (:p1, 2), Ground)],
        )
        @test_throws ArgumentError hbsolve(ws, wp, sources, (8,), (16,),
            shortcircuit)
    end

    @testset "scattering blocks inside a JTWPA" begin
        # a 20 cell JTWPA with every ground capacitor expressed as a one
        # port scattering block; hbsolve including pump, signal, and idlers
        # must match the lumped version
        Lj = IctoLj(3.4e-6)
        Cg = 45.0e-15
        Cj = 55e-15
        Nj = 20
        # `lumped = true` uses two terminal capacitors (whose second
        # terminals join the ground group); otherwise grounded one port
        # scattering blocks (single scalar terminal, reference auto-tied)
        function twpa(lumped::Bool)
            shuntfor(C) = lumped ? Capacitor(C) :
                ScatteringParameters(capS(C, 50.0); nports = 1, grounded = true)
            components = Any[:p1 => Port(1; termination = nothing), :r1 => Resistor(50.0),
                :cin => shuntfor(Cg/2)]
            connections = Any[]
            ground = Any[(:p1, 2), (:r1, 2)]
            lumped && push!(ground, (:cin, 2))
            prev = nothing
            for i in 1:Nj
                jid = Symbol(:jj, i)
                cid = Symbol(:cj, i)
                push!(components, jid => JosephsonJunction(Lj),
                    cid => Capacitor(Cj))
                left = prev === nothing ?
                    Any[(:p1, 1), (:r1, 1), (:cin, 1)] : Any[(prev, 2)]
                push!(left, (jid, 1), (cid, 1))
                if i > 1
                    gid = Symbol(:cg, i)
                    push!(components, gid => shuntfor(Cg))
                    push!(left, (gid, 1))
                    lumped && push!(ground, (gid, 2))
                end
                push!(connections, Tuple(left))
                prev = jid
                push!(connections, ((jid, 2), (cid, 2)))
            end
            push!(components, :cout => shuntfor(Cg/2),
                :r2 => Resistor(50.0), :p2 => Port(2; termination = nothing))
            push!(connections, ((prev, 2), (:cout, 1), (:r2, 1), (:p2, 1)))
            lumped && push!(ground, (:cout, 2))
            push!(ground, (:r2, 2), (:p2, 2), Ground)
            push!(connections, ground)
            return Circuit(components, connections)
        end
        lumped = twpa(true)
        stamped = twpa(false)
        wst = 2*pi*(3.0:0.5:9.0)*1e9
        wpt = (2*pi*7.12*1e9,)
        sourcest = [(mode=(1,), port=1, current=1.85e-6)]
        sol1 = hbsolve(wst, wpt, sourcest, (8,), (16,), lumped)
        sol2 = hbsolve(wst, wpt, sourcest, (8,), (16,), stamped)
        @test isapprox(sol1.linearized.S((0,), 2, (0,), 1, :),
            sol2.linearized.S((0,), 2, (0,), 1, :); rtol = 1e-8)
        @test isapprox(sol1.linearized.S((0,), 1, (0,), 1, :),
            sol2.linearized.S((0,), 1, (0,), 1, :); rtol = 1e-6)
    end

    @testset "transmission line stub" begin
        # an open circuited transmission line stub to ground, expressed two
        # ways: the two port TransmissionLine block with port 2 left open
        # (its floating node carries no current, which is the open boundary
        # condition), and a one port block with the analytic open stub
        # reflection S11 = exp(-2*im*w*tau). both must agree exactly.
        Z0 = 50.0
        len = 3e-3
        tau = len/JosephsonCircuits.speed_of_light
        tl = TransmissionLine(Z0, len; grounded = true)
        stub1 = ScatteringParameters(w -> fill(exp(-2*im*w*tau), 1, 1);
            nports = 1, grounded = true, negative_frequency = Native())
        function stubcircuit(stub, extra...)
            components = Any[:p1 => Port(1; termination = nothing), :r1 => Resistor(50.0),
                :l => Inductor(2e-9), :stub => stub]
            for (id, def) in extra
                push!(components, id => def)
            end
            return Circuit(components,
                [((:p1, 1), (:r1, 1), (:l, 1)),
                 ((:l, 2), (:stub, 1)),
                 ((:p1, 2), (:r1, 2), Ground)],
            )
        end
        wsl = 2*pi*(1.0:0.5:15.0)*1e9
        sol1 = hblinsolve(wsl, stubcircuit(tl))
        sol2 = hblinsolve(wsl, stubcircuit(stub1))
        @test isapprox(sol1.S((0,), 1, (0,), 1, :),
            sol2.S((0,), 1, (0,), 1, :); rtol = 1e-8)
        # exactly at the quarter wave resonance S11 = -1: the stub is an
        # ideal short, which has no admittance representation, but the
        # hybrid stamps handle it exactly: the port sees the inductor
        # shorted to ground
        wres = 2*pi*JosephsonCircuits.speed_of_light/(4*len)
        solres = hblinsolve([wres], stubcircuit(stub1))
        L = 2e-9
        zl = im*wres*L/50.0
        @test isapprox(solres.S((0,), 1, (0,), 1, 1),
            (zl - 1)/(zl + 1); rtol = 1e-8)
        @test isapprox(solres.S((0,), 1, (0,), 1, 1),
            hblinsolve([wres], stubcircuit(tl)).S((0,), 1, (0,), 1, 1);
            rtol = 1e-8)
    end

    @testset "through line across its half wave resonances" begin
        # a matched line between the port and a load: at theta = n*pi,
        # det(I+S) = 0 and the admittance representation does not exist,
        # but the response is trivial (|S21| = 1). sweep across several
        # half wave frequencies including one exactly on resonance.
        Z0 = 50.0
        len = 5e-3
        tau = len/JosephsonCircuits.speed_of_light
        fhalf = 1/(2*tau) # first half wave frequency
        line = TransmissionLine(Z0, len; grounded = true)
        c = Circuit(
            [:p1 => Port(1; termination = nothing), :r1 => Resistor(50.0), :line => line,
             :r2 => Resistor(50.0), :p2 => Port(2; termination = nothing)],
            [((:p1, 1), (:r1, 1), (:line, 1)),
             ((:line, 2), (:r2, 1), (:p2, 1)),
             ((:p1, 2), (:r1, 2), (:r2, 2), (:p2, 2), Ground)],
        )
        wsweep = 2*pi*[0.5*fhalf, fhalf, 1.5*fhalf, 2*fhalf, 2.5*fhalf]
        sol = hblinsolve(wsweep, c)
        S21 = sol.S((0,), 2, (0,), 1, :)
        S11 = sol.S((0,), 1, (0,), 1, :)
        @test isapprox(S21, exp.(-im*wsweep*tau); rtol = 1e-8)
        @test all(abs.(S11) .< 1e-8)
    end

    @testset "hbnlsolve pump path and a dissipative block" begin
        # the pump solve alone with a scattering block matches the lumped
        # circuit
        C2 = 1000.0e-15
        lumped = jpacircuit(Capacitor(C2))
        stamped = Circuit(
            [:p1 => Port(1; termination = nothing), :r1 => Resistor(50.0),
             :cc => Capacitor(100.0e-15),
             :jj => JosephsonJunction(1000.0e-12),
             :c2 => ScatteringParameters(capS(C2, 50.0); nports = 1,
                 grounded = true)],
            [((:p1, 1), (:r1, 1), (:cc, 1)),
             ((:cc, 2), (:jj, 1), (:c2, 1)),
             ((:jj, 2), (:r1, 2), (:p1, 2), Ground)],
        )
        out1 = hbnlsolve(wp, (8,), sources, lumped)
        out2 = hbnlsolve(wp, (8,), sources, stamped)
        @test isapprox(out1.nodeflux[:], out2.nodeflux[:]; rtol = 1e-8)

        # a dissipative block adds its vacuum noise, so a circuit which
        # holds one alongside a dissipative lumped component still closes
        # the commutation relations
        lossy = Circuit(
            [:p1 => Port(1; termination = nothing), :r1 => Resistor(50.0),
             :att => ScatteringParameters([0.0 0.5; 0.5 0.0]; grounded = true),
             :rt => Resistor(50.0)],
            [((:p1, 1), (:r1, 1), (:att, 1)),
             ((:att, 2), (:rt, 1)),
             ((:rt, 2), (:p1, 2), (:r1, 2), Ground)],
        )
        out = hblinsolve(2*pi*[5e9], lossy; keyedarrays = false)
        @test isapprox(out.CM[1,1], 1.0; atol = 1e-12)
        @test_logs hblinsolve(2*pi*[5e9], lossy)   # and warns about nothing
    end

    @testset "tabulated blocks are evaluated on the backend" begin
        # A callable provider is an arbitrary Julia function, so its values
        # are computed on the host. Tabulated and constant data are not: the
        # evaluation is a search and an interpolation, which a kernel does,
        # so those blocks never touch the host per frequency. The kernel runs
        # on CPU() unchanged, so this covers its arithmetic everywhere.
        Cval = 1000.0e-15
        ftab = 2*pi*collect(range(0.5e9, 30e9, length = 200))
        capval(w) = (1 - im*w*Cval*50.0)/(1 + im*w*Cval*50.0)
        for (lbl, scale) in (("lossless", 1.0), ("lossy", 0.85))
            tab = reshape([scale*capval(w) for w in ftab], 1, 1, :)
            blk = ScatteringParameters((ftab, tab); nports = 1, grounded = true,
                extrapolation = :constant)
            circuit = Circuit(
                Any[:p1 => Port(1; termination = nothing), :r1 => Resistor(50.0),
                    :cc => Capacitor(100.0e-15),
                    :jj => JosephsonJunction(1000.0e-12), :c2 => blk,
                    :rl => Resistor(1.0e5)],
                Any[((:p1,1), (:r1,1), (:cc,1)),
                    ((:cc,2), (:jj,1), (:c2,1), (:rl,1)),
                    ((:jj,2), (:r1,2), (:p1,2), Ground), ((:rl,2), Ground)])
            nl = hbnlsolve(wp, (16,), sources, circuit; keyedarrays = false)
            psc = JosephsonCircuits.compile(circuit)
            cg = JosephsonCircuits.calccircuitgraph(psc)
            sf = JosephsonCircuits.truncfreqs(
                JosephsonCircuits.calcfreqsdft((4,)); dc = true, odd = false,
                even = true, maxintermodorder = Inf)
            wsweep = 2*pi*collect(range(4.13e9, 5.27e9, length = 7))
            d = JosephsonCircuits.hblinsolve(wsweep, psc, cg,
                Dict{Symbol,Number}(), sf; nonlinear = nl, debuglsys = true)
            ssys = d.lsys.scattering
            @test JosephsonCircuits.candeviceevaluate(ssys)

            A = copy(d.lsys.Asparse)
            nz = JosephsonCircuits.SparseArrays.nnz(A)
            perm = JosephsonCircuits.cscvaluepermutation(A)
            host = Matrix{ComplexF64}(undef, nz, length(wsweep))
            for (i, w) in enumerate(wsweep)
                JosephsonCircuits.assemblesystemmatrix!(A, d.lsys,
                    w .+ d.wpumpmodes)
                host[:, i] .= JosephsonCircuits.SparseArrays.nonzeros(A)
            end
            plan, _, _ = JosephsonCircuits.planfrequencysweep(d.lsys,
                JosephsonCircuits.CPU(); adjoint = false)
            got = Matrix{ComplexF64}(undef, nz, length(wsweep))
            JosephsonCircuits.assemblesweep!(got, plan, wsweep)
            st = JosephsonCircuits.plandevicescattering(ssys,
                JosephsonCircuits.sweepdestinations(A, ssys.Aindex, false),
                nz, length(wsweep), JosephsonCircuits.CPU(), d.Nmodes)
            dp = JosephsonCircuits.plandeviceproviders(ssys, length(wsweep),
                JosephsonCircuits.CPU(), d.wpumpmodes, ssys.scale)
            @test !isnothing(dp)
            JosephsonCircuits.stagedeviceproviders!(st.values, dp, wsweep, 1,
                length(wsweep))
            JosephsonCircuits.applyscatteringstamps!(got, st)
            @test got[invperm(perm), :] == host
        end

        # a callable provider is not something a kernel can evaluate, and
        # must be recognised as such rather than silently mis-evaluated
        cb = ScatteringParameters(capS(Cval, 50.0); nports = 1, grounded = true)
        ccircuit = Circuit(
            Any[:p1 => Port(1; termination = nothing), :r1 => Resistor(50.0),
                :cc => Capacitor(100.0e-15),
                :jj => JosephsonJunction(1000.0e-12), :c2 => cb],
            Any[((:p1,1), (:r1,1), (:cc,1)), ((:cc,2), (:jj,1), (:c2,1)),
                ((:jj,2), (:r1,2), (:p1,2), Ground)])
        cnl = hbnlsolve(wp, (16,), sources, ccircuit; keyedarrays = false)
        cpsc = JosephsonCircuits.compile(ccircuit)
        ccg = JosephsonCircuits.calccircuitgraph(cpsc)
        csf = JosephsonCircuits.truncfreqs(
            JosephsonCircuits.calcfreqsdft((4,)); dc = true, odd = false,
            even = true, maxintermodorder = Inf)
        cd = JosephsonCircuits.hblinsolve(2*pi*[4.13e9, 5.27e9], cpsc, ccg,
            Dict{Symbol,Number}(), csf; nonlinear = cnl, debuglsys = true)
        @test !JosephsonCircuits.candeviceevaluate(cd.lsys.scattering)
    end

    @testset "the three callable forms give identical results" begin
        # A callable which returns a fresh matrix is the natural way to write
        # a block by hand and is the default. On a line whose every cell is a
        # block those matrices dominate the allocation of a sweep, so a
        # callable may instead write into a destination, or return one entry
        # at a time. The last is the only form a kernel can call, and so the
        # only one which lets a callable block be evaluated on a backend.
        # All three must give exactly the same answer.
        Cval = 1000.0e-15
        ret = ScatteringParameters(capS(Cval, 50.0); nports = 1, grounded = true)
        inp = ScatteringParameters((d, w) -> (d[1,1] =
                (1 - im*w*Cval*50.0)/(1 + im*w*Cval*50.0); d);
            nports = 1, grounded = true, form = :inplace)
        ent = ScatteringParameters((p, q, w) ->
                (1 - im*w*Cval*50.0)/(1 + im*w*Cval*50.0);
            nports = 1, grounded = true, form = :entry)
        wtest = 2*pi*[0.0, 4.13e9, -5.27e9]
        B1 = zeros(ComplexF64,1,1,3); C1 = zeros(ComplexF64,1,1,3)
        JosephsonCircuits.evaluatehybrid!(B1, C1, ret, wtest)
        for other in (inp, ent)
            B2 = similar(B1); C2m = similar(C1)
            JosephsonCircuits.evaluatehybrid!(B2, C2m, other, wtest)
            @test B1 == B2
            @test C1 == C2m
        end

        # and through the whole solve
        mk(b) = Circuit(
            Any[:p1 => Port(1; termination = nothing), :r1 => Resistor(50.0),
                :cc => Capacitor(100.0e-15),
                :jj => JosephsonJunction(1000.0e-12), :c2 => b],
            Any[((:p1,1), (:r1,1), (:cc,1)), ((:cc,2), (:jj,1), (:c2,1)),
                ((:jj,2), (:r1,2), (:p1,2), Ground)])
        s1 = hbsolve(ws, wp, sources, (8,), (8,), mk(ret); keyedarrays = false)
        for other in (inp, ent)
            s2 = hbsolve(ws, wp, sources, (8,), (8,), mk(other);
                keyedarrays = false)
            @test s1.linearized.S == s2.linearized.S
        end

        # the form says how a function is called, so it belongs to a
        # callable; anything else is a mistake worth naming rather than a
        # missing method, and an unknown form likewise
        @test_throws ArgumentError ScatteringParameters([0.0 0.5; 0.5 0.0];
            form = :entry)
        @test_throws ArgumentError ScatteringParameters(
            ([1.0e9, 2.0e9], zeros(ComplexF64,1,1,2)); nports = 1,
            grounded = true, form = :inplace)
        @test_throws ArgumentError ScatteringParameters(capS(Cval, 50.0);
            nports = 1, grounded = true, form = :elementwise)
    end

    @testset "entry-wise callables are evaluated on the backend" begin
        # An entry-wise callable is the one form a kernel can call. It needs
        # the closures to live in a device array, which needs them to capture
        # only numbers and to share a type: blocks built from one helper do,
        # which is how a generated circuit is written. Anything else falls
        # back to the host.
        Cval = 1000.0e-15
        mkentry(C) = (p, q, w) -> (1 - im*w*C*50.0)/(1 + im*w*C*50.0)
        blk = ScatteringParameters(mkentry(Cval); nports = 1, grounded = true,
            form = :entry)
        circuit = Circuit(
            Any[:p1 => Port(1; termination = nothing), :r1 => Resistor(50.0),
                :cc => Capacitor(100.0e-15),
                :jj => JosephsonJunction(1000.0e-12), :c2 => blk,
                :rl => Resistor(1.0e5)],
            Any[((:p1,1), (:r1,1), (:cc,1)),
                ((:cc,2), (:jj,1), (:c2,1), (:rl,1)),
                ((:jj,2), (:r1,2), (:p1,2), Ground), ((:rl,2), Ground)])
        nl = hbnlsolve(wp, (16,), sources, circuit; keyedarrays = false)
        psc = JosephsonCircuits.compile(circuit)
        cg = JosephsonCircuits.calccircuitgraph(psc)
        sf = JosephsonCircuits.truncfreqs(
            JosephsonCircuits.calcfreqsdft((4,)); dc = true, odd = false,
            even = true, maxintermodorder = Inf)
        wsweep = 2*pi*collect(range(4.13e9, 5.27e9, length = 7))
        d = JosephsonCircuits.hblinsolve(wsweep, psc, cg,
            Dict{Symbol,Number}(), sf; nonlinear = nl, debuglsys = true)
        ssys = d.lsys.scattering
        @test JosephsonCircuits.candeviceevaluate(ssys)

        A = copy(d.lsys.Asparse)
        nz = JosephsonCircuits.SparseArrays.nnz(A)
        perm = JosephsonCircuits.cscvaluepermutation(A)
        host = Matrix{ComplexF64}(undef, nz, length(wsweep))
        for (i, w) in enumerate(wsweep)
            JosephsonCircuits.assemblesystemmatrix!(A, d.lsys,
                w .+ d.wpumpmodes)
            host[:, i] .= JosephsonCircuits.SparseArrays.nonzeros(A)
        end
        plan, _, _ = JosephsonCircuits.planfrequencysweep(d.lsys,
            JosephsonCircuits.CPU(); adjoint = false)
        got = Matrix{ComplexF64}(undef, nz, length(wsweep))
        JosephsonCircuits.assemblesweep!(got, plan, wsweep)
        st = JosephsonCircuits.plandevicescattering(ssys,
            JosephsonCircuits.sweepdestinations(A, ssys.Aindex, false),
            nz, length(wsweep), JosephsonCircuits.CPU(), d.Nmodes)
        dp = JosephsonCircuits.plandeviceproviders(ssys, length(wsweep),
            JosephsonCircuits.CPU(), d.wpumpmodes, ssys.scale)
        @test !isnothing(dp)
        JosephsonCircuits.stagedeviceproviders!(st.values, dp, wsweep, 1,
            length(wsweep))
        JosephsonCircuits.applyscatteringstamps!(got, st)
        @test got[invperm(perm), :] == host

        # A two port block whose scattering matrix is not symmetric, so the
        # kernel reading S[q,p] where it should read S[p,q] is a failure
        # rather than a no-op. A one port block cannot see that at all.
        mkasym(a) = (p, q, w) -> (p == q ? 0.1*a : (p < q ? 0.2*a : 0.7*a)) *
            (1 - im*w*1e-12*50.0)/(1 + im*w*1e-12*50.0)
        two = ScatteringParameters(mkasym(1.0); nports = 2, grounded = false,
            form = :entry)
        tcircuit = Circuit(
            Any[:p1 => Port(1; termination = nothing), :r1 => Resistor(50.0), :cc => two,
                :jj => JosephsonJunction(1000.0e-12),
                :c2 => Capacitor(1000.0e-15)],
            Any[((:p1,1), (:r1,1), (:cc,1,1)),
                ((:cc,2,1), (:jj,1), (:c2,1)),
                ((:cc,1,2), (:cc,2,2), (:jj,2), (:c2,2), (:r1,2), (:p1,2),
                 Ground)])
        tnl = hbnlsolve(wp, (16,), sources, tcircuit; keyedarrays = false)
        tpsc = JosephsonCircuits.compile(tcircuit)
        tcg = JosephsonCircuits.calccircuitgraph(tpsc)
        td = JosephsonCircuits.hblinsolve(wsweep, tpsc, tcg,
            Dict{Symbol,Number}(), sf; nonlinear = tnl, debuglsys = true)
        tssys = td.lsys.scattering
        @test JosephsonCircuits.candeviceevaluate(tssys)
        tA = copy(td.lsys.Asparse)
        tnz = JosephsonCircuits.SparseArrays.nnz(tA)
        tperm = JosephsonCircuits.cscvaluepermutation(tA)
        thost = Matrix{ComplexF64}(undef, tnz, length(wsweep))
        for (i, w) in enumerate(wsweep)
            JosephsonCircuits.assemblesystemmatrix!(tA, td.lsys,
                w .+ td.wpumpmodes)
            thost[:, i] .= JosephsonCircuits.SparseArrays.nonzeros(tA)
        end
        tplan, _, _ = JosephsonCircuits.planfrequencysweep(td.lsys,
            JosephsonCircuits.CPU(); adjoint = false)
        tgot = Matrix{ComplexF64}(undef, tnz, length(wsweep))
        JosephsonCircuits.assemblesweep!(tgot, tplan, wsweep)
        tst = JosephsonCircuits.plandevicescattering(tssys,
            JosephsonCircuits.sweepdestinations(tA, tssys.Aindex, false),
            tnz, length(wsweep), JosephsonCircuits.CPU(), td.Nmodes)
        tdp = JosephsonCircuits.plandeviceproviders(tssys, length(wsweep),
            JosephsonCircuits.CPU(), td.wpumpmodes, tssys.scale)
        JosephsonCircuits.stagedeviceproviders!(tst.values, tdp, wsweep, 1,
            length(wsweep))
        JosephsonCircuits.applyscatteringstamps!(tgot, tst)
        @test tgot[invperm(tperm), :] == thost

        # a closure which captures something that is not a number cannot live
        # in a device array, so it stays on the host
        buf = [1.0]
        heavy = ScatteringParameters((p, q, w) ->
                (1 - im*w*buf[1]*1e-12*50.0)/(1 + im*w*buf[1]*1e-12*50.0);
            nports = 1, grounded = true, form = :entry)
        hcircuit = Circuit(
            Any[:p1 => Port(1; termination = nothing), :r1 => Resistor(50.0),
                :cc => Capacitor(100.0e-15),
                :jj => JosephsonJunction(1000.0e-12), :c2 => heavy],
            Any[((:p1,1), (:r1,1), (:cc,1)), ((:cc,2), (:jj,1), (:c2,1)),
                ((:jj,2), (:r1,2), (:p1,2), Ground)])
        hnl = hbnlsolve(wp, (16,), sources, hcircuit; keyedarrays = false)
        hpsc = JosephsonCircuits.compile(hcircuit)
        hcg = JosephsonCircuits.calccircuitgraph(hpsc)
        hd = JosephsonCircuits.hblinsolve(wsweep[1:3], hpsc, hcg,
            Dict{Symbol,Number}(), sf; nonlinear = hnl, debuglsys = true)
        @test !JosephsonCircuits.candeviceevaluate(hd.lsys.scattering)
    end

    @testset "a lossy non-reciprocal block: the adjoint is a transpose" begin
        # Every other block here is reciprocal, where S equals its own
        # transpose, so the forward system and the transposed one agree in
        # more ways than they do in general and a confusion between them is
        # invisible. An isolator is not reciprocal: it transmits from port
        # one to port two and blocks the reverse, and scaling it by less than
        # one makes it dissipative as a real isolator is. With one in the
        # circuit the two differ, and the adjoint the noise reads has to be
        # the transposed system: the noise a dissipative block adds is a
        # source in its auxiliary port current rows, and the response at an
        # output port to it is that source contracted against the transposed
        # solution driven at the port.
        a = 0.7
        for (lbl, blk) in (
                ("constant", ScatteringParameters(ComplexF64[0 0; a 0];
                    nports = 2, grounded = false)),
                ("entry", ScatteringParameters((p, q, w) ->
                        (p == 2 && q == 1) ? complex(a) : complex(0.0);
                    nports = 2, grounded = false, form = :entry)))
            circuit = Circuit(
                Any[:p1 => Port(1; termination = nothing), :r1 => Resistor(50.0), :iso => blk,
                    :jj => JosephsonJunction(1000.0e-12),
                    :cg => Capacitor(500.0e-15),
                    :p2 => Port(2; termination = nothing), :r2 => Resistor(50.0),
                    :cl => Capacitor(200.0e-15), :rl => Resistor(2.0e4)],
                Any[((:p1,1), (:r1,1), (:iso,1,1)),
                    ((:iso,2,1), (:jj,1), (:cg,1), (:p2,1), (:r2,1),
                     (:cl,1)),
                    ((:cl,2), (:rl,1)),
                    ((:iso,1,2), (:iso,2,2), (:jj,2), (:cg,2), (:r1,2),
                     (:p1,2), (:r2,2), (:p2,2), (:rl,2), Ground)])
            nl = hbnlsolve(wp, (16,), sources, circuit; keyedarrays = false)
            psc = JosephsonCircuits.compile(circuit)
            cg = JosephsonCircuits.calccircuitgraph(psc)
            sf = JosephsonCircuits.truncfreqs(
                JosephsonCircuits.calcfreqsdft((2,)); dc = true, odd = false,
                even = true, maxintermodorder = Inf)
            wsweep = 2*pi*collect(range(4.13e9, 5.27e9, length = 5))
            d = JosephsonCircuits.hblinsolve(wsweep, psc, cg,
                Dict{Symbol,Number}(), sf; nonlinear = nl, debuglsys = true)
            lsys = d.lsys
            @test !isnothing(lsys.scattering)

            # the circuit must actually be non-reciprocal and must actually
            # need the adjoint, or this tests nothing it is here for
            arrays = JosephsonCircuits.LinearizedArrays(; requestS = true,
                requestSnoise = false, requestSsensitivity = false,
                requestQE = true, requestCM = true, requestnodeflux = false,
                requestnodefluxadjoint = false, requestvoltage = false,
                requestvoltageadjoint = false, Nports = 2,
                Nmodes = d.Nmodes,
                Nnoisechannels = length(d.noiseportimpedanceindices),
                Ncomponents = 0, Nnodes = d.Nnodes,
                Nfrequencies = length(wsweep))
            @test JosephsonCircuits.needsadjointsolve(arrays,
                d.noiseportimpedanceindices)
            sol = hbsolve(wsweep, wp, sources, (16,), (2,), circuit)
            @test maximum(abs, sol.linearized.S((0,),2,(0,),1,:)) > 0.1
            @test maximum(abs, sol.linearized.S((0,),1,(0,),2,:)) < 1e-10

            # the forward system and its transpose, which for this block are
            # genuinely different matrices, each reproduced by the device
            # assembly: the adjoint direction of a sweep stamps the block
            # into the transposed structure, so its destinations differ from
            # the forward ones while its values do not
            A = copy(lsys.Asparse)
            nz = JosephsonCircuits.SparseArrays.nnz(A)
            perm = JosephsonCircuits.cscvaluepermutation(A)
            JosephsonCircuits.assemblesystemmatrix!(A, lsys,
                wsweep[1] .+ d.wpumpmodes)
            # a non-reciprocal block makes the system genuinely
            # non-symmetric, so its transpose is a different matrix
            @test A != JosephsonCircuits.SparseArrays.sparse(transpose(A))

            for isadjoint in (false, true)
                dest = JosephsonCircuits.sweepdestinations(lsys.Asparse,
                    lsys.scattering.Aindex, isadjoint)
                host = Matrix{ComplexF64}(undef, nz, length(wsweep))
                for (i, w) in enumerate(wsweep)
                    JosephsonCircuits.assemblesystemmatrix!(A, lsys,
                        w .+ d.wpumpmodes)
                    host[:, i] .=
                        JosephsonCircuits.SparseArrays.nonzeros(A)
                end
                plan, _, _ = JosephsonCircuits.planfrequencysweep(lsys,
                    JosephsonCircuits.CPU(); adjoint = isadjoint)
                got = Matrix{ComplexF64}(undef, nz, length(wsweep))
                JosephsonCircuits.assemblesweep!(got, plan, wsweep)
                st = JosephsonCircuits.plandevicescattering(lsys.scattering,
                    dest, nz, length(wsweep), JosephsonCircuits.CPU(),
                    d.Nmodes)
                dp = JosephsonCircuits.plandeviceproviders(lsys.scattering,
                    length(wsweep), JosephsonCircuits.CPU(), d.wpumpmodes,
                    lsys.scattering.scale)
                if isnothing(dp)
                    JosephsonCircuits.stagescatteringstamps!(st, wsweep, 1,
                        length(wsweep), d.wpumpmodes)
                else
                    JosephsonCircuits.stagedeviceproviders!(st.values, dp,
                        wsweep, 1, length(wsweep))
                end
                JosephsonCircuits.applyscatteringstamps!(got, st)
                # a sweep of the system itself stores its values in
                # compressed sparse row order; one of the transpose stores
                # them in the system's own order, because compressed sparse
                # row of the transpose is compressed sparse column of the
                # matrix, which is how the host holds it
                @test (isadjoint ? got : got[invperm(perm), :]) == host
            end

            # the commutation relations, which the isolator's own vacuum
            # noise closes and a confusion of the two systems does not
            @test all(x -> isapprox(abs(x), 1.0; atol = 1e-9),
                sol.linearized.CM)
        end
    end

    @testset "the device sweep assembly reproduces the host assembler" begin
        # hblinsolve on a backend does not call assemblesystemmatrix!: it
        # assembles a batch of frequencies with a kernel, because the terms
        # which are a constant quadratic in the signal frequency can all be
        # done at once. The scattering blocks are not such a term, since each
        # is evaluated through its own provider, so they are added afterwards
        # as a gather-add. This checks that the two together reproduce the
        # host assembler exactly.
        #
        # The kernels run on CPU() unchanged, so this covers the device path's
        # arithmetic on every machine; what it cannot cover is the device
        # solver itself.
        Cval = 1000.0e-15
        lossyS(C, Z0, a) = w -> fill(a*(1 - im*w*C*Z0)/(1 + im*w*C*Z0), 1, 1)
        for (lbl, blk) in (
                ("lossless", ScatteringParameters(capS(Cval, 50.0); nports = 1,
                    grounded = true)),
                ("lossy", ScatteringParameters(lossyS(Cval, 50.0, 0.85);
                    nports = 1, grounded = true)))
            circuit = Circuit(
                Any[:p1 => Port(1; termination = nothing), :r1 => Resistor(50.0),
                    :cc => Capacitor(100.0e-15),
                    :jj => JosephsonJunction(1000.0e-12), :c2 => blk,
                    :rl => Resistor(1.0e5)],
                Any[((:p1,1), (:r1,1), (:cc,1)),
                    ((:cc,2), (:jj,1), (:c2,1), (:rl,1)),
                    ((:jj,2), (:r1,2), (:p1,2), Ground), ((:rl,2), Ground)])
            nl = hbnlsolve(wp, (16,), sources, circuit; keyedarrays = false)
            psc = JosephsonCircuits.compile(circuit)
            cg = JosephsonCircuits.calccircuitgraph(psc)
            sf = JosephsonCircuits.truncfreqs(
                JosephsonCircuits.calcfreqsdft((4,)); dc = true, odd = false,
                even = true, maxintermodorder = Inf)
            wsweep = 2*pi*collect(range(4.13e9, 5.27e9, length = 7))
            d = JosephsonCircuits.hblinsolve(wsweep, psc, cg,
                Dict{Symbol,Number}(), sf; nonlinear = nl, debuglsys = true)
            lsys = d.lsys
            @test !isnothing(lsys.scattering)
            A = copy(lsys.Asparse)
            nz = JosephsonCircuits.SparseArrays.nnz(A)
            perm = JosephsonCircuits.cscvaluepermutation(A)
            # both directions of the sweep: the forward system and its
            # transpose, which is the adjoint the noise and quantum
            # efficiency calculations are read from
            for isadjoint in (false, true)
                dest = JosephsonCircuits.sweepdestinations(A,
                    lsys.scattering.Aindex, isadjoint)
                host = Matrix{ComplexF64}(undef, nz, length(wsweep))
                for (i, w) in enumerate(wsweep)
                    JosephsonCircuits.assemblesystemmatrix!(A, lsys,
                        w .+ d.wpumpmodes)
                    host[:, i] .= JosephsonCircuits.SparseArrays.nonzeros(A)
                end
                plan, _, _ = JosephsonCircuits.planfrequencysweep(lsys,
                    JosephsonCircuits.CPU(); adjoint = isadjoint)
                got = Matrix{ComplexF64}(undef, nz, length(wsweep))
                JosephsonCircuits.assemblesweep!(got, plan, wsweep)
                st = JosephsonCircuits.plandevicescattering(lsys.scattering,
                    dest, nz, length(wsweep), JosephsonCircuits.CPU(),
                    d.Nmodes)
                JosephsonCircuits.stagescatteringstamps!(st, wsweep, 1,
                    length(wsweep), d.wpumpmodes)
                JosephsonCircuits.applyscatteringstamps!(got, st)
                # the plan stores its values in compressed sparse row order
                # for the system and in the system's own order for the
                # transpose, whose compressed sparse row it already is
                @test (isadjoint ? got : got[invperm(perm), :]) == host
            end
        end
    end

    @testset "a dissipative block adds its vacuum noise" begin
        # A block which absorbs must add noise, or its output would violate
        # the commutation relations. The added noise wave has the vacuum
        # covariance I - S S', which enters the hybrid stamp as a source in
        # the auxiliary port current rows; the sharp test of the whole
        # construction is that |S|^2 + |Snoise|^2 comes back to one, which is
        # what calccm! computes.
        Z0 = 50.0
        wsn = 2*pi*[5.0e9]
        oneport(b) = Circuit(
            Any[:p1 => Port(1; termination = nothing), :r1 => Resistor(50.0),
                :cc => Capacitor(100.0e-15),
                :l1 => Inductor(1000.0e-12), :x => b],
            Any[((:p1,1), (:r1,1), (:cc,1)),
                ((:cc,2), (:l1,1), (:x,1)),
                ((:l1,2), (:r1,2), (:p1,2), Ground)])

        # every loss level, from lossless through fully absorbing to a sign
        # reversed reflection
        for Sv in (1.0, 0.99, 0.9, 0.7, 0.5, 0.0, -0.5, -0.9)
            blk = ScatteringParameters(fill(complex(Sv), 1, 1); grounded = true)
            out = hblinsolve(wsn, oneport(blk); keyedarrays = false,
                returnSnoise = true, returnCM = true)
            @test isapprox(out.CM[1,1], 1.0; atol = 1e-13)
            # a lossless block contributes no noise at all
            if abs(Sv) == 1
                @test iszero(sum(abs2, out.Snoise[:,1,1]))
            else
                @test sum(abs2, out.Snoise[:,1,1]) > 0
            end
        end

        # the same physical load written at three reference impedances is
        # one element, and its noise cannot depend on the bookkeeping. The
        # resistor, whose noise the package already computes, is the
        # reference.
        R = 200.0
        ref = hblinsolve(wsn, Circuit(
                Any[:p1 => Port(1; termination = nothing), :r1 => Resistor(50.0),
                    :cc => Capacitor(100.0e-15),
                    :l1 => Inductor(1000.0e-12), :x => Resistor(R)],
                Any[((:p1,1), (:r1,1), (:cc,1)),
                    ((:cc,2), (:l1,1), (:x,1)),
                    ((:l1,2), (:r1,2), (:p1,2), (:x,2), Ground)]);
            keyedarrays = false, returnSnoise = true)
        for zr in (50.0, 200.0, 377.0)
            blk = ScatteringParameters(fill(complex((R-zr)/(R+zr)), 1, 1);
                zref = zr, grounded = true)
            out = hblinsolve(wsn, oneport(blk); keyedarrays = false,
                returnSnoise = true)
            @test isapprox(sum(abs2, out.Snoise[:,1,1]),
                sum(abs2, ref.Snoise[:,1,1]); rtol = 1e-10)
        end

        # a two port: a series resistor, whose noise the package computes
        # from the resistor itself, against the same element as a block
        twoportR(x) = Circuit(
            Any[:p1 => Port(1; termination = nothing), :r1 => Resistor(50.0), :x => x,
                :l1 => Inductor(1000.0e-12), :p2 => Port(2; termination = nothing),
                :r2 => Resistor(50.0)],
            Any[((:p1,1), (:r1,1), (:x,1)),
                ((:x,2), (:l1,1), (:p2,1), (:r2,1)),
                ((:l1,2), (:r1,2), (:p1,2), (:r2,2), (:p2,2), Ground)])
        twoportB(b) = Circuit(
            Any[:p1 => Port(1; termination = nothing), :r1 => Resistor(50.0), :x => b,
                :l1 => Inductor(1000.0e-12), :p2 => Port(2; termination = nothing),
                :r2 => Resistor(50.0)],
            Any[((:p1,1), (:r1,1), (:x,1,1)),
                ((:x,2,1), (:l1,1), (:p2,1), (:r2,1)),
                ((:x,1,2), (:x,2,2), (:l1,2), (:r1,2), (:p1,2), (:r2,2),
                 (:p2,2), Ground)])
        for R in (50.0, 200.0)
            z = R/Z0; sd = z/(z+2); so = 2/(z+2)
            a = hblinsolve(wsn, twoportR(Resistor(R)); keyedarrays = false,
                returnSnoise = true)
            b = hblinsolve(wsn,
                twoportB(ScatteringParameters(ComplexF64[sd so; so sd];
                    grounded = false));
                keyedarrays = false, returnSnoise = true, returnCM = true)
            # the series element is transparent in one direction, so I - S is
            # singular and only a formulation which never inverts it works
            for i in 1:2
                @test isapprox(sum(abs2, b.Snoise[:,i,1]),
                    sum(abs2, a.Snoise[:,i,1]); rtol = 1e-10)
                @test isapprox(b.CM[i,1], 1.0; atol = 1e-12)
            end
        end

        # a lossy isolator: non-reciprocal, so the covariance is not the
        # covariance of the transposed block and the adjoint is not the
        # forward solution. A reciprocal block cannot see either mistake.
        for a in (1.0, 0.8, 0.5, 0.2)
            for Siso in (ComplexF64[0 0; a 0], ComplexF64[0 a; 0 0])
                out = hblinsolve(wsn, twoportB(ScatteringParameters(Siso; grounded = false));
                    keyedarrays = false, returnCM = true)
                @test isapprox(out.CM[1,1], 1.0; atol = 1e-12)
                @test isapprox(out.CM[2,1], 1.0; atol = 1e-12)
            end
        end

        # a covariance which is Hermitian but not symmetric, which
        # distinguishes the factor L from its conjugate
        for Sm in (ComplexF64[0.3im 0; 0.8 0.2],
                   ComplexF64[0.2 0.5im; 0.1 0.6],
                   ComplexF64[0.1+0.4im 0.2; 0.7im 0.3])
            V = I - Sm*Sm'
            @test !isapprox(V, transpose(V))
            out = hblinsolve(wsn, twoportB(ScatteringParameters(Sm; grounded = false));
                keyedarrays = false, returnCM = true)
            @test isapprox(out.CM[1,1], 1.0; atol = 1e-12)
            @test isapprox(out.CM[2,1], 1.0; atol = 1e-12)
        end
    end

    @testset "block noise channels: which blocks carry them and how they are named" begin
        # a lossless block's channels are identically zero, so a block whose
        # stored data can be shown to be unitary is given none. A callable is
        # opaque, so it is given channels whether it needs them or not, which
        # costs work and never correctness.
        Z0 = 50.0
        @test JosephsonCircuits.provablylossless(
            ScatteringParameters(ComplexF64[0 1; 1 0]))
        @test !JosephsonCircuits.provablylossless(
            ScatteringParameters(ComplexF64[0 0.5; 0.5 0]))
        ftab = 2*pi*collect(range(1e9, 10e9, length = 8))
        unit = reshape([complex(cis(w*1e-11)) for w in ftab], 1, 1, :)
        @test JosephsonCircuits.provablylossless(
            ScatteringParameters((ftab, unit); nports = 1, grounded = true))
        @test !JosephsonCircuits.provablylossless(
            ScatteringParameters((ftab, 0.9*unit); nports = 1, grounded = true))
        @test !JosephsonCircuits.provablylossless(
            ScatteringParameters(w -> fill(complex(1.0), 1, 1); nports = 1,
                grounded = true))

        # the rows of Snoise: the dissipative lumped components first, then
        # one channel per port of each block which carries them
        c = Circuit(
            Any[:p1 => Port(1; termination = nothing), :r1 => Resistor(50.0),
                :lossless => ScatteringParameters(ComplexF64[0 1; 1 0];
                    grounded = false),
                :att => ScatteringParameters(ComplexF64[0 0.5; 0.5 0];
                    grounded = false),
                :rt => Resistor(75.0)],
            Any[((:p1,1), (:r1,1), (:lossless,1,1)),
                ((:lossless,2,1), (:att,1,1)),
                ((:att,2,1), (:rt,1)),
                ((:lossless,1,2), (:lossless,2,2), (:att,1,2), (:att,2,2),
                 (:rt,2), (:p1,2), (:r1,2), Ground)])
        out = hblinsolve(2*pi*[5e9], c; returnSnoise = true)
        names = collect(JosephsonCircuits.AxisKeys.axiskeys(out.Snoise, 2))
        # the resistor, then the two ports of the dissipative block, and
        # nothing at all for the lossless one
        @test names == ["rt", "att/port1#1", "att/port1#2"]
        @test isapprox(out.CM[1,1], 1.0; atol = 1e-12)
    end

    @testset "a block declared lossless carries no noise channels" begin
        # `Passive` gives a block channels unless its data can be shown to be
        # unitary, which a callable's cannot. A lossless callable therefore
        # carries channels which are identically zero, and computing them was
        # measured at a third of the run time of a five hundred cell line on
        # the host and half of it on a backend. `Lossless` is how to say so.
        Z0 = 50.0
        capf(w) = fill((1 - im*w*1000.0e-15*Z0)/(1 + im*w*1000.0e-15*Z0), 1, 1)
        mk2(f, noise) = Circuit(
            Any[:p1 => Port(1; termination = nothing), :r1 => Resistor(Z0),
                :cc => Capacitor(100.0e-15),
                :x => ScatteringParameters(f; nports = 1, grounded = true,
                    noise = noise)],
            Any[((:p1,1), (:r1,1), (:cc,1)), ((:cc,2), (:x,1)),
                ((:r1,2), (:p1,2), Ground)])
        mk(noise) = Circuit(
            Any[:p1 => Port(1; termination = nothing), :r1 => Resistor(Z0),
                :cc => Capacitor(100.0e-15),
                :jj => JosephsonJunction(1000.0e-12),
                :cj => ScatteringParameters(capf; nports = 1, grounded = true,
                    noise = noise)],
            Any[((:p1,1), (:r1,1), (:cc,1)),
                ((:cc,2), (:jj,1), (:cj,1)),
                ((:jj,2), (:r1,2), (:p1,2), Ground)])
        for (noise, haschannels) in ((Passive(), true), (Lossless(), false))
            psc = JosephsonCircuits.compile(mk(noise))
            ssys = JosephsonCircuits.scatteringstampsystem(psc, 1;
                auxoffset = 0, Ntotal = 64, scale = 1.0)
            @test isnothing(JosephsonCircuits.planscatteringnoise(ssys)) !=
                haschannels
        end
        # the block really is unitary, so declaring it changes no answer
        wsl = 2*pi*collect(range(4.5e9, 5.0e9, length = 4))
        wpl = (2*pi*2*4.75e9,)
        srcl = [(mode=(1,), port=1, current=0.3e-6)]
        sol(noise) = hbsolve(wsl, wpl, srcl, (4,), (8,), mk(noise);
            returnSnoise = true).linearized
        a, b = sol(Passive()), sol(Lossless())
        @test a.S == b.S
        @test isapprox(a.QE, b.QE; atol = 1e-12)
        @test isapprox(a.CM, b.CM; atol = 1e-12)

        # where the data can be checked, the assertion is held to it
        @test_throws ArgumentError ScatteringParameters(ComplexF64[0 0.5; 0.5 0];
            noise = Lossless())
        @test_throws ArgumentError ScatteringParameters(
            (2*pi*[1e9, 2e9], reshape(ComplexF64[0.9, 0.9], 1, 1, :));
            nports = 1, grounded = true, noise = Lossless())
        # and unitary stored data is accepted
        @test ScatteringParameters(ComplexF64[0 1; 1 0];
            noise = Lossless()).nports == 2

        # a callable cannot be held to the declaration when it is built, so
        # it is looked at when a sweep is asked for: sampling can show that a
        # block dissipates, which is the direction that matters here
        lossyf(w) = fill(0.8*(1 - im*w*1e-12*Z0)/(1 + im*w*1e-12*Z0), 1, 1)
        wl = 2*pi*collect(range(4e9, 6e9, length = 8))
        @test_throws ArgumentError hblinsolve(wl,
            mk2(lossyf, Lossless()))
        # the same block without the declaration is fine, and dissipative
        out = hblinsolve(wl, mk2(lossyf, Passive()); keyedarrays = false,
            returnSnoise = true, returnCM = true)
        @test maximum(abs, out.Snoise) > 0
        @test all(x -> isapprox(x, 1.0; atol = 1e-12), out.CM)
        # and a truly lossless callable passes the same check
        clean = hblinsolve(wl, mk2(capf, Lossless()); keyedarrays = false,
            returnSnoise = true, returnCM = true)
        @test all(isfinite, clean.S)
        @test all(x -> isapprox(x, 1.0; atol = 1e-12), clean.CM)
    end

    @testset "a pumped circuit with a dissipative block" begin
        # the noise source of each mode is scaled by that mode's frequency,
        # which a single mode circuit cannot show. A parametric amplifier
        # with a lossy non-reciprocal block in its line has several, at
        # frequencies which differ by the pump.
        c = Circuit(
            Any[:p1 => Port(1; termination = nothing), :r1 => Resistor(50.0),
                :x => ScatteringParameters(ComplexF64[0 0; 0.8 0]; grounded = false),
                :c1 => Capacitor(100.0e-15),
                :jj => JosephsonJunction(1.0e-6),
                :cj => Capacitor(100.0e-15)],
            Any[((:p1,1), (:r1,1), (:x,1,1)),
                ((:x,2,1), (:c1,1)),
                ((:c1,2), (:jj,1), (:cj,1)),
                ((:x,1,2), (:x,2,2), (:jj,2), (:cj,2), (:r1,2), (:p1,2),
                 Ground)])
        sol = hbsolve(2*pi*[4.5e9], (2*pi*5.0e9,),
            [(mode=(1,), port=1, current=1.0e-6)], (4,), (8,), c)
        # one per mode: the commutation relations are signed by the mode
        # frequency, so the conjugated modes come back to minus one
        @test all(x -> isapprox(abs(x), 1.0; atol = 1e-12), sol.linearized.CM)
        @test all(x -> 0 <= x <= 1 + 1e-12, sol.linearized.QE ./
            sol.linearized.QEideal)
    end

    @testset "the noise covariance factorization" begin
        # any factor of the covariance describes the same noise, but the host
        # and the kernel must agree on which one, or their noise scattering
        # matrices differ row by row while their sums over rows agree. Both
        # run psdcholesky!, so this is the one thing both depend on.
        function factor(V)
            n = size(V, 1)
            L = Complex{Float64}[V[p,c] for c in 1:n for p in 1:n]
            JosephsonCircuits.psdcholesky!(L, 0, n)
            return reshape(L, n, n)
        end
        for V in (
                # positive definite
                ComplexF64[2 1-im; 1+im 3],
                # singular: the covariance of a series element, whose noise
                # leaves by one port only
                ComplexF64[1 -1; -1 1],
                # a zero row and column, which is what a port that reflects
                # perfectly gives
                ComplexF64[0 0; 0 1],
                # identically zero: a lossless block
                ComplexF64[0 0; 0 0],
                # three ports, rank two
                ComplexF64[2 1 1; 1 1 0; 1 0 1])
            L = factor(V)
            @test isapprox(L*L', V; atol = 1e-12)
            @test istril(L)
        end
        # rounding can make a lossless block's covariance slightly indefinite,
        # which must not become a negative square root
        L = factor(ComplexF64[-1e-18 0; 0 -2e-18])
        @test all(iszero, L)
    end

    @testset "the scattering block noise channels are formed on the backend" begin
        # The channels are a contraction of the auxiliary port current rows
        # of the adjoint solution, which is where that solution already is on
        # a backend, so they are formed there rather than bringing it back.
        # The kernels run on CPU() unchanged, so this covers their arithmetic
        # on every machine; what it cannot cover is the device solver itself.
        Cval = 1000.0e-15
        ftab = 2*pi*collect(range(0.5e9, 30e9, length = 200))
        tabvals = zeros(ComplexF64, 2, 2, length(ftab))
        for (k, w) in enumerate(ftab)
            z = 1/(im*w*Cval*50.0)
            tabvals[1,1,k] = tabvals[2,2,k] = 0.9*z/(z+2)
            tabvals[1,2,k] = tabvals[2,1,k] = 0.9*2/(z+2)
        end
        for (lbl, blk) in (
                # a constant matrix is stored as a one point table, so it and
                # the tabulated block share the interpolating kernel
                ("constant", ScatteringParameters(ComplexF64[0.3im 0; 0.8 0.2];
                    grounded = false)),
                ("tabulated", ScatteringParameters((ftab, tabvals); nports = 2,
                    grounded = false, extrapolation = :constant)),
                ("entry callable", ScatteringParameters(
                    (p, q, w) -> p == q ? complex(0.2) : complex(0.6);
                    nports = 2, grounded = false, form = :entry)))
            circuit = Circuit(
                Any[:p1 => Port(1; termination = nothing), :r1 => Resistor(50.0), :x => blk,
                    :l1 => Inductor(1000.0e-12), :p2 => Port(2; termination = nothing),
                    :r2 => Resistor(50.0)],
                Any[((:p1,1), (:r1,1), (:x,1,1)),
                    ((:x,2,1), (:l1,1), (:p2,1), (:r2,1)),
                    ((:x,1,2), (:x,2,2), (:l1,2), (:r1,2), (:p1,2), (:r2,2),
                     (:p2,2), Ground)])
            psc = JosephsonCircuits.compile(circuit)
            cg = JosephsonCircuits.calccircuitgraph(psc)
            sf = JosephsonCircuits.truncfreqs(
                JosephsonCircuits.calcfreqsdft((0,)); dc = true, odd = false,
                even = true, maxintermodorder = Inf)
            wsweep = 2*pi*[5.0e9]
            d = JosephsonCircuits.hblinsolve(wsweep, psc, cg,
                Dict{Symbol,Number}(), sf; debuglsys = true)
            lsys = d.lsys
            ssys = lsys.scattering
            noiseplan = JosephsonCircuits.planscatteringnoise(ssys)
            @test !isnothing(noiseplan)
            @test JosephsonCircuits.candeviceevaluate(ssys)

            # the adjoint solution, which is the transposed one
            A = copy(lsys.Asparse)
            wmodes = wsweep[1] .+ d.wpumpmodes
            JosephsonCircuits.assemblesystemmatrix!(A, lsys, wmodes)
            phiadj = Matrix(transpose(Matrix(A))) \ Matrix(d.bnm)
            nrhs = size(phiadj, 2)
            Nmodes = d.Nmodes

            # the host contraction
            nrows = noiseplan.Nchannels*Nmodes
            want = zeros(ComplexF64, nrows, nrhs)
            JosephsonCircuits.scatteringnoisewaves!(want, noiseplan, ssys,
                phiadj, wmodes, 0)

            # the same thing through the kernels
            backend = JosephsonCircuits.CPU()
            bp = JosephsonCircuits.plandeviceblocknoise(ssys, noiseplan,
                Nmodes, backend)
            dp = JosephsonCircuits.plandeviceproviders(ssys, 1, backend,
                d.wpumpmodes, ssys.scale)
            @test !isnothing(dp)
            got = zeros(ComplexF64, nrows, nrhs)
            if isnothing(dp.funcs)
                JosephsonCircuits.blocknoisefactorkernel!(backend, 64)(
                    bp.factors, bp.blockindex, bp.factoroff, dp.nports,
                    dp.freqoff, dp.nfreq, dp.freqs, dp.valoff, dp.vals,
                    dp.conjsym, dp.extrapcode, wmodes, Nmodes, bp.nentries;
                    ndrange = bp.nentries*Nmodes)
            else
                JosephsonCircuits.blocknoiseentryfactorkernel!(backend, 64)(
                    bp.factors, bp.blockindex, bp.factoroff, dp.nports,
                    dp.funcs, dp.conjsym, wmodes, Nmodes, bp.nentries;
                    ndrange = bp.nentries*Nmodes)
            end
            JosephsonCircuits.blocknoisecontractkernel!(backend, 64)(got,
                phiadj, bp.factors, bp.blockindex, bp.factoroff, bp.auxbase,
                dp.nports, bp.channelentry, bp.channellocal, wmodes, Nmodes,
                bp.nchannels, 0, nrhs;
                ndrange = bp.nchannels*Nmodes*nrhs)
            @test isapprox(got, want; rtol = 1e-12)
            @test maximum(abs, want) > 0     # or this compares two zeros

            # A sweep on a backend forms the outputs of the frequencies of a
            # batch on several workers at once, and the covariance factors
            # are written and read at one frequency, so a worker cannot share
            # them. Everything which says where each block is, it can.
            other = JosephsonCircuits.withfactors(bp)
            @test other.factors !== bp.factors
            @test length(other.factors) == length(bp.factors)
            for f in (:blockindex, :factoroff, :auxbase, :channelentry,
                      :channellocal)
                @test getfield(other, f) === getfield(bp, f)
            end
            @test (other.nentries, other.nchannels, other.nmodes) ==
                (bp.nentries, bp.nchannels, bp.nmodes)
        end
    end

    @testset "a dissipative block with the Native negative frequency rule" begin
        # ConjugateSymmetry evaluates a block at |w| and conjugates below
        # zero; Native takes the data as given at the signed frequency. The
        # noise path has to follow whichever rule the block declares, because
        # its covariance is I - S S' at the frequency the stamp used. On data
        # which is conjugate symmetric of its own accord the two rules
        # describe the same block, so they must agree exactly, noise
        # included.
        Z0 = 50.0
        symS(C, a) = w -> fill(a*(1 - im*w*C*Z0)/(1 + im*w*C*Z0), 1, 1)
        oneport(b) = Circuit(
            Any[:p1 => Port(1; termination = nothing), :r1 => Resistor(Z0),
                :cc => Capacitor(100.0e-15),
                :l1 => Inductor(1000.0e-12), :x => b],
            Any[((:p1,1), (:r1,1), (:cc,1)),
                ((:cc,2), (:l1,1), (:x,1)),
                ((:l1,2), (:r1,2), (:p1,2), Ground)])
        wsn = 2*pi*collect(range(4.0e9, 6.0e9, length = 6))
        wpn = (2*pi*2*5.0e9,)
        srcn = [(mode=(1,), port=1, current=0.3e-6)]
        for a in (1.0, 0.85)   # lossless, then dissipative
            f = symS(1000.0e-15, a)
            native = hbsolve(wsn, wpn, srcn, (4,), (8,),
                oneport(ScatteringParameters(f; nports = 1, grounded = true,
                    negative_frequency = Native())); returnSnoise = true)
            conj = hbsolve(wsn, wpn, srcn, (4,), (8,),
                oneport(ScatteringParameters(f; nports = 1, grounded = true));
                returnSnoise = true)
            @test native.linearized.S == conj.linearized.S
            @test native.linearized.Snoise == conj.linearized.Snoise
            @test native.linearized.QE == conj.linearized.QE
            @test native.linearized.CM == conj.linearized.CM
            @test all(x -> isapprox(abs(x), 1.0; atol = 1e-12),
                native.linearized.CM)
            if a < 1
                @test maximum(abs, native.linearized.Snoise) > 0
            end
        end
        # data which is not conjugate symmetric, so the two rules describe
        # different blocks and only Native is meaningful: a complex baseband
        # resonance, whose covariance still closes the commutation relations
        asym(w0, kappa, a) = w -> fill(complex(a/(1 + im*(w - w0)/kappa)), 1, 1)
        for a in (1.0, 0.85)
            out = hbsolve(wsn, wpn, srcn, (4,), (8,),
                oneport(ScatteringParameters(asym(2*pi*5.2e9, 2*pi*1e9, a);
                    nports = 1, grounded = true,
                    negative_frequency = Native())); returnSnoise = true)
            @test all(x -> isapprox(abs(x), 1.0; atol = 1e-12),
                out.linearized.CM)
        end
    end

    @testset "blocks with more than two ports" begin
        # The covariance of an n port block is n^2 entries built from 2n^3
        # interpolations, and its factorization is a triangular sweep; at one
        # and two ports neither has much to get wrong. These have a
        # covariance which is full rank and not diagonal, one which is rank
        # deficient in two directions at once, and a non-reciprocal three
        # port.
        M4 = ComplexF64[1 1 0 0; 0 1 1 0; 0 0 1 1im; 1im 0 0 1]
        smax = maximum(JosephsonCircuits.LinearAlgebra.svdvals(M4))
        fourport(S) = Circuit(
            Any[:p1 => Port(1; termination = nothing), :r1 => Resistor(50.0), :p2 => Port(2; termination = nothing),
                :r2 => Resistor(50.0), :l1 => Inductor(1000.0e-12),
                :x => ScatteringParameters(S; grounded = true),
                :ct => Capacitor(300.0e-15), :rt => Resistor(75.0)],
            Any[((:p1,1), (:r1,1), (:x,1)),
                ((:x,2), (:l1,1), (:p2,1), (:r2,1)),
                ((:x,3), (:ct,1)),
                ((:x,4), (:rt,1)),
                ((:l1,2), (:ct,2), (:rt,2), (:r1,2), (:p1,2), (:r2,2),
                 (:p2,2), Ground)])
        # a circulator with its third port terminated is an isolator
        threeport(S) = Circuit(
            Any[:p1 => Port(1; termination = nothing), :r1 => Resistor(50.0), :p2 => Port(2; termination = nothing),
                :r2 => Resistor(50.0), :l1 => Inductor(1000.0e-12),
                :x => ScatteringParameters(S; grounded = true),
                :rt => Resistor(50.0)],
            Any[((:p1,1), (:r1,1), (:x,1)),
                ((:x,2), (:l1,1), (:p2,1), (:r2,1)),
                ((:x,3), (:rt,1)),
                ((:l1,2), (:rt,2), (:r1,2), (:p1,2), (:r2,2), (:p2,2),
                 Ground)])
        wmp = 2*pi*collect(range(4.0e9, 6.0e9, length = 8))
        circulator(a) = ComplexF64[0 0 a; a 0 0; 0 a 0]
        for (lbl, circuit, nchan) in (
                # every eigenvalue of I - S S' positive, and off diagonal
                ("four port, full rank", fourport(0.75/smax*M4), 4),
                # two of them zero, so the factorization takes its zero pivot
                # branch twice at four ports
                ("four port, rank deficient", fourport(1.0/smax*M4), 4),
                # non-reciprocal and lossy at three ports
                ("three port circulator", threeport(circulator(0.7)), 3),
                # unitary, so it is provably lossless and carries no channels
                ("three port lossless", threeport(circulator(1.0)), 0))
            out = hblinsolve(wmp, circuit; keyedarrays = false,
                returnSnoise = true, returnCM = true)
            # one lumped noise port (the termination) plus the block channels
            @test size(out.Snoise, 1) == 1 + nchan
            @test all(x -> isapprox(x, 1.0; atol = 1e-12), out.CM)
            @test maximum(abs, out.S) > 0
        end

        # the same channels through the kernels, which is what a backend
        # runs; four ports is where the triangular indexing of the
        # factorization and the covariance both have something to get wrong
        for S in (0.75/smax*M4, circulator(0.7))
            n = size(S, 1)
            circuit = n == 4 ? fourport(S) : threeport(S)
            psc = JosephsonCircuits.compile(circuit)
            cg = JosephsonCircuits.calccircuitgraph(psc)
            sf = JosephsonCircuits.truncfreqs(
                JosephsonCircuits.calcfreqsdft((0,)); dc = true, odd = false,
                even = true, maxintermodorder = Inf)
            wsweep = 2*pi*[5.0e9]
            d = JosephsonCircuits.hblinsolve(wsweep, psc, cg,
                Dict{Symbol,Number}(), sf; debuglsys = true)
            lsys = d.lsys
            ssys = lsys.scattering
            noiseplan = JosephsonCircuits.planscatteringnoise(ssys)
            @test noiseplan.Nchannels == n
            A = copy(lsys.Asparse)
            wmodes = wsweep[1] .+ d.wpumpmodes
            JosephsonCircuits.assemblesystemmatrix!(A, lsys, wmodes)
            phiadj = Matrix(transpose(Matrix(A))) \ Matrix(d.bnm)
            nrhs = size(phiadj, 2)
            Nmodes = d.Nmodes
            want = zeros(ComplexF64, noiseplan.Nchannels*Nmodes, nrhs)
            JosephsonCircuits.scatteringnoisewaves!(want, noiseplan, ssys,
                phiadj, wmodes, 0)
            backend = JosephsonCircuits.CPU()
            bp = JosephsonCircuits.plandeviceblocknoise(ssys, noiseplan,
                Nmodes, backend)
            dp = JosephsonCircuits.plandeviceproviders(ssys, 1, backend,
                d.wpumpmodes, ssys.scale)
            got = zeros(ComplexF64, size(want)...)
            JosephsonCircuits.blocknoisefactorkernel!(backend, 64)(
                bp.factors, bp.blockindex, bp.factoroff, dp.nports,
                dp.freqoff, dp.nfreq, dp.freqs, dp.valoff, dp.vals,
                dp.conjsym, dp.extrapcode, wmodes, Nmodes, bp.nentries;
                ndrange = bp.nentries*Nmodes)
            JosephsonCircuits.blocknoisecontractkernel!(backend, 64)(got,
                phiadj, bp.factors, bp.blockindex, bp.factoroff, bp.auxbase,
                dp.nports, bp.channelentry, bp.channellocal, wmodes, Nmodes,
                bp.nchannels, 0, nrhs;
                ndrange = bp.nchannels*Nmodes*nrhs)
            @test isapprox(got, want; rtol = 1e-12)
            @test maximum(abs, want) > 0
        end
    end

    @testset "sensitivities of a circuit containing a scattering block" begin
        # The hybrid rows of a block make the system non-symmetric, but the
        # adjoint the contraction needs is the transposed system, which is
        # what hblinsolve solves for every circuit. Nothing about the
        # sensitivity of a lumped component changes because a block is
        # elsewhere in the circuit, so this must agree with the equivalent
        # lumped circuit and with finite differences of the block circuit's
        # own scattering parameters.
        Z0 = 50.0
        blockjpa(Cc) = Circuit(
            Any[:p1 => Port(1; termination = nothing), :r1 => Resistor(Z0), :cc => Capacitor(Cc),
                :jj => JosephsonJunction(1000.0e-12),
                :c2 => ScatteringParameters(capS(1000.0e-15, Z0); nports = 1,
                    grounded = true)],
            Any[((:p1,1), (:r1,1), (:cc,1)),
                ((:cc,2), (:jj,1), (:c2,1)),
                ((:jj,2), (:r1,2), (:p1,2), Ground)])
        lumpedjpa(Cc) = Circuit(
            Any[:p1 => Port(1; termination = nothing), :r1 => Resistor(Z0), :cc => Capacitor(Cc),
                :jj => JosephsonJunction(1000.0e-12),
                :c2 => Capacitor(1000.0e-15)],
            Any[((:p1,1), (:r1,1), (:cc,1)),
                ((:cc,2), (:jj,1), (:c2,1)),
                ((:jj,2), (:c2,2), (:r1,2), (:p1,2), Ground)])
        wsens = 2*pi*collect(4.6e9:0.01e9:4.9e9)
        wpsens = (2*pi*4.75001e9,)
        srcsens = [(mode=(1,), port=1, current=0.00565e-6)]
        Cc0 = 100.0e-15
        # a capacitor and the junction: the junction goes through the
        # operating point's own residual term, which is a different path
        sens(c, m) = hbsolve(wsens, wpsens, srcsens, (8,), (16,), c;
            keyedarrays = false, ftol = 1e-13,
            sensitivitynames = ["cc", "jj"],
            returnSsensitivity = true, sensitivityoperatingpoint = true,
            sensitivitymode = m).linearized.Ssensitivity

        fwd = sens(blockjpa(Cc0), :forward)
        @test maximum(abs, fwd) > 1        # or the comparisons are vacuous
        # the two contraction orders
        @test isapprox(fwd, sens(blockjpa(Cc0), :reverse); rtol = 1e-10,
            norm = v -> maximum(abs, v))
        # the equivalent lumped circuit
        @test isapprox(fwd, sens(lumpedjpa(Cc0), :forward); rtol = 1e-10,
            norm = v -> maximum(abs, v))
        # central finite differences of a relative perturbation
        S(c) = hbsolve(wsens, wpsens, srcsens, (8,), (16,), c;
            keyedarrays = false, ftol = 1e-13).linearized.S
        h = 1e-6
        fd = (S(blockjpa(Cc0*(1+h))) .- S(blockjpa(Cc0*(1-h))))./(2*h)
        @test isapprox(fwd[:,:,1,:], fd; rtol = 1e-6,
            norm = v -> maximum(abs, v))
        @test maximum(abs, fwd[:,:,2,:]) > 1   # the junction column is live

        # An operating point is a host object whichever backend solved for
        # it: what differentiates it evaluates transforms one direction at a
        # time, so it needs a working system rather than a snapshot, and a
        # system on a backend would put a transfer in every caller where no
        # test without a device could see it. A solve on a backend gets a
        # host twin instead, so the contraction below is the host code these
        # tests run.
        op = hbsolve(wsens, wpsens, srcsens, (8,), (16,), blockjpa(Cc0);
            keyedarrays = false, sensitivitynames = ["cc"],
            returnSsensitivity = true,
            sensitivityoperatingpoint = true).nonlinear.operatingpoint
        @test op.sys.phimatrix isa Array
        @test op.sys.phitd isa Array
        @test op.sys.sintd isa Array
        @test op.jacobian isa
            JosephsonCircuits.SparseArrays.SparseMatrixCSC{Float64,Int}
    end

    @testset "a block whose data the backend cannot evaluate" begin
        # The noise channels of a dissipative block are formed on the backend
        # through the same providers its stamps are, so a provider the
        # backend cannot evaluate costs the whole adjoint solution coming
        # back at every frequency. That is a performance cliff with no other
        # symptom, so hblinsolve warns about it; this is the condition it
        # warns on.
        Z0 = 50.0
        matrixform = ScatteringParameters(w -> fill(complex(0.5), 1, 1);
            nports = 1, grounded = true)
        entryform = ScatteringParameters((p, q, w) -> complex(0.5); nports = 1,
            grounded = true, form = :entry)
        for (blk, ok) in ((matrixform, false), (entryform, true))
            circuit = Circuit(
                Any[:p1 => Port(1; termination = nothing), :r1 => Resistor(Z0),
                    :cc => Capacitor(100.0e-15), :x => blk],
                Any[((:p1,1), (:r1,1), (:cc,1)), ((:cc,2), (:x,1)),
                    ((:r1,2), (:p1,2), Ground)])
            psc = JosephsonCircuits.compile(circuit)
            cg = JosephsonCircuits.calccircuitgraph(psc)
            sf = JosephsonCircuits.truncfreqs(
                JosephsonCircuits.calcfreqsdft((0,)); dc = true, odd = false,
                even = true, maxintermodorder = Inf)
            d = JosephsonCircuits.hblinsolve(2*pi*[5e9], psc, cg,
                Dict{Symbol,Number}(), sf; debuglsys = true)
            @test JosephsonCircuits.candeviceevaluate(d.lsys.scattering) == ok
            # either way the block is dissipative and carries noise channels
            @test !isnothing(
                JosephsonCircuits.planscatteringnoise(d.lsys.scattering))
            # and either way the answer is the same, on the host
            out = hblinsolve(2*pi*[5e9], circuit; keyedarrays = false,
                returnCM = true)
            @test isapprox(out.CM[1,1], 1.0; atol = 1e-12)
        end
    end

    @testset "the temperature of a dissipative element" begin
        # A channel at temperature T carries coth(hbar*w/(2*k*T)) times its
        # vacuum noise. That is a statement about the state of the channel,
        # so it belongs where the noise power is asked for and not in
        # `Snoise`, which is a scattering matrix. The commutation relations
        # are a statement about the transformation, so they stay at one at
        # every temperature, which is what makes them a cross check.
        Z0 = 50.0
        w0 = 2*pi*5.0e9
        wst = [w0]
        # a resistor and a dissipative block, so both kinds of channel are
        # here and can be told apart. `nothing` means the component states no
        # temperature and takes the one the analysis is run at.
        mk(rT, bT) = Circuit(
            Any[:p1 => Port(1; termination = nothing), :r1 => Resistor(Z0),
                :cc => Capacitor(100.0e-15), :l1 => Inductor(1000.0e-12),
                :x => ScatteringParameters(fill(complex(0.5), 1, 1);
                    grounded = true,
                    noise = isnothing(bT) ? Passive() : ThermalEquilibrium(bT)),
                :cl => Capacitor(200.0e-15),
                :rt => Resistor(2.0e4; temperature = rT)],
            Any[((:p1,1), (:r1,1), (:cc,1)),
                ((:cc,2), (:l1,1), (:x,1), (:cl,1)),
                ((:cl,2), (:rt,1)),
                ((:l1,2), (:rt,2), (:r1,2), (:p1,2), Ground)])
        run(c; kw...) = hblinsolve(wst, c; keyedarrays = false,
            returnSnoise = true, returnQE = true, returnCM = true, kw...)
        plain = mk(nothing, nothing)
        cold = run(plain)
        @test isapprox(cold.CM[1,1], 1.0; atol = 1e-12)

        # every channel warm: the scattering parameters do not move, the
        # commutation relations do not move, and the quantum efficiency
        # falls by exactly the occupation of the channels
        for T in (0.1, 1.0, 4.0)
            warm = run(plain; temperature = T)
            f = JosephsonCircuits.thermaloccupation(w0, T)
            @test warm.Snoise == cold.Snoise
            @test isapprox(warm.CM[1,1], 1.0; atol = 1e-12)
            @test warm.QE[1,1,1] < cold.QE[1,1,1]
            # the quantum efficiency is the signal over the signal plus the
            # noise, and only the noise term carries the occupation
            sig = sum(abs2, cold.S[1,:,1])
            nz = sum(abs2, cold.Snoise[:,1,1])
            @test isapprox(warm.QE[1,1,1],
                abs2(cold.S[1,1,1])/(sig + f*nz); rtol = 1e-10)
        end
        # zero is the default and changes nothing
        @test run(plain; temperature = 0.0).CM == cold.CM

        # a component states its own temperature and takes that instead;
        # the quantum efficiency is what moves, so it is what is compared
        qe(o) = o.QE[1,1,1]
        onlyr = run(mk(1.0, nothing))
        onlyb = run(mk(nothing, 1.0))
        both = run(plain; temperature = 1.0)
        @test qe(onlyr) < qe(cold)
        @test qe(onlyb) < qe(cold)
        @test qe(both) < qe(onlyb) < qe(onlyr) < qe(cold)
        # stating both is the same as the analysis default being that
        @test isapprox(qe(run(mk(1.0, 1.0))), qe(both); rtol = 1e-12)
        # and a stated temperature overrides the default downwards too
        @test isapprox(qe(run(mk(0.0, 0.0); temperature = 1.0)), qe(cold);
            rtol = 1e-12)
        # ThermalEquilibrium(0) is Passive
        @test isapprox(qe(run(mk(nothing, 0.0))), qe(cold); rtol = 1e-12)
        # and every one of them leaves the commutation relations at one
        for o in (onlyr, onlyb, both, run(mk(1.0, 1.0)))
            @test isapprox(o.CM[1,1], 1.0; atol = 1e-12)
        end

        # a netlist of tuples states no temperatures, so everything in one
        # takes the analysis default; that format is unchanged
        @variables Rp Cc Lj Cj
        legacy = [("P1","1","0",1), ("R1","1","0",Rp), ("C1","1","2",Cc),
                  ("Lj1","2","0",Lj), ("C2","2","0",Cj), ("R2","2","0",2.0e4)]
        defs = Dict(Lj => 1000.0e-12, Cc => 100.0e-15, Cj => 1000.0e-15,
            Rp => 50.0)
        lcold = hblinsolve(wst, legacy, defs; keyedarrays = false,
            returnCM = true, returnQE = true)
        lwarm = hblinsolve(wst, legacy, defs; keyedarrays = false,
            returnCM = true, returnQE = true, temperature = 1.0)
        @test isapprox(lcold.CM[1,1], 1.0; atol = 1e-12)
        @test isapprox(lwarm.CM[1,1], 1.0; atol = 1e-12)
        @test lwarm.QE[1,1,1] < lcold.QE[1,1,1]
        # and a parsed netlist of tuples carries no temperatures at all
        @test isempty(JosephsonCircuits.compile(
            legacy).componenttemperatures)
    end

    @testset "the added noise covariance" begin
        # `Cnoise` is the `Y` of the Gaussian channel whose `X` is `S`: the
        # circuit takes an input covariance to `S sigma S' + Cnoise`. It is
        # where the temperature enters, `Snoise` supplying the transformation
        # and the occupation of each channel the state.
        Z0 = 50.0
        w0 = 2*pi*5.0e9
        wsc = [w0]
        mkc(bT) = Circuit(
            Any[:p1 => Port(1; termination = nothing), :r1 => Resistor(Z0),
                :x => ScatteringParameters(ComplexF64[0.2 0.6; 0.6 0.2];
                    grounded = false,
                    noise = isnothing(bT) ? Passive() : ThermalEquilibrium(bT)),
                :l1 => Inductor(1000.0e-12), :p2 => Port(2; termination = nothing),
                :r2 => Resistor(Z0), :cl => Capacitor(200.0e-15),
                :rt => Resistor(2.0e4)],
            Any[((:p1,1), (:r1,1), (:x,1,1)),
                ((:x,2,1), (:l1,1), (:p2,1), (:r2,1), (:cl,1)),
                ((:cl,2), (:rt,1)),
                ((:x,1,2), (:x,2,2), (:l1,2), (:r1,2), (:p1,2), (:r2,2),
                 (:p2,2), (:rt,2), Ground)])
        run(; kw...) = hblinsolve(wsc, mkc(nothing); keyedarrays = false,
            returnSnoise = true, returnCnoise = true, returnQE = true,
            returnCM = true, kw...)
        for T in (0.0, 1.0, 4.0)
            o = run(temperature = T)
            C = o.Cnoise[:,:,1]
            Sn = o.Snoise[:,:,1]
            f = JosephsonCircuits.thermaloccupation(w0, T)
            # the definition
            want = [sum(f*Sn[c,i]*conj(Sn[c,j]) for c in axes(Sn,1))
                for i in axes(Sn,2), j in axes(Sn,2)]
            @test isapprox(C, want; rtol = 1e-12)
            # a covariance is Hermitian and positive semidefinite
            @test isapprox(C, C'; atol = 1e-12*maximum(abs, C))
            @test minimum(real.(eigvals(Hermitian(Matrix(C))))) > -1e-10
            # and its diagonal is the noise term the quantum efficiency uses
            for i in axes(C, 1)
                @test isapprox(real(C[i,i]), f*sum(abs2, Sn[:,i]);
                    rtol = 1e-12)
            end
            # with one mode at a positive frequency the commutation relations
            # are the signal plus that diagonal, at zero temperature
            if iszero(T)
                for i in axes(C, 1)
                    @test isapprox(o.CM[i,1],
                        sum(abs2, o.S[i,:,1]) + real(C[i,i]); rtol = 1e-10)
                end
            end
        end
        # the covariance scales with the occupation while Snoise does not
        cold = run()
        warm = run(temperature = 1.0)
        @test warm.Snoise == cold.Snoise
        @test isapprox(warm.Cnoise,
            JosephsonCircuits.thermaloccupation(w0, 1.0)*cold.Cnoise;
            rtol = 1e-12)
        # a block states its own temperature and only its channels warm
        stated = hblinsolve(wsc, mkc(4.0); keyedarrays = false,
            returnCnoise = true)
        @test real(stated.Cnoise[1,1,1]) > real(cold.Cnoise[1,1,1])

        # and it costs nothing unless it is asked for
        @test isempty(hblinsolve(wsc, mkc(nothing);
            keyedarrays = false).Cnoise)
    end

    @testset "unsupported cases still error clearly" begin
        # arbitrary noise covariance permits active blocks: not supported
        active = ScatteringParameters([0.0 2.0; 2.0 0.0];
            noise = NoiseCovariance([1.0 0.0; 0.0 1.0]), grounded = true)
        c = Circuit([:a => active, :r => Resistor(50.0), :p => Port(1; termination = nothing)],
            [((:p, 1), (:r, 1), (:a, 1)), ((:a, 2), Ground),
             ((:p, 2), (:r, 2), Ground)])
        @test_throws ComponentNotSupportedError compile(c)
        # a scattering block has no scalar value to perturb, so a
        # sensitivity with respect to one is rejected; sensitivities with
        # respect to the lumped components of a circuit which contains a
        # block are supported, and tested above
        C2 = 1000.0e-15
        stamped = Circuit(
            [:p1 => Port(1; termination = nothing), :r1 => Resistor(50.0),
             :cc => Capacitor(100.0e-15),
             :jj => JosephsonJunction(1000.0e-12),
             :c2 => ScatteringParameters(capS(C2, 50.0); nports = 1,
                 grounded = true)],
            [((:p1, 1), (:r1, 1), (:cc, 1)),
             ((:cc, 2), (:jj, 1), (:c2, 1)),
             ((:jj, 2), (:r1, 2), (:p1, 2), Ground)],
        )
        @test_throws ArgumentError hblinsolve(2*pi*[5e9], stamped;
            sensitivitynames = ["c2/port1"], returnSsensitivity = true)
        # gaussian channels remain unsupported
        η = 0.5
        cg = Circuit([:ch => GaussianChannel(sqrt(η)*Matrix(1.0I, 2, 2),
                (1-η)/2*Matrix(1.0I, 2, 2); nmodes = 1, grounded = true)],
            [((:ch, 1), Ground)])
        @test_throws ComponentNotSupportedError compile(cg)
    end

    @testset "one block definition used as several instances" begin
        # Elaboration deduplicates definitions by object identity so that a
        # definition can be shared, and the stamp system used to group blocks
        # by that same identity, which merged two instances of one shared
        # definition into a single block. It failed with a complaint that a
        # port appeared twice.
        capS(C) = w -> fill((1 - im*w*C*50.0)/(1 + im*w*C*50.0), 1, 1)
        mk(shared) = begin
            b1 = ScatteringParameters(capS(1e-12); nports = 1)
            b2 = shared ? b1 : ScatteringParameters(capS(1e-12); nports = 1)
            Circuit(
                Any[:p1 => Port(1; termination = nothing),
                    :r1 => Resistor(50.0), :cc => Capacitor(100e-15),
                    :b1 => b1, :l => Inductor(1e-9), :b2 => b2,
                    :jj => JosephsonJunction(1e-9), :gnd => Ground()],
                Any[[(:p1,1),(:r1,1),(:cc,1)], [(:cc,2),(:b1,1),(:l,1)],
                    [(:l,2),(:b2,1),(:jj,1)],
                    [(:p1,2),(:r1,2),(:jj,2),(:gnd,1)]])
        end
        # two instances, at their own nodes, whichever way they were built
        for shared in (true, false)
            cc = JosephsonCircuits.compile(mk(shared))
            @test length(cc.scatteringblocks) == 2
            @test [b.path for b in cc.scatteringblocks] == ["b1", "b2"]
            sys = JosephsonCircuits.scatteringstampsystem(cc, 2;
                auxoffset = 0, Ntotal = 1000, scale = 1.0)
            @test length(sys.blocks) == 2
            @test allunique([b.signalnodes for b in sys.blocks])
        end
        # and sharing the definition is the same circuit as not sharing it
        ws = 2*pi*(4.0:0.1:6.0)*1e9
        @test Array(hblinsolve(ws, mk(true)).S) ==
            Array(hblinsolve(ws, mk(false)).S)
    end

    @testset "stamp system from the compiled blocks" begin
        # A compiled block is one instance with its own terminal map, so the
        # stamp system can be built from those directly rather than from the
        # per port entries the parsed representation is limited to. The two
        # must agree in every field.
        JC = JosephsonCircuits
        one = ScatteringParameters(
            w -> fill((1 - im*w*1e-12*50.0)/(1 + im*w*1e-12*50.0), 1, 1);
            nports = 1)
        two = ScatteringParameters([0.0 1.0; 1.0 0.0]; grounded = false)
        c = Circuit(
            Any[:p1 => Port(1; termination = nothing), :r1 => Resistor(50.0),
                :b1 => one, :b2 => one, :t => two, :l => Inductor(1e-9),
                :jj => JosephsonJunction(1e-9), :gnd => Ground()],
            Any[[(:p1,1),(:r1,1),(:b1,1),(:t,1,1)], [(:t,2,1),(:l,1)],
                [(:l,2),(:b2,1),(:jj,1)],
                [(:p1,2),(:r1,2),(:jj,2),(:t,1,2),(:t,2,2),(:gnd,1)]])
        cc = JC.compile(c)
        for Nmodes in (1, 2, 4)
            a = JC.scatteringstampsystem(cc, Nmodes; auxoffset = 7,
                Ntotal = 200, scale = 1.5)
            b = JC.scatteringstampsystem(cc.scatteringblocks, Nmodes;
                auxoffset = 7, Ntotal = 200, scale = 1.5)
            for f in fieldnames(typeof(a))
                f === :blocks && continue
                @test getfield(a, f) == getfield(b, f)
            end
            @test [x.signalnodes for x in a.blocks] ==
                [x.signalnodes for x in b.blocks]
            @test [x.refnodes for x in a.blocks] ==
                [x.refnodes for x in b.blocks]
            @test [x.auxbase for x in a.blocks] ==
                [x.auxbase for x in b.blocks]
            @test [x.name for x in a.blocks] == [x.name for x in b.blocks]
        end
    end

end

# The zero frequency pencil, which is what makes a block visible at direct
# current. `evaluatehybrid!` writes `i = 0` there; the descriptor reads the
# block's own relation instead, and the classification of what it leaves
# undetermined falls out of the rank of `C(0)` rather than out of a list of
# block types.
@testset verbose=true "scattering blocks at zero frequency" begin
    JC = JosephsonCircuits
    mk(S; n = 2) = ScatteringParameters(S; nports = n, grounded = true,
        noise = Lossless())
    descr(lab, blk) = JC.dcblockdescriptor(JC.StampedScatteringBlock(
        blk, collect(1:blk.nports), zeros(Int, blk.nports), 0, lab))

    @testset "a resistive block has a determined current" begin
        # a 100 ohm series resistor between 50 ohm ports has S(0) with
        # S11 = Z/(Z+2*Z0) = 1/2 and S21 = 2*Z0/(Z+2*Z0) = 1/2
        d = descr("R", mk(w -> JC.ABCDtoS(JC.ABCD_seriesZ(100.0 + 0im))))
        r = sqrt(50.0)
        @test d.B0 ≈ [0.5 -0.5; -0.5 0.5] ./ r
        @test d.C0 ≈ [1.5 0.5; 0.5 1.5] .* r
        @test d.freecurrents == 0        # I + S(0) is invertible
    end

    @testset "an open block reproduces the row it replaces" begin
        # S(0) = I gives B0 = 0 and C0 = 2*sqrt(R), so the row is i = 0,
        # which is what the zero frequency special case hard codes. An open
        # circuit is the one block for which that special case was right.
        d = descr("open", mk(w -> Matrix{Complex{Float64}}(I, 2, 2)))
        @test all(iszero, d.B0)
        @test d.C0 ≈ 2*sqrt(50.0)*Matrix(I, 2, 2)
        @test d.freecurrents == 0
    end

    @testset "a short constrains the voltage and frees the current" begin
        # S(0) = -1 gives C0 = 0: the row is B0*V = 0, so the port voltage
        # is pinned and the current is not determined by it at all
        d = descr("short", mk(JC.S_short!(ones(Complex{Float64}, 1, 1)); n = 1))
        @test all(iszero, d.C0)
        @test d.B0 ≈ fill(2/sqrt(50.0), 1, 1)
        @test d.freecurrents == 1
    end

    @testset "an ideal through leaves one current direction free" begin
        # S(0) = [0 1; 1 0]: I + S(0) is singular, so one combination of the
        # port currents is undetermined while the voltages are tied together
        d = descr("through", mk(w -> JC.ABCDtoS(JC.ABCD_seriesZ(0.0 + 0im))))
        @test rank(d.C0) == 1
        @test d.freecurrents == 1
        @test d.B0 ≈ [1.0 -1.0; -1.0 1.0] ./ sqrt(50.0)
    end

    @testset "a block with no real zero frequency limit is refused" begin
        # a limit which exists but cannot be evaluated at zero
        @test_throws ArgumentError descr("cap",
            mk(w -> JC.ABCDtoS(JC.ABCD_seriesZ(1/(im*w*1e-12)))))
        # and one which is genuinely complex there
        @test_throws ArgumentError descr("reactive",
            mk(w -> Complex{Float64}[0 im; im 0]))
    end
end
