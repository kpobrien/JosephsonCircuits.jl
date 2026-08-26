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
        blk = ScatteringBlock(f; nports = 1, grounded = true)
        ws = [0.0, 2*pi*5e9, -2*pi*5e9]
        B = zeros(Complex{Float64}, 1, 1, 3)
        Cm = zeros(Complex{Float64}, 1, 1, 3)
        _, _, lossy = JosephsonCircuits.evaluatehybrid!(B, Cm, blk, ws)
        @test !lossy
        # zero frequency rows are i = 0
        @test B[1,1,1] == 0 && Cm[1,1,1] == 1
        @test isapprox(B[1,1,2]/Cm[1,1,2], im*ws[2]*C; rtol = 1e-12)
        # the negative frequency rule gives the conjugate
        @test isapprox(B[1,1,3]/Cm[1,1,3], conj(B[1,1,2]/Cm[1,1,2]);
            rtol = 1e-12)
        # an ideal short: B = 2/sqrt(Z0) forces v = 0, C = 0
        short = ScatteringBlock(w -> fill(-1.0 + 0.0im, 1, 1);
            nports = 1, grounded = true)
        _, _, _ = JosephsonCircuits.evaluatehybrid!(B, Cm, short,
            [2*pi*5e9, 0.0, 2*pi*1e9])
        @test isapprox(B[1,1,1], 2/sqrt(Z0); rtol = 1e-12)
        @test Cm[1,1,1] == 0
        # a lossy attenuator is flagged
        att = ScatteringBlock([0.0 0.5; 0.5 0.0])
        B2 = zeros(Complex{Float64}, 2, 2, 1)
        C2m = zeros(Complex{Float64}, 2, 2, 1)
        _, _, lossy2 = JosephsonCircuits.evaluatehybrid!(B2, C2m, att,
            [2*pi*5e9])
        @test lossy2
    end

    # the JPA from the README with the shunt capacitor C2 expressed as a
    # one port grounded scattering block. the stamped system must reproduce
    # the lumped circuit through the full nonlinear + linearized solve.
    function jpacircuit(shunt)
        return Circuit(
            [:p1 => Port(1), :r1 => Resistor(50.0),
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
        blk = ScatteringBlock(capS(C2, 50.0); nports = 1, grounded = true)
        stamped = Circuit(
            [:p1 => Port(1), :r1 => Resistor(50.0),
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
        blk = ScatteringBlock(seriesS; nports = 2)
        lumped = jpacircuit(Capacitor(1000.0e-15))
        stamped = Circuit(
            [:p1 => Port(1), :r1 => Resistor(50.0),
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
            [:p1 => Port(1), :r1 => Resistor(50.0),
             :cc => Capacitor(100.0e-15),
             :jj => JosephsonJunction(1000.0e-12),
             :c2 => ScatteringBlock((wgrid, Sgrid); grounded = true,
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
            [:p1 => Port(1), :r1 => Resistor(50.0),
             :cc => Capacitor(100.0e-15),
             :jj => JosephsonJunction(1000.0e-12),
             :c2 => ScatteringBlock((wshort, Sshort); grounded = true)],
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
                ScatteringBlock(capS(C, 50.0); nports = 1, grounded = true)
            components = Any[:p1 => Port(1), :r1 => Resistor(50.0),
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
                :r2 => Resistor(50.0), :p2 => Port(2))
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
        stub1 = ScatteringBlock(w -> fill(exp(-2*im*w*tau), 1, 1);
            nports = 1, grounded = true, negative_frequency = Native())
        function stubcircuit(stub, extra...)
            components = Any[:p1 => Port(1), :r1 => Resistor(50.0),
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
            [:p1 => Port(1), :r1 => Resistor(50.0), :line => line,
             :r2 => Resistor(50.0), :p2 => Port(2)],
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

    @testset "hbnlsolve pump path and lossy warning" begin
        # the pump solve alone with a scattering block matches the lumped
        # circuit
        C2 = 1000.0e-15
        lumped = jpacircuit(Capacitor(C2))
        stamped = Circuit(
            [:p1 => Port(1), :r1 => Resistor(50.0),
             :cc => Capacitor(100.0e-15),
             :jj => JosephsonJunction(1000.0e-12),
             :c2 => ScatteringBlock(capS(C2, 50.0); nports = 1,
                 grounded = true)],
            [((:p1, 1), (:r1, 1), (:cc, 1)),
             ((:cc, 2), (:jj, 1), (:c2, 1)),
             ((:jj, 2), (:r1, 2), (:p1, 2), Ground)],
        )
        out1 = hbnlsolve(wp, (8,), sources, lumped)
        out2 = hbnlsolve(wp, (8,), sources, stamped)
        @test isapprox(out1.nodeflux[:], out2.nodeflux[:]; rtol = 1e-8)

        # a dissipative block warns once that its noise is not included
        lossy = Circuit(
            [:p1 => Port(1), :r1 => Resistor(50.0),
             :att => ScatteringBlock([0.0 0.5; 0.5 0.0]; grounded = true),
             :rt => Resistor(50.0)],
            [((:p1, 1), (:r1, 1), (:att, 1)),
             ((:att, 2), (:rt, 1)),
             ((:rt, 2), (:p1, 2), (:r1, 2), Ground)],
        )
        @test_logs (:warn, r"dissipative") match_mode=:any hblinsolve(
            2*pi*[5e9], lossy)
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
            blk = ScatteringBlock((ftab, tab); nports = 1, grounded = true,
                extrapolation = :constant)
            circuit = Circuit(
                Any[:p1 => Port(1), :r1 => Resistor(50.0),
                    :cc => Capacitor(100.0e-15),
                    :jj => JosephsonJunction(1000.0e-12), :c2 => blk,
                    :rl => Resistor(1.0e5)],
                Any[((:p1,1), (:r1,1), (:cc,1)),
                    ((:cc,2), (:jj,1), (:c2,1), (:rl,1)),
                    ((:jj,2), (:r1,2), (:p1,2), Ground), ((:rl,2), Ground)])
            nl = hbnlsolve(wp, (16,), sources, circuit; keyedarrays = false)
            psc = JosephsonCircuits.parsesortcircuit(circuit)
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
        cb = ScatteringBlock(capS(Cval, 50.0); nports = 1, grounded = true)
        ccircuit = Circuit(
            Any[:p1 => Port(1), :r1 => Resistor(50.0),
                :cc => Capacitor(100.0e-15),
                :jj => JosephsonJunction(1000.0e-12), :c2 => cb],
            Any[((:p1,1), (:r1,1), (:cc,1)), ((:cc,2), (:jj,1), (:c2,1)),
                ((:jj,2), (:r1,2), (:p1,2), Ground)])
        cnl = hbnlsolve(wp, (16,), sources, ccircuit; keyedarrays = false)
        cpsc = JosephsonCircuits.parsesortcircuit(ccircuit)
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
        ret = ScatteringBlock(capS(Cval, 50.0); nports = 1, grounded = true)
        inp = ScatteringBlock((d, w) -> (d[1,1] =
                (1 - im*w*Cval*50.0)/(1 + im*w*Cval*50.0); d);
            nports = 1, grounded = true, form = :inplace)
        ent = ScatteringBlock((p, q, w) ->
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
            Any[:p1 => Port(1), :r1 => Resistor(50.0),
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
        @test_throws ArgumentError ScatteringBlock([0.0 0.5; 0.5 0.0];
            form = :entry)
        @test_throws ArgumentError ScatteringBlock(
            ([1.0e9, 2.0e9], zeros(ComplexF64,1,1,2)); nports = 1,
            grounded = true, form = :inplace)
        @test_throws ArgumentError ScatteringBlock(capS(Cval, 50.0);
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
        blk = ScatteringBlock(mkentry(Cval); nports = 1, grounded = true,
            form = :entry)
        circuit = Circuit(
            Any[:p1 => Port(1), :r1 => Resistor(50.0),
                :cc => Capacitor(100.0e-15),
                :jj => JosephsonJunction(1000.0e-12), :c2 => blk,
                :rl => Resistor(1.0e5)],
            Any[((:p1,1), (:r1,1), (:cc,1)),
                ((:cc,2), (:jj,1), (:c2,1), (:rl,1)),
                ((:jj,2), (:r1,2), (:p1,2), Ground), ((:rl,2), Ground)])
        nl = hbnlsolve(wp, (16,), sources, circuit; keyedarrays = false)
        psc = JosephsonCircuits.parsesortcircuit(circuit)
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
        two = ScatteringBlock(mkasym(1.0); nports = 2, form = :entry)
        tcircuit = Circuit(
            Any[:p1 => Port(1), :r1 => Resistor(50.0), :cc => two,
                :jj => JosephsonJunction(1000.0e-12),
                :c2 => Capacitor(1000.0e-15)],
            Any[((:p1,1), (:r1,1), (:cc,1,1)),
                ((:cc,2,1), (:jj,1), (:c2,1)),
                ((:cc,1,2), (:cc,2,2), (:jj,2), (:c2,2), (:r1,2), (:p1,2),
                 Ground)])
        tnl = hbnlsolve(wp, (16,), sources, tcircuit; keyedarrays = false)
        tpsc = JosephsonCircuits.parsesortcircuit(tcircuit)
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
        heavy = ScatteringBlock((p, q, w) ->
                (1 - im*w*buf[1]*1e-12*50.0)/(1 + im*w*buf[1]*1e-12*50.0);
            nports = 1, grounded = true, form = :entry)
        hcircuit = Circuit(
            Any[:p1 => Port(1), :r1 => Resistor(50.0),
                :cc => Capacitor(100.0e-15),
                :jj => JosephsonJunction(1000.0e-12), :c2 => heavy],
            Any[((:p1,1), (:r1,1), (:cc,1)), ((:cc,2), (:jj,1), (:c2,1)),
                ((:jj,2), (:r1,2), (:p1,2), Ground)])
        hnl = hbnlsolve(wp, (16,), sources, hcircuit; keyedarrays = false)
        hpsc = JosephsonCircuits.parsesortcircuit(hcircuit)
        hcg = JosephsonCircuits.calccircuitgraph(hpsc)
        hd = JosephsonCircuits.hblinsolve(wsweep[1:3], hpsc, hcg,
            Dict{Symbol,Number}(), sf; nonlinear = hnl, debuglsys = true)
        @test !JosephsonCircuits.candeviceevaluate(hd.lsys.scattering)
    end

    @testset "a lossy non-reciprocal block: the adjoint is not a transpose" begin
        # Every other block here is reciprocal, where S equals its own
        # transpose, so the conjugated pump system and the transposed one
        # agree in more ways than they do in general and a confusion between
        # them is invisible. An isolator is not reciprocal: it transmits from
        # port one to port two and blocks the reverse, and scaling it by less
        # than one makes it dissipative as a real isolator is. With one in the
        # circuit the two constructions differ, and the assembly which the
        # device path builds for the adjoint has to be the conjugated pump
        # one.
        a = 0.7
        for (lbl, blk) in (
                ("constant", ScatteringBlock(ComplexF64[0 0; a 0];
                    nports = 2)),
                ("entry", ScatteringBlock((p, q, w) ->
                        (p == 2 && q == 1) ? complex(a) : complex(0.0);
                    nports = 2, form = :entry)))
            circuit = Circuit(
                Any[:p1 => Port(1), :r1 => Resistor(50.0), :iso => blk,
                    :jj => JosephsonJunction(1000.0e-12),
                    :cg => Capacitor(500.0e-15),
                    :p2 => Port(2), :r2 => Resistor(50.0),
                    :cl => Capacitor(200.0e-15), :rl => Resistor(2.0e4)],
                Any[((:p1,1), (:r1,1), (:iso,1,1)),
                    ((:iso,2,1), (:jj,1), (:cg,1), (:p2,1), (:r2,1),
                     (:cl,1)),
                    ((:cl,2), (:rl,1)),
                    ((:iso,1,2), (:iso,2,2), (:jj,2), (:cg,2), (:r1,2),
                     (:p1,2), (:r2,2), (:p2,2), (:rl,2), Ground)])
            nl = hbnlsolve(wp, (16,), sources, circuit; keyedarrays = false)
            psc = JosephsonCircuits.parsesortcircuit(circuit)
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
                Nnoiseports = length(d.noiseportimpedanceindices),
                Ncomponents = 0, Nnodes = d.Nnodes,
                Nfrequencies = length(wsweep))
            @test JosephsonCircuits.needsadjointsolve(arrays,
                d.noiseportimpedanceindices)
            sol = hbsolve(wsweep, wp, sources, (16,), (2,), circuit)
            @test maximum(abs, sol.linearized.S((0,),2,(0,),1,:)) > 0.1
            @test maximum(abs, sol.linearized.S((0,),1,(0,),2,:)) < 1e-10

            # the forward and the conjugated pump system, which for this
            # block are genuinely different matrices
            A = copy(lsys.Asparse)
            nz = JosephsonCircuits.SparseArrays.nnz(A)
            perm = JosephsonCircuits.cscvaluepermutation(A)
            dest = JosephsonCircuits.sweepdestinations(A,
                lsys.scattering.Aindex, false)
            JosephsonCircuits.assemblesystemmatrix!(A, lsys,
                wsweep[1] .+ d.wpumpmodes)
            forwardnz = copy(JosephsonCircuits.SparseArrays.nonzeros(A))
            JosephsonCircuits.assemblesystemmatrix!(A, lsys,
                wsweep[1] .+ d.wpumpmodes; conjugatepump = true)
            @test forwardnz != JosephsonCircuits.SparseArrays.nonzeros(A)

            for conj in (false, true)
                host = Matrix{ComplexF64}(undef, nz, length(wsweep))
                for (i, w) in enumerate(wsweep)
                    JosephsonCircuits.assemblesystemmatrix!(A, lsys,
                        w .+ d.wpumpmodes; conjugatepump = conj)
                    host[:, i] .=
                        JosephsonCircuits.SparseArrays.nonzeros(A)
                end
                plan, _, _ = JosephsonCircuits.planfrequencysweep(lsys,
                    JosephsonCircuits.CPU(); adjoint = false,
                    conjugatepump = conj)
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
                @test got[invperm(perm), :] == host
            end
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
                ("lossless", ScatteringBlock(capS(Cval, 50.0); nports = 1,
                    grounded = true)),
                ("lossy", ScatteringBlock(lossyS(Cval, 50.0, 0.85);
                    nports = 1, grounded = true)))
            circuit = Circuit(
                Any[:p1 => Port(1), :r1 => Resistor(50.0),
                    :cc => Capacitor(100.0e-15),
                    :jj => JosephsonJunction(1000.0e-12), :c2 => blk,
                    :rl => Resistor(1.0e5)],
                Any[((:p1,1), (:r1,1), (:cc,1)),
                    ((:cc,2), (:jj,1), (:c2,1), (:rl,1)),
                    ((:jj,2), (:r1,2), (:p1,2), Ground), ((:rl,2), Ground)])
            nl = hbnlsolve(wp, (16,), sources, circuit; keyedarrays = false)
            psc = JosephsonCircuits.parsesortcircuit(circuit)
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
            dest = JosephsonCircuits.sweepdestinations(A,
                lsys.scattering.Aindex, false)

            # the forward system, and the conjugated pump system which is the
            # adjoint when scattering blocks are present: their hybrid rows
            # break the similarity that otherwise lets the adjoint come from a
            # transposed solve
            for conj in (false, true)
                host = Matrix{ComplexF64}(undef, nz, length(wsweep))
                for (i, w) in enumerate(wsweep)
                    JosephsonCircuits.assemblesystemmatrix!(A, lsys,
                        w .+ d.wpumpmodes; conjugatepump = conj)
                    host[:, i] .= JosephsonCircuits.SparseArrays.nonzeros(A)
                end
                plan, _, _ = JosephsonCircuits.planfrequencysweep(lsys,
                    JosephsonCircuits.CPU(); adjoint = false,
                    conjugatepump = conj)
                got = Matrix{ComplexF64}(undef, nz, length(wsweep))
                JosephsonCircuits.assemblesweep!(got, plan, wsweep)
                st = JosephsonCircuits.plandevicescattering(lsys.scattering,
                    dest, nz, length(wsweep), JosephsonCircuits.CPU(),
                    d.Nmodes)
                JosephsonCircuits.stagescatteringstamps!(st, wsweep, 1,
                    length(wsweep), d.wpumpmodes)
                JosephsonCircuits.applyscatteringstamps!(got, st)
                # the plan stores its values in compressed sparse row order
                @test got[invperm(perm), :] == host
            end
        end
    end

    @testset "unsupported cases still error clearly" begin
        # arbitrary noise covariance permits active blocks: not supported
        active = ScatteringBlock([0.0 2.0; 2.0 0.0];
            noise = NoiseCovariance([1.0 0.0; 0.0 1.0]), grounded = true)
        c = Circuit([:a => active, :r => Resistor(50.0), :p => Port(1)],
            [((:p, 1), (:r, 1), (:a, 1)), ((:a, 2), Ground),
             ((:p, 2), (:r, 2), Ground)])
        @test_throws ComponentNotSupportedError parsesortcircuit(c)
        # sensitivities with respect to scattering blocks are rejected
        C2 = 1000.0e-15
        stamped = Circuit(
            [:p1 => Port(1), :r1 => Resistor(50.0),
             :cc => Capacitor(100.0e-15),
             :jj => JosephsonJunction(1000.0e-12),
             :c2 => ScatteringBlock(capS(C2, 50.0); nports = 1,
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
        @test_throws ComponentNotSupportedError parsesortcircuit(cg)
    end
end
