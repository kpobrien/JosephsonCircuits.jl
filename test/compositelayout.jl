using JosephsonCircuits
using LinearAlgebra
using Random
using Test

@testset verbose=true "canonical state layout" begin
    JC = JosephsonCircuits

    @testset "layout of a state with one zero frequency mode" begin
        # 80 nodes, 77 modes, mode 1 is the zero frequency one and the only
        # self conjugate mode: the shape of the two tone hard case
        isdc = [i == 1 for i in 1:77]
        ml = JC.ModeLayout(isdc, 80*77)
        L = JC.compositelayout(ml, isdc)

        @test ml.rdim == 12240
        @test L.ndc == 80                  # one per node
        @test JC.canonicaldim(L) == ml.rdim
        @test L.nvdc == 0
        @test JC.isinternal(L)
        @test JC.nwindow(L) == 80

        # each node contributes 1 + 2*76 = 153 internal entries with the
        # zero frequency one first, so the window names exactly those
        @test L.dcpos == [k*153 + 1 for k in 0:79]
        @test JC.windowindices(L) == L.dcpos
        @test [JC.windowindex(L, k) for k in 1:80] == L.dcpos
    end

    @testset "a self conjugate mode which is not zero frequency stays alternating" begin
        # mode 40 is self conjugate (a Nyquist mode) but not direct current
        isreal = [i == 1 || i == 40 for i in 1:77]
        isdc   = [i == 1 for i in 1:77]
        ml = JC.ModeLayout(isreal, 80*77)
        L = JC.compositelayout(ml, isdc)
        @test L.ndc == 80                  # only the zero frequency mode
        @test JC.canonicaldim(L) == ml.rdim
        @test JC.isinternal(L)
    end

    @testset "the mode tuples pick out the zero frequency mode" begin
        modes = [(0,0), (1,0), (0,1), (1,1)]
        ml = JC.ModeLayout([all(iszero, m) for m in modes], 5*4)
        @test JC.compositelayout(ml, modes).ndc == 5
    end

    @testset "the internal state is the first block, and the voltages follow" begin
        isdc = [i == 1 for i in 1:9]
        ml = JC.ModeLayout(isdc, 7*9)
        L = JC.compositelayout(ml, isdc; nvdc = 2)
        @test JC.canonicaldim(L) == ml.rdim + 2
        @test JC.voltagerange(L) == (ml.rdim + 1):(ml.rdim + 2)
        @test JC.windowindices(L) == vcat(L.dcpos, ml.rdim + 1, ml.rdim + 2)
        @test !JC.isinternal(L)
        Random.seed!(11)
        r = randn(ml.rdim)
        u = zeros(ml.rdim + 2); u[end-1:end] .= (3.0, 4.0)
        JC.gathercanonical!(u, r, L)
        @test u[1:ml.rdim] == r            # a copy, so bit exact
        @test u[end-1:end] == [3.0, 4.0]   # and the voltages are left alone
        back = similar(r)
        JC.scattercanonical!(back, u, L)
        @test back == r
        @test JC.internalpart(u, L) == r
    end

    @testset "rejected layouts" begin
        isreal = [i == 1 for i in 1:9]
        ml = JC.ModeLayout(isreal, 7*9)
        # a zero frequency mode has no imaginary part, so it must be self
        # conjugate; marking a conjugate pair as direct current is a bug
        @test_throws ArgumentError JC.compositelayout(ml, [i == 2 for i in 1:9])
        @test_throws DimensionMismatch JC.compositelayout(ml, trues(3))
        @test_throws ArgumentError JC.compositelayout(ml, isreal; nvdc = -1)

        L = JC.compositelayout(ml, isreal)
        @test_throws DimensionMismatch JC.gathercanonical!(zeros(3), zeros(ml.rdim), L)
        @test_throws DimensionMismatch JC.gathercanonical!(zeros(ml.rdim), zeros(3), L)
        @test_throws DimensionMismatch JC.scattercanonical!(zeros(ml.rdim), zeros(3), L)
    end

    @testset "residual and product are the same operator in the new basis" begin
        # a small two tone circuit, taken to a point that is not the origin
        circuit = Circuit(
            [:p1 => Port(1), :cc => Capacitor(100e-15),
             :jj => JosephsonJunction(1000e-12), :cj => Capacitor(1000e-15),
             :gnd => Ground()],
            [[(:p1, 1), (:cc, 1)],
             [(:cc, 2), (:jj, 1), (:cj, 1)],
             [(:p1, 2), (:jj, 2), (:cj, 2), (:gnd, 1)]])
        srcs = [(mode = (1,), port = 1, current = 0.5e-6)]
        # three iterations is a truncated solve on purpose: what this needs
        # is a point away from the origin, not the solution, so the solver
        # reports that it did not converge and that is asserted here rather
        # than printed
        s = @test_logs((:warn,), match_mode=:any,
            JC.hbnlsolve((2*pi*4.75e9,), (4,), srcs, circuit;
                dc = true, odd = true, even = true,
                returnoperatingpoint = true, iterations = 3))
        sys = s.operatingpoint.sys
        ml = s.operatingpoint.modelayout
        L = JC.compositelayout(ml, s.frequencies.modes)
        @test JC.isinternal(L)
        @test L.ndc > 0

        Random.seed!(3)
        x = randn(L.rdim); v = randn(L.rdim)
        u = similar(x); vc = similar(v)
        JC.gathercanonical!(u, x, L); JC.gathercanonical!(vc, v, L)
        work = JC.CanonicalWork(L, x)

        # the residual in the new basis is the old one, gathered
        Fi = zeros(L.rdim)
        JC.setpoint!(sys, x); JC.residual!(Fi, sys)
        want = similar(Fi); JC.gathercanonical!(want, Fi, L)
        got = zeros(L.rdim); uu = copy(u)
        JC.canonicalresidual(
            (F, J, xx) -> (JC.setpoint!(sys, xx);
                isnothing(F) || JC.residual!(F, sys); nothing), work)(got, nothing, uu)
        @test got == want
        @test uu == u                      # the point is handed back unchanged

        # and so is the Jacobian vector product
        Ji = zeros(L.rdim)
        JC.setpoint!(sys, x); JC.jacobianvectorproduct!(Ji, sys, v)
        wantj = similar(Ji); JC.gathercanonical!(wantj, Ji, L)
        gotj = zeros(L.rdim)
        JC.setpoint!(sys, x)
        JC.canonicaljvp((o, vv) -> JC.jacobianvectorproduct!(o, sys, vv),
            work)(gotj, vc)
        @test gotj == wantj
    end

    @testset "the layout is used exactly when it is needed" begin
        # A circuit which injects no direct current does not build the
        # explicit block, and one which does builds it. That is the whole
        # condition now: there is no flag, and no second path to compare
        # against, so what is asserted is that each case solves and that the
        # block appears only in the second.
        circuit = Circuit(
            [:p1 => Port(1), :cc => Capacitor(100e-15),
             :jj => JosephsonJunction(1000e-12), :cj => Capacitor(1000e-15),
             :gnd => Ground()],
            [[(:p1, 1), (:cc, 1)],
             [(:cc, 2), (:jj, 1), (:cj, 1)],
             [(:p1, 2), (:jj, 2), (:cj, 2), (:gnd, 1)]])
        kw = (; dc = true, odd = true, even = true, ftol = 1e-12)

        # alternating current only: the voltages are the zero they sit at
        ac = JC.hbnlsolve((2*pi*4.75e9,), (8,),
            [(mode = (1,), port = 1, current = 1.2e-6)], circuit; kw...)
        @test ac.solverinfo.converged
        @test all(iszero, ac.dcnodevoltage)

        # with direct current in the drive the block is built and reports
        dc = JC.hbnlsolve((2*pi*4.75e9,), (8,),
            [(mode = (1,), port = 1, current = 1.2e-6),
             (mode = (0,), port = 1, current = 1.0e-7)], circuit;
            kw..., rtol = 1e-12)
        @test dc.solverinfo.converged
        @test !isnothing(dc.dcnodevoltage)
    end

    @testset "the methods which can carry the block, and the one which cannot" begin
        circuit = Circuit(
            [:p1 => Port(1; Z0 = 50.0), :cc => Capacitor(100e-15),
             :jj => JosephsonJunction(1000e-12), :gnd => Ground()],
            [[(:p1, 1), (:cc, 1)], [(:cc, 2), (:jj, 1)],
             [(:p1, 2), (:jj, 2), (:gnd, 1)]])
        srcs = [(mode = (0,), port = 1, current = 1.0e-7)]

        # the two which solve the real system take the explicit block
        for m in (NewtonKrylov(), Newton())
            s = JC.hbnlsolve((2*pi*4.75e9,), (4,), srcs, circuit;
                dc = true, odd = true, method = m, rtol = 1e-12)
            @test s.solverinfo.converged
            @test !isnothing(s.dcnodevoltage)
        end

        # `:quasinewton` solves the complex holomorphic system, which has no
        # place for the real direct current unknowns, and says so
        @test_throws ArgumentError JC.hbnlsolve((2*pi*4.75e9,), (4,), srcs,
            circuit; dc = true, odd = true, method = QuasiNewton())
    end

    @testset "the assembled Jacobian survives a growing internal pattern" begin
        # The canonical Jacobian's pattern is the internal pattern under a
        # permutation plus the direct current block, and none of it moves.
        # It used to be rebuilt from the *values* at each point, and sparse
        # addition prunes exact zeros, so the pattern it produced was the
        # one the starting point happened to have. A junction driven hard
        # enough to develop harmonics fills in mode coupling entries which
        # were zero at the origin, the pattern grew, and the solve stopped
        # with the internal error which guarded against exactly that.
        #
        # The plan is built from the pattern alone, so it holds everywhere.
        # The check is that `:newton` reaches the same point `:newtonkrylov`
        # does, on a drive strong enough for the fill in to happen. The
        # drive is a few times the critical current, so the junction is
        # strongly nonlinear, but not so far past it that the operating
        # point is one of several and the two paths pick different ones.
        circuit = Circuit(
            [:p1 => Port(1; Z0 = 50.0), :cc => Capacitor(100e-15),
             :jj => JosephsonJunction(1000e-12), :cj => Capacitor(1000e-15)],
            [[(:p1,1),(:cc,1)], [(:cc,2),(:jj,1),(:cj,1)],
             [(:p1,2),(:jj,2),(:cj,2), Ground]])
        srcs = [(mode = (1,), port = 1, current = 1.0e-6),
                (mode = (0,), port = 1, current = 1.0e-7)]
        kw = (; dc = true, odd = true, even = true, keyedarrays = false,
              rtol = 1e-12)
        a = JC.hbnlsolve((2*pi*4.75e9,), (8,), srcs, circuit;
            kw..., method = Newton())
        b = JC.hbnlsolve((2*pi*4.75e9,), (8,), srcs, circuit;
            kw..., method = NewtonKrylov())
        @test a.solverinfo.converged
        @test b.solverinfo.converged
        # the scattering parameters and the direct current operating point
        # are the physical content and agree to roundoff
        @test maximum(abs, a.S .- b.S) < 1e-10
        @test isapprox(a.dcnodevoltage, b.dcnodevoltage; rtol = 1e-8)
        # the node fluxes may differ by a whole number of flux quanta, which
        # is the additive static gauge: a junction phase is defined modulo
        # 2*pi and the two methods are free to land on different branches
        turns = (a.nodeflux .- b.nodeflux) ./ (2*pi)
        @test all(x -> isapprox(x, round(real(x)); atol = 1e-8), turns)
    end

    @testset "the assembled Jacobian is the matrix free one" begin
        # A circuit with an explicit direct current block: a scattering
        # block which is a resistor, driven at zero frequency. The assembled
        # canonical Jacobian has to reproduce the matrix free product for
        # every unit vector, which is a sharper check than differencing the
        # residual and is what lets the direct solve methods use the same
        # formulation as the Krylov one.
        blk = ScatteringParameters(
            w -> JC.ABCDtoS(JC.ABCD_seriesZ(100.0 + 0im));
            nports = 2, grounded = true, noise = Lossless())
        c = Circuit(
            [:p1 => Port(1; Z0 = 1.0e9), :b => blk, :c1 => Capacitor(1e-12)],
            [[(:p1,1),(:b,1),(:c1,1)], [(:b,2), Ground],
             [(:p1,2),(:c1,2), Ground]])
        srcs = [(mode = (0,), port = 1, current = 1.0e-6)]
        d = JC.hbnlsolve((2*pi*5e9,), (1,), srcs, c, Dict{Any,Any}();
            dc = true, odd = true, keyedarrays = false, returnsystem = true)

        sys, ml, Nmodes = d.sys, d.modelayout, d.Nmodes
        Nnodes = length(d.dcplan.componentof) + 1
        tr = JC.transportrows(d.dcplan, d.bnmsource, Nmodes)
        L = JC.compositelayout(ml, d.frequencies.modes;
            nvdc = JC.nvoltages(tr))
        psc = JC.compile(c)
        Naux = ml.dim - d.Nnodal
        nauxsc = JC.countscatteringports(psc)*Nmodes
        ssys = JC.scatteringstampsystem(psc.scatteringblocks, Nmodes;
            auxoffset = d.Nnodal + Naux - nauxsc,
            Ntotal = d.Nnodal + Naux, scale = d.Lmean)
        br = JC.dcblockrows(ssys.blocks, d.dcplan.componentof, Nmodes,
            d.dcplan.modeindex, Nnodes - 1, d.Lmean)
        work = JC.CanonicalWork(L, zeros(ml.rdim); transport = tr,
            blockrows = br, nnodaldc = Nnodes - 1)
        @test L.nvdc > 0                    # the block is really explicit

        n = JC.canonicaldim(L)
        JC.setpoint!(sys, zeros(ml.rdim))
        jvp = JC.canonicaljvp(
            (o, v) -> JC.jacobianvectorproduct!(o, sys, v), work)
        free = zeros(n, n); col = zeros(n); out = zeros(n)
        for k in 1:n
            fill!(col, 0.0); col[k] = 1.0
            jvp(out, col); free[:,k] .= out
        end
        # the coupling the preconditioner subtracts is the Jacobian's own:
        # the columns of the direct current coordinates on the nodal zero
        # frequency rows, which is the resistor current the voltages drive
        # and the two terminals each block current enters
        H = JC.dccoupling(work)
        idx = JC.dcsubsystemindices(work)
        nodal = L.dcpos[1:size(H, 1)]
        @test free[nodal, idx] ≈ Matrix(H)
        # and no other zero frequency row sees the direct current
        # coordinates from outside the subsystem: the block is triangular
        below = setdiff(JC.windowindices(L), nodal, idx)
        @test all(iszero, free[below, idx])

        Jint = copy(d.Jr)
        JC.jacobian!(Jint, sys)
        assembled = JC.canonicaljacobian(Jint, work)
        @test size(assembled) == (n, n)
        @test Matrix(assembled) == free      # exactly, not to a tolerance
    end
end
