using JosephsonCircuits, LinearAlgebra, SparseArrays, Random, Test

@testset verbose=true "the mode layout and its real form" begin
    JC = JosephsonCircuits

    # the number of real slots of complex index i
    width(L, i) = Int(L.ptr[i+1] - L.ptr[i])

    # the real form of a sparse matrix built from coordinates, one real block
    # per stored entry, as an independent reference for `complex_to_real`
    function complex_to_real_ref(A, rl, cl)
        I, J, V = Int[], Int[], Float64[]
        for j in 1:size(A,2), idx in nzrange(A, j)
            i, a = rowvals(A)[idx], nonzeros(A)[idx]
            r0, wr = rl.ptr[i], width(rl, i)
            c0, wc = cl.ptr[j], width(cl, j)
            push!(I, r0); push!(J, c0); push!(V, real(a))
            wr == 2 && (push!(I, r0+1); push!(J, c0); push!(V, imag(a)))
            if wc == 2
                push!(I, r0); push!(J, c0+1); push!(V, -imag(a))
                wr == 2 && (push!(I, r0+1); push!(J, c0+1); push!(V, real(a)))
            end
        end
        return sparse(I, J, V, rl.rdim, cl.rdim)
    end

    canon(xc, L) = [width(L, i) == 1 ? Complex(real(xc[i]), 0.0) : xc[i] for i in 1:L.dim]

    CASES = (([true,false,false,false,false], [true,false,false,false,false], 8, 8, 0.1),
                   ([true,false,false,false,false], [true,false,false,false,false], 12, 7, 0.2),
                   ([true,true,false], [false,false,false], 10, 6, 0.3),
                   ([false,false], [true,true], 9, 11, 0.25),
                   ([false], [false], 20, 15, 0.2),          # no real modes at all
                   ([true], [true], 13, 5, 0.4),             # every mode real
                   (rand(Bool, 6), rand(Bool, 6), 5, 4, 0.15))

    @testset "ModeLayout" begin
        L = JC.ModeLayout([true,false,false,false,false], 10)
        @test L.rdim == 18 && count(L.isreal) == 1 && L.nmodes == 5
        @test L.ptr == [1,2,4,6,8,10,11,13,15,17,19]
        @test L.inv == [1,2,2,3,3,4,4,5,5,6,7,7,8,8,9,9,10,10]
        @test JC.ModeLayout([true,false,true,false], 8).rdim == 12
        @test JC.ModeLayout(falses(4), 8).rdim == 16
        @test JC.ModeLayout(trues(4), 8).rdim == 8
        @test_throws DimensionMismatch JC.ModeLayout([true,false], 7)
    end

    @testset "complex_to_real of a sparse matrix" begin
        for (rmask, cmask, mmul, nmul, p) in CASES
            m, n = length(rmask)*mmul, length(cmask)*nmul
            rl, cl = JC.ModeLayout(rmask, m), JC.ModeLayout(cmask, n)
            A = sprand(ComplexF64, m, n, p)
            Ar = JC.complex_to_real(A, rl, cl)
            R = complex_to_real_ref(A, rl, cl)
            @test size(Ar) == (rl.rdim, cl.rdim)
            # the values and the stored structure, a real block of every
            # stored entry kept whatever its value
            @test Ar == R
            @test SparseArrays.getcolptr(Ar) == SparseArrays.getcolptr(R)
            @test rowvals(Ar) == rowvals(R)
            # the real form of the product
            xc = canon(rand(ComplexF64, n), cl)
            @test Ar * JC.complex_to_real(xc, cl.isreal) ≈
                JC.complex_to_real(A * xc, rl.isreal)
        end
        @test JC.complex_to_real(sprand(ComplexF32, 6, 6, 0.4),
            JC.ModeLayout([true,false], 6), JC.ModeLayout([true,false], 6)) isa
            SparseMatrixCSC{Float32,Int}
        # a layout of another dimension is refused
        mask = [true,false,false,false,false]
        A = sprand(ComplexF64, 20, 15, 0.2)
        rl, cl = JC.ModeLayout(mask, 20), JC.ModeLayout(mask, 15)
        @test_throws DimensionMismatch JC.complex_to_real(A, JC.ModeLayout(mask, 25), cl)
        @test_throws DimensionMismatch JC.complex_to_real(A, rl, JC.ModeLayout(mask, 25))
    end

    disr(mask, i) = mask[(i - 1) % length(mask) + 1]
    dcanon(x, mask) = [disr(mask, i) ? Complex(real(x[i]), 0.0) : x[i] for i in eachindex(x)]

    @testset "realdim / complexdim" begin
        @test JC.realdim(10, [true,false,false,false,false]) == 18
        @test JC.realdim(8, falses(4)) == 16
        @test JC.realdim(8, trues(4)) == 8
        @test JC.complexdim(JC.realdim(30, [true,false,false]), [true,false,false]) == 30
        @test_throws DimensionMismatch JC.realdim(7, [true,false])
        @test_throws DimensionMismatch JC.complexdim(7, [true,false,false])
        @test_throws ArgumentError JC.realdim(0, Bool[])
    end

    @testset "vectors" begin
        L = [true,false,false,false,false]
        @test JC.complex_to_real(ComplexF64[1, 2+3im, 4+5im, 6+7im, 8+9im], L) == Float64[1,2,3,4,5,6,7,8,9]
        for (mask, d) in ((( [true,false,false,false,false]), 20), ([false], 7), ([true], 6),
                          ([true,true,false], 12), (rand(Bool, 8), 32))
            xr = rand(JC.realdim(d, mask))
            @test JC.complex_to_real(JC.real_to_complex(xr, mask), mask) == xr
            xc = rand(ComplexF64, d)
            @test JC.real_to_complex(JC.complex_to_real(xc, mask), mask) == dcanon(xc, mask)
            @test JC.real_to_complex!(fill(ComplexF64(NaN, NaN), d), xr, mask) == JC.real_to_complex(xr, mask)
            @test_throws DimensionMismatch JC.complex_to_real!(Vector{Float64}(undef, length(xr)+1), xc, mask)
            @test_throws DimensionMismatch JC.real_to_complex!(Vector{ComplexF64}(undef, d+length(mask)), xr, mask)
        end
    end
end

@testset "the gather and scatter kernels are inverse permutations" begin
    # the device side of the canonical layout's permuted copies, run on
    # the CPU backend here
    v = randn(Random.default_rng(), 9)
    index = [4, 1, 9, 2, 7]
    got = zeros(2, 3)
    JosephsonCircuits.gathervalues!(got, v, reshape(vcat(index, 5), 2, 3))
    @test vec(got) == v[vcat(index, 5)]
    @test_throws DimensionMismatch JosephsonCircuits.gathervalues!(
        zeros(1, 1), v, index)
    w = fill(NaN, 9)
    JosephsonCircuits.scattervalues!(w, v[index], index)
    @test w[index] == v[index] && all(isnan, w[setdiff(1:9, index)])
    @test_throws DimensionMismatch JosephsonCircuits.scattervalues!(
        w, v[index], index[1:end-1])
end

# The canonical state layout: which entries of the state a mode owns,
# and the gather and scatter between the full state and the canonical
# one the solvers iterate on.
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
        @test iszero(L.nvdc)
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
        @test iszero(L.nvdc)
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
        @test !iszero(L.nvdc)
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
        # The canonical Jacobian's pattern is the internal pattern as it is
        # plus the direct current block, and none of it moves.
        # A plan built from the values at one point would carry only the
        # entries nonzero there, since sparse addition prunes exact zeros,
        # and a junction driven hard enough to develop harmonics fills in
        # mode coupling entries which are zero at the origin. The plan is
        # built from the pattern alone, so it holds everywhere.
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
            Ntotal = d.Nnodal + Naux, scale = d.Lscale,
            modeoffsets = zeros(Nmodes))
        br = JC.dcblockrows(ssys.blocks, d.dcplan.componentof, Nmodes,
            d.dcplan.modeindex, Nnodes - 1, d.Lscale)
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
        assembled = JC.canonicaljacobian!(JC.canonicaljacobianplan(Jint, work),
            Jint)
        @test size(assembled) == (n, n)
        @test Matrix(assembled) == free      # exactly, not to a tolerance
    end
end
