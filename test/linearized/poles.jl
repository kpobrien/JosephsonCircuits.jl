using Test, JosephsonCircuits, LinearAlgebra, SparseArrays, Random

# a backend whose eigensolve gives the host's multipliers and no vectors,
# as a device's gives none where its left vectors fail their residuals
struct RefusingBackend <: JosephsonCircuits.KernelAbstractions.Backend end
function JosephsonCircuits.mapspectrum!(M::Matrix{Float64}, ilo::Int, ihi::Int, ::RefusingBackend)
    values, _ = JosephsonCircuits.mapspectrum!(copy(M), ilo, ihi, JosephsonCircuits.KernelAbstractions.CPU())
    return values, ks -> nothing
end

@testset "HB poles" begin
    JC = JosephsonCircuits
    # Compare sets without depending on the eigensolver's ordering.
    distance(a, b) = maximum(x -> minimum(abs.(x .- b)), a; init = 0.0)
    closepoles(a, b; tol = 1e-9) = length(a) == length(b) &&
        max(distance(a, b), distance(b, a)) < tol

    @testset "RLC poles, units, ports, and the sparse inverse" begin
        for w0 in (1.0, 2pi*8e9)
            C, L, R = 0.5/w0, 2/w0, 2.0
            c = Circuit([(:r, 1, 0, Resistor(R)),
                (:c, 1, 0, Capacitor(C)), (:l, 1, 0, Inductor(L))])
            expected = w0 .* [-0.5 + sqrt(3)/2*im, -0.5 - sqrt(3)/2*im]
            p = hbstability(c)
            @test closepoles(p.poles/w0, expected/w0)
            @test p.converged
            @test p.infinite == 0
            @test all(iszero, p.edgeweights)
            @test size(p.nodevoltage) == (1, 1, 2)
            @test size(p.junctionflux) == (1, 0, 2)
            for scale in (0.3w0, 4w0)
                q = hbstability(compile(c); frequencyscale = scale)
                @test closepoles(q.poles/w0, expected/w0)
            end
            Random.seed!(421)
            q = hbstability(c; method = ShiftInvert(w0*(0.1 + im); nev = 2))
            @test q.converged
            @test closepoles(q.poles/w0, expected/w0)
            @test all(==(1), q.shiftindices)
            @test only(q.searches).converged == 2
            # Port termination is physical damping; reference impedance
            # alone does not add damping to an unterminated port.
            port = Circuit([(:p, 1, 0, Port(1; Z0 = R)),
                (:c, 1, 0, Capacitor(C)), (:l, 1, 0, Inductor(L))])
            @test closepoles(hbstability(port).poles/w0, expected/w0)

            sys = JC.hbpolesystem(compile(c), Dict(), nothing, (0,), w0)
            A, B = JC.polecompanion(sys)
            z = 0.15 + 0.3im
            n = size(sys.Q0, 1)
            op = JC.PoleShiftInvert(JC.kluordered(sys.Q0 + z*sys.Q1 + z^2*sys.Q2), sys.Q1, sys.Q2, z,
                zeros(ComplexF64, n), zeros(ComplexF64, n))
            x = randn(ComplexF64, size(A, 1))
            y = similar(x)
            mul!(y, op, x)
            @test y ≈ (A-z*B)\(B*x) rtol = 1e-12
            mul!(y, op, 2x)
            @test y ≈ (A-z*B)\(B*(2x)) rtol = 1e-12
        end
    end

    @testset "DC coordinates and algebraic modes" begin
        rc = Circuit([(:r, 1, 0, Resistor(2.0)), (:c, 1, 0, Capacitor(0.5))])
        p = hbstability(rc)
        @test p.poles ≈ [-1.0]
        @test p.infinite == 1
        @test p.converged
        q = hbstability(rc; method = ShiftInvert(0.1im; nev = 1))
        @test q.converged
        @test q.poles ≈ [-1.0]
        # A capacitor with no discharge path has a physical zero mode.
        # It must survive the removal of the arbitrary constant flux.
        cap = Circuit([(:c, 1, 0, Capacitor(0.5))])
        p = hbstability(cap)
        @test p.poles == [0.0]
        @test norm(p.nodevoltage) > 0
        @test p.converged
        q = hbstability(cap; method = ShiftInvert(0.1im; nev = 1))
        @test q.converged
        @test abs(only(q.poles)) < 1e-12
        # An algebraic voltage divider with one capacitor has one pole,
        # not two slow poles invented by regularizing the capacitance.
        divider = Circuit([(:r1, 1, 0, Resistor(2.0)),
            (:r2, 1, 2, Resistor(3.0)), (:c, 2, 0, Capacitor(0.5))])
        p = hbstability(divider)
        @test p.poles ≈ [-0.4]
        @test p.infinite == 3
        @test p.converged
        resistor = Circuit([(:r, 1, 0, Resistor(2.0))])
        p = hbstability(resistor)
        @test isempty(p.poles)
        @test p.infinite == 2
        Random.seed!(18)
        q = hbstability(resistor; method = ShiftInvert(0.1im))
        @test isempty(q.poles)
        # Negative conductance is an allowed active, time-local element.
        active = Circuit([(:r, 1, 0, Resistor(-2.0)), (:c, 1, 0, Capacitor(0.5))])
        @test hbstability(active).poles ≈ [1.0]
    end

    @testset "floating inductive island and mutual inductance" begin
        # Common voltage decays with rate G/C; the differential mode is
        # an RLC oscillator with stiffness 2/L. No artificial flux zero.
        for reversed in (false, true)
            a, b = reversed ? (2, 1) : (1, 2)
            c = Circuit([(:r1, 1, 0, Resistor(4.0)),
                (:r2, 2, 0, Resistor(4.0)),
                (:c1, 1, 0, Capacitor(1.0)), (:c2, 2, 0, Capacitor(1.0)),
                (:l, a, b, Inductor(2.0))])
            expected = ComplexF64[-0.25, -0.125+sqrt(1-0.125^2)*im,
                -0.125-sqrt(1-0.125^2)*im]
            p = hbstability(c)
            @test closepoles(p.poles, expected)
            @test p.converged
        end
        c = Circuit([(:c1, 1, 0, Capacitor(1.0)),
            (:c2, 2, 0, Capacitor(1.0)), (:l1, 1, 0, Inductor(2.0)),
            (:l2, 2, 0, Inductor(2.0)), (:k, :l1, :l2, MutualInductor(0.3))])
        p = hbstability(c)
        expected = im .* [1/sqrt(2*1.3), -1/sqrt(2*1.3),
            1/sqrt(2*0.7), -1/sqrt(2*0.7)]
        @test closepoles(p.poles, expected)
        @test p.converged
        @test p.infinite == 4
        perfect = Circuit([(:c1, 1, 0, Capacitor(1.0)),
            (:c2, 2, 0, Capacitor(1.0)), (:l1, 1, 0, Inductor(2.0)),
            (:l2, 2, 0, Inductor(2.0)), (:k, :l1, :l2, MutualInductor(1.0))])
        p = hbstability(perfect)
        @test closepoles(p.poles, [0.5im, -0.5im])
        @test p.infinite == 6
        @test p.converged
    end

    @testset "pumped junction against an independent monodromy integration" begin
        c = Circuit([(:p, 1, 0, Port(1; Z0 = 20.0)),
            (:j, 1, 0, JosephsonJunction(1.0)), (:c, 1, 0, Capacitor(1.0))])
        omega = 0.8
        pump = hbnlsolve((omega,), (12,),
            [(mode = (1,), port = 1, current = 0.04*JC.phi0)], c;
            method = Newton(), atol = 1e-12, keyedarrays = false)
        @test pump.solverinfo.converged
        # Integrate d/dt [dphi; dv] = [0 1; -cos(phi*) -G] [dphi; dv]
        # directly, using a hand-written RK4 and Fourier reconstruction.
        phase(t) = sum(2real(pump.nodeflux[k]*cis(only(mode)*omega*t))
            for (k, mode) in enumerate(pump.modes))
        function monodromy(nsteps)
            dt = 2pi/omega/nsteps
            U = Matrix{Float64}(I, 2, 2)
            rhs(t, X) = [0.0 1.0; -cos(phase(t)) -0.05]*X
            for k in 0:nsteps-1
                t = k*dt
                a = rhs(t, U)
                b = rhs(t+dt/2, U+dt/2*a)
                d = rhs(t+dt/2, U+dt/2*b)
                e = rhs(t+dt, U+dt*d)
                U += dt/6*(a+2b+2d+e)
            end
            return eigvals(U)
        end
        expected = monodromy(2048)
        @test closepoles(monodromy(4096), expected; tol = 1e-10)
        reference = nothing
        for H in (4, 6)
            p = hbstability(c; nonlinear = pump, Nmodulationharmonics = (H,), method = DenseSpectrum())
            @test p.converged
            # Select the physical representatives in the pump's central
            # strip, rather than folding all inaccurate boundary aliases.
            indices = findall(s -> -omega/2 < imag(s) <= omega/2, p.poles)
            @test length(indices) == 2
            multipliers = exp.(p.poles[indices]*2pi/omega)
            @test closepoles(multipliers, expected; tol = 2e-8)
            @test maximum(p.edgeweights[indices]) < 1e-4
            # each pole moved to the harmonic which carries most of it, the
            # harmonics the truncation shapes aside, is one of the two
            # physical modes: its aliases fall together
            moved = [s + im*n*omega for (s, n) in zip(p.poles, p.harmonics) if abs(n) <= H - 2]
            @test length(moved) == 2*(2H - 3)
            @test all(z -> min(abs(z - moved[1]), abs(z - conj(moved[1]))) < 1e-8, moved)
            if !isnothing(reference)
                @test closepoles(p.poles[indices], reference; tol = 1e-8)
            end
            reference = p.poles[indices]
            # the node's two exponents sum to -G/C, Liouville's formula for
            # its one capacitance and conductance, so a pair decays at -1/40
            @test all(s -> isapprox(real(s), -1/40; atol = 1e-8), p.poles[indices])
        end
        # The period map: each mode once, its multipliers against the
        # hand-written monodromy's, falling as the fourth power of the
        # step, and the least damped mode's profile against the dense
        # spectrum's.
        mapped(steps) = hbstability(c; nonlinear = pump, Nmodulationharmonics = (4,),
            method = Monodromy(; nev = :all, steps))
        coarse, fine = mapped(64), mapped(128)
        @test fine.converged && length(fine.poles) == 2
        multipliers(r) = exp.(r.poles*2pi/omega)
        @test closepoles(multipliers(fine), expected; tol = 5e-7)
        @test distance(multipliers(fine), expected) < distance(multipliers(coarse), expected)/12
        # each rate's error, sixteen fifteenths of its change at twice the
        # steps, every mode carried there, against the change of the finer
        # map's own rates, and falling with the step as they do
        @test coarse.rateerrors ≈ (16/15) .* abs.(real.(fine.poles) .- real.(coarse.poles)) rtol = 1e-6
        @test all(fine.rateerrors .< coarse.rateerrors ./ 12)
        dense = hbstability(c; nonlinear = pump, Nmodulationharmonics = (4,), method = DenseSpectrum())
        k = argmin(abs.(dense.poles .- fine.poles[1]))
        a, b = vec(fine.junctionflux[:, :, 1]), vec(dense.junctionflux[:, :, k])
        @test abs(dot(a, b))/(norm(a)*norm(b)) ≈ 1 atol = 1e-10
        @test length(hbstability(c; nonlinear = pump, Nmodulationharmonics = (4,), method = Monodromy(nev = 1)).poles) == 1
        # an inductive divider beside the junction, whose middle node no
        # capacitor, resistor or junction touches: the map on the divider's
        # constraint has the junction's two modes alone, whose rates sum to
        # -G/C as Liouville's formula has it for the one node, within the
        # step's error
        divider = Circuit([(:p, 1, 0, Port(1; Z0 = 20.0)), (:j, 1, 0, JosephsonJunction(1.0)),
            (:c, 1, 0, Capacitor(1.0)), (:l1, 1, 2, Inductor(2.0)), (:l2, 2, 0, Inductor(3.0))])
        dividerpump = hbnlsolve((omega,), (12,), [(mode = (1,), port = 1, current = 0.04*JC.phi0)], divider;
            method = Newton(), atol = 1e-12, keyedarrays = false)
        split = hbstability(divider; nonlinear = dividerpump, method = Monodromy(nev = :all))
        @test length(split.poles) == 2 && isapprox(sum(real, split.poles), -1/20; rtol = 1e-5)
        # a junction on a node no capacitor touches, behind an inductor: the
        # orbit meets the node's constraint only to the truncation of its
        # harmonics, and the map starts the transient on it projected; its
        # two modes, split at the pump frequency, have rates summing to
        # -G/C, and the dense spectrum holds them modulo the pump frequency
        behind = Circuit([(:p, 1, 0, Port(1; Z0 = 20.0)), (:c, 1, 0, Capacitor(1.0)),
            (:l, 1, 2, Inductor(0.5)), (:j, 2, 0, JosephsonJunction(1.0))])
        behindpump = hbnlsolve((omega,), (8,), [(mode = (1,), port = 1, current = 0.02*JC.phi0)], behind;
            keyedarrays = false)
        projected = hbstability(behind; nonlinear = behindpump, method = Monodromy(nev = :all))
        behinddense = hbstability(behind; nonlinear = behindpump, method = DenseSpectrum())
        alias(z) = complex(real(z), mod(imag(z) + omega/2, omega) - omega/2)
        @test length(projected.poles) == 2 && isapprox(sum(real, projected.poles), -1/20; rtol = 1e-5)
        @test all(s -> minimum(abs.(alias.(behinddense.poles) .- alias(s))) < 1e-5, projected.poles)
        # Behind a weak port the junction is nearly lossless, its rate a
        # millionth of its frequency: over 512 steps the tangent's stages,
        # solved exactly at each step, give the dense spectrum's.
        weak = Circuit([(:p, 1, 0, Port(1; Z0 = 1e6)),
            (:j, 1, 0, JosephsonJunction(1.0)), (:c, 1, 0, Capacitor(1.0))])
        weakpump = hbnlsolve((omega,), (12,), [(mode = (1,), port = 1, current = 0.04*JC.phi0)], weak;
            method = Newton(), atol = 1e-12, keyedarrays = false)
        rate = maximum(real, hbstability(weak; nonlinear = weakpump, method = DenseSpectrum()).poles)
        exact = hbstability(weak; nonlinear = weakpump, Nmodulationharmonics = (4,), method = Monodromy(steps = 512, nev = 1))
        @test abs(real(only(exact.poles)) - rate) < 1e-6*abs(rate)
        Random.seed!(917)
        iterative = hbstability(c; nonlinear = pump, Nmodulationharmonics = (6,),
            method = ShiftInvert([0.01+0.2im, 0.01-0.2im]; nev = 2))
        @test iterative.converged
        @test length(iterative.searches) == 2
        @test distance(reference, iterative.poles) < 1e-8
        @test all(iterative.residuals .< 1e-9)

        # Small-signal scattering and the pole pencil have the same
        # homogeneous nullspace, up to a constant row scaling.
        sys = JC.hbpolesystem(compile(c), Dict(), pump, (4,), omega)
        dbg = hblinsolve([0.31], c; nonlinear = pump,
            Nmodulationharmonics = (4,), threewavemixing = true,
            fourwavemixing = true, debuglsys = true)
        A = copy(dbg.lsys.Asparse)
        JC.assemblesystemmatrix!(A, dbg.lsys, 0.31)
        z = 0.31im/omega
        Q = sys.Q0+z*sys.Q1+z^2*sys.Q2
        scale = Matrix(Q)/Matrix(A)
        @test norm(scale-Diagonal(diag(scale))) < 1e-10

        # A converged HB solution can be unstable: the inverted junction
        # at zero drive has negative stiffness and an analytic growing pole.
        inverted = hbnlsolve((omega,), (0,),
            [(mode = (0,), port = 1, current = 0.0)], c;
            dc = true, x0 = ComplexF64[pi], method = Newton(), atol = 1e-12)
        @test inverted.solverinfo.converged
        p = hbstability(c; nonlinear = inverted, Nmodulationharmonics = (0,), method = DenseSpectrum())
        @test closepoles(p.poles, ComplexF64[-0.025+sqrt(1+0.025^2), -0.025-sqrt(1+0.025^2)])
        @test maximum(real, p.poles) > 0
    end

    @testset "the period map's eigensolve against LAPACK's" begin
        # mapeigen!, the Schur form and the chosen multipliers' vectors
        # alone, against geevx, which forms every vector and condition, on
        # maps turned into every coordinate by a reflection: a complex pair
        # at the top whose partner nev = 1 leaves out, an isolated
        # multiplier beside a defective block, and a nilpotent map, whose
        # multipliers none resolves. The multipliers are geevx's, and so are
        # the chosen ones' vectors and conditions; with nev = :all each
        # multiplier is tested against its own bound and the count left
        # unresolved is geevx's, with nev = 1 the largest alone are.
        function geevx(M)
            _, wr, wi, VL, VR, _, _, _, abnrm, rconde, _ = LAPACK.geevx!('S', 'V', 'V', 'E', copy(M))
            return complex.(wr, wi), VL, VR, abnrm, rconde
        end
        reflection(m) = (w = collect(1.0:m); I - 2*w*w'/(w'*w))
        turned(D) = reflection(size(D, 1))*D*reflection(size(D, 1))
        jordan(m) = (J = diagm(vcat(1.01, fill(0.2, m - 1))); foreach(k -> J[k, k + 1] = 1.0, 2:m - 1); J)
        pair = zeros(4, 4)
        pair[1:2, 1:2] .= 0.9 .* [cos(0.7) -sin(0.7); sin(0.7) cos(0.7)]
        pair[3, 3], pair[4, 4] = 0.5, 0.1
        # the map, the multipliers nev = 1 chooses, and tests, and those
        # nev = :all leaves unresolved
        for (M, chosen, tested, unresolved) in ((turned(pair), 1, 1, 0), (jordan(24), 1, 1, 23),
                (turned(jordan(8)), 1, 1, 0), (turned(diagm(1 => ones(5))), 0, 6, 6))
            λ, VL, VR, abnrm, rconde = geevx(M)
            one, every = JC.mapeigen!(copy(M), 1, JC.CPU()), JC.mapeigen!(copy(M), typemax(Int), JC.CPU())
            @test isapprox(one.values, λ; rtol = 1e-13) && isapprox(every.values, λ; rtol = 1e-13)
            @test length(one.selected) == chosen && one.tested == tested
            @test every.unresolved == count(k -> abs(λ[k]) <= eps()*abnrm/rconde[k], eachindex(λ)) == unresolved
            @test one.unresolved == count(k -> abs(λ[k]) <= eps()*abnrm, eachindex(λ)) + tested - chosen
            wi, columns = imag.(λ), collect(eachindex(λ))
            for (j, k) in enumerate(one.selected)
                parallel(a, b) = abs(dot(a, b)) ≈ norm(a)*norm(b)
                @test parallel(JC.eigenvector(one.right, k, one.column, wi), JC.eigenvector(VR, k, columns, wi))
                @test parallel(JC.eigenvector(one.left, k, one.column, wi), JC.eigenvector(VL, k, columns, wi))
                @test one.conditions[j] ≈ rconde[k] rtol = 1e-12
            end
        end
        # where a backend gives no vectors, the host's Schur path takes
        # over on the balanced map it leaves, and gives the host's result
        for M in (turned(pair), jordan(24))
            @test JC.mapeigen!(copy(M), 1, RefusingBackend()) == JC.mapeigen!(copy(M), 1, JC.CPU())
        end
    end

    @testset "unsupported models and search input" begin
        c = Circuit([(:r, 1, 0, Resistor(2.0)), (:c, 1, 0, Capacitor(0.5))])
        @test_throws ArgumentError ShiftInvert(ComplexF64[])
        @test_throws ArgumentError ShiftInvert([NaN])
        @test_throws ArgumentError ShiftInvert(1.0; nev = 0)
        @test_throws ArgumentError ShiftInvert(1.0; nev = 3, krylovdim = 2)
        @test_throws ArgumentError ShiftInvert(1.0; nev = 3, krylovdim = 3)
        @test_throws ArgumentError DenseSpectrum(tol = 0)
        @test_throws ArgumentError Monodromy(steps = 0)
        @test_throws ArgumentError Monodromy(nev = 0)
        @test_throws ArgumentError Monodromy(nev = :some)
        # a profile of 17 harmonics in a period of 8 samples
        @test_throws ArgumentError hbstability(c; pumpfrequency = 1.0, Nmodulationharmonics = (8,), method = Monodromy(steps = 8))
        @test_throws ArgumentError hbstability(c; method = :other)
        @test_throws ArgumentError hbstability(c; factorization = QRfactorization())
        @test_throws ArgumentError DenseSpectrum(maxunknowns = 0)
        # a search for more poles than the circuit has: its one pole, and the
        # companion's infinite eigenvalues it comes upon set aside
        rc = hbstability(c; method = ShiftInvert(0.1im))
        @test rc.converged && rc.poles ≈ [-1.0] && rc.infinite > 0
        # two nodes of time constants 1 and 1e-8, the fast one growing or
        # decaying: the search comes upon all four of the companion's
        # eigenvalues, and the remote pole it does not resolve is rejected,
        # not counted among the dense spectrum's two infinite ones
        for r2 in (-1.0, 1.0)
            twin = Circuit([(:r1, 1, 0, Resistor(1.0)), (:c1, 1, 0, Capacitor(1.0)),
                (:r2, 2, 0, Resistor(r2)), (:c2, 2, 0, Capacitor(1e-8))])
            dense = hbstability(twin; method = DenseSpectrum())
            Random.seed!(123)
            searched = hbstability(twin; method = ShiftInvert(0.1im))
            @test searched.infinite <= dense.infinite
            @test all(s -> minimum(abs.(dense.poles .- s)) <= 1e-8*abs(s), searched.poles)
            @test searched.converged == (length(searched.poles) == length(dense.poles))
        end
        # an unpumped ladder of 1001 nodes exceeds the dense spectrum the
        # default takes without a pump, which names the method it can give
        ladder = Circuit(vcat([(Symbol(:c, k), k, 0, Capacitor(1.0)) for k in 1:1001],
            [(Symbol(:l, k), k, k + 1, Inductor(1.0)) for k in 1:1000], [(:r, 1, 0, Resistor(1.0))]))
        refusal = try hbstability(ladder) catch e; e end
        @test refusal isa ArgumentError && !occursin("Monodromy", sprint(showerror, refusal))
        # 1201 harmonics of one node exceed the default's thousand unknowns
        @test_throws ArgumentError hbstability(c; pumpfrequency = 1.0, Nmodulationharmonics = (600,), method = DenseSpectrum())
        @test_throws ArgumentError hbstability(c; frequencyscale = Inf)
        @test_throws ArgumentError hbstability(c; Nmodulationharmonics = (2, 2))
        for cap in (Capacitor(1-0.01im), Capacitor(FrequencyDependent(w -> 1.0)))
            bad = Circuit([(:c, 1, 0, cap), (:r, 1, 0, Resistor(2.0))])
            @test_throws ArgumentError hbstability(bad)
        end
        bad = Circuit([(:p, 1, 0, Port(1)),
            (:s, 1, ScatteringParameters(zeros(1, 1); zref = 50.0))])
        @test isempty(hbstability(bad).poles)
    end
end
