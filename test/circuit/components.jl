using JosephsonCircuits
using LinearAlgebra
using Test

@testset verbose=true "component models" begin

    @testset "matrix providers and evaluation" begin
        # constant provider with conjugate symmetry
        S = [0.0 0.8im; 0.8im 0.0]
        blk = ScatteringParameters(S)
        ws = [-2*pi*5e9, 2*pi*5e9]
        dest = zeros(Complex{Float64}, 2, 2, 2)
        JosephsonCircuits.evaluatescattering!(dest, blk, ws)
        @test dest[:,:,2] == S
        @test dest[:,:,1] == conj.(S)

        # tabulated provider: interpolation and extrapolation policies
        freqs = [1.0, 2.0, 3.0]
        vals = zeros(Complex{Float64}, 1, 1, 3)
        vals[1,1,:] .= [0.1, 0.3, 0.5]
        tblk = ScatteringParameters((freqs, vals))
        d = zeros(Complex{Float64}, 1, 1, 1)
        JosephsonCircuits.evaluatescattering!(d, tblk, [1.5])
        @test d[1,1,1] ≈ 0.2
        @test_throws ArgumentError JosephsonCircuits.evaluatescattering!(
            d, tblk, [4.0])
        tconst = ScatteringParameters((freqs, vals); extrapolation = :constant)
        JosephsonCircuits.evaluatescattering!(d, tconst, [4.0])
        @test d[1,1,1] ≈ 0.5
        tlin = ScatteringParameters((freqs, vals); extrapolation = :linear)
        JosephsonCircuits.evaluatescattering!(d, tlin, [4.0])
        @test d[1,1,1] ≈ 0.7
        # unitary knots do not prove a lossless interpolant: S = 1 and
        # S = -1 interpolate to a perfect absorber between them, and a
        # linearly extrapolated table is unbounded
        JC = JosephsonCircuits
        @test !JC.provablylossless(ScatteringParameters(([1.0, 3.0],
            reshape(ComplexF64[1, -1], 1, 1, 2))))
        @test JC.provablylossless(ScatteringParameters(([1.0, 3.0],
            reshape(ComplexF64[1, 1], 1, 1, 2))))
        @test JC.provablylossless(ScatteringParameters(([1.0, 3.0],
            reshape(ComplexF64[1, 1], 1, 1, 2)); extrapolation = :constant))
        @test !JC.provablylossless(ScatteringParameters(([1.0, 3.0],
            reshape(ComplexF64[1, 1], 1, 1, 2)); extrapolation = :linear))
        @test_throws ArgumentError ScatteringParameters(([1.0, 3.0],
            reshape(ComplexF64[1, -1], 1, 1, 2)); noise = Lossless())
        # strictly increasing frequencies required, so a repeated knot,
        # which would put a zero interval under the interpolation, is out
        @test_throws ArgumentError ScatteringParameters(([2.0, 1.0],
            vals[:,:,1:2]))
        @test_throws ArgumentError ScatteringParameters(([1.0, 1.0],
            vals[:,:,1:2]))
        @test_throws ArgumentError ScatteringParameters((freqs, vals);
            interpolation = :quartic)

        # cubic interpolation: exact at the knots, and it follows a
        # rotating phase where the chords cut across it
        fdense = collect(range(1.0, 10.0, 50))
        Sphase = zeros(Complex{Float64}, 1, 1, 50)
        for (k, w) in enumerate(fdense)
            Sphase[1,1,k] = 0.9*cis(-2*w)
        end
        pcub = ScatteringParameters((fdense, Sphase))
        plin = ScatteringParameters((fdense, Sphase); interpolation = :linear)
        JosephsonCircuits.evaluatescattering!(d, pcub, [fdense[7]])
        @test d[1,1,1] ≈ Sphase[1,1,7] atol = 1e-14
        errc = errl = 0.0
        for k in 1:49
            w = (fdense[k] + fdense[k+1])/2
            JosephsonCircuits.evaluatescattering!(d, pcub, [w])
            errc = max(errc, abs(d[1,1,1] - 0.9*cis(-2*w)))
            JosephsonCircuits.evaluatescattering!(d, plin, [w])
            errl = max(errl, abs(d[1,1,1] - 0.9*cis(-2*w)))
        end
        @test errc < errl/10
        # `:linear` extrapolation continues with the interpolant's end
        # slope; `:constant` holds the end value
        pext = ScatteringParameters((fdense, Sphase); extrapolation = :linear)
        JosephsonCircuits.evaluatescattering!(d, pext, [10.0 + 1e-3])
        endslope = pext.provider.endslopes[1,1,2]
        @test d[1,1,1] ≈ Sphase[1,1,end] + endslope*1e-3 atol = 1e-14
        @test abs(endslope - (-2*im)*0.9*cis(-20.0)) < 0.05
        pconst = ScatteringParameters((fdense, Sphase);
            extrapolation = :constant)
        JosephsonCircuits.evaluatescattering!(d, pconst, [11.0])
        @test d[1,1,1] == Sphase[1,1,end]

        # a table too short for a cubic takes the highest order it
        # determines: three samples their parabola, two their line
        f3 = [1.0, 2.0, 4.0]
        v3 = zeros(Complex{Float64}, 1, 1, 3)
        par(w) = 0.3 + 0.1*w - 0.02*w^2 + im*0.01*w^2
        for (k, w) in enumerate(f3)
            v3[1,1,k] = par(w)
        end
        p3 = ScatteringParameters((f3, v3))
        JosephsonCircuits.evaluatescattering!(d, p3, [3.1])
        @test d[1,1,1] ≈ par(3.1) atol = 1e-14
        p2 = ScatteringParameters(([1.0, 3.0],
            reshape(Complex{Float64}[0.2, 0.6], 1, 1, 2)))
        JosephsonCircuits.evaluatescattering!(d, p2, [2.0])
        @test d[1,1,1] ≈ 0.4 atol = 1e-14

        # callable provider requires nports
        f(w) = [0.0 exp(-im*w*1e-12); exp(-im*w*1e-12) 0.0]
        @test_throws ArgumentError ScatteringParameters(f)
        cblk = ScatteringParameters(f; nports = 2)
        d2 = zeros(Complex{Float64}, 2, 2, 1)
        JosephsonCircuits.evaluatescattering!(d2, cblk, [2*pi*5e9])
        @test d2[2,1,1] ≈ exp(-im*2*pi*5e9*1e-12)

        # passivity validation; an active block declares its noise, which
        # is held to the minimum the commutation relations require, here
        # |I - S S'| = 3 I
        @test_throws ArgumentError ScatteringParameters([0.0 2.0; 2.0 0.0])
        active = ScatteringParameters([0.0 2.0; 2.0 0.0];
            noise = NoiseCovariance([3.0 0.0; 0.0 3.0]))
        @test active.nports == 2
        @test active.noise.provider isa JosephsonCircuits.ConstantMatrixProvider
        @test_throws ArgumentError ScatteringParameters([0.0 2.0; 2.0 0.0];
            noise = NoiseCovariance([1.0 0.0; 0.0 1.0]))
        @test_throws ArgumentError NoiseCovariance([3.0 0.0; 0.0 3.0]; atol = -1.0)
        # noise covariance must be Hermitian
        @test_throws ArgumentError ScatteringParameters([0.0 2.0; 2.0 0.0];
            noise = NoiseCovariance([3.0 1.0; 0.0 3.0]))
        # thermal equilibrium noise model carries the temperature
        blkT = ScatteringParameters(S; noise = ThermalEquilibrium(20e-3))
        @test blkT.noise.temperature == 20e-3
        # a temperature is finite and nonnegative, and so is a table's
        # frequency, which a sorted check alone lets through
        @test_throws ArgumentError ThermalEquilibrium(-1.0)
        @test_throws ArgumentError Resistor(50.0; temperature = NaN)
        @test_throws ArgumentError Capacitor(1e-12; temperature = -2.0)
        @test Inductor(1e-9; temperature = 0).temperature === 0.0
        @test_throws ArgumentError ScatteringParameters(([1.0, NaN, 3.0], vals))
        @test_throws ArgumentError ScatteringParameters(([1.0, 2.0, Inf], vals))

        # reference impedance handling
        @test ScatteringParameters(S; zref = 30.0).zref == [30.0, 30.0]
        @test ScatteringParameters(S; zref = [50.0, 75.0]).zref == [50.0, 75.0]
        @test_throws DimensionMismatch ScatteringParameters(S; zref = [50.0])
        @test_throws ArgumentError ScatteringParameters(S; zref = -50.0)
    end

    @testset "transmission line" begin
        Z0, len = 50.0, 1e-3
        tl = TransmissionLine(Z0, len)
        @test tl isa ScatteringParameters
        @test tl.zref == [Z0, Z0]
        @test tl.negative_frequency isa Native
        w = 2*pi*5e9
        d = zeros(Complex{Float64}, 2, 2, 2)
        JosephsonCircuits.evaluatescattering!(d, tl, [w, -w])
        delay = len/JosephsonCircuits.speed_of_light
        @test d[2,1,1] ≈ exp(-im*w*delay)
        @test d[1,1,1] == 0
        # native evaluation satisfies the conjugation identity by construction
        @test d[:,:,2] ≈ conj.(d[:,:,1])
    end

    @testset "gaussian channels" begin
        # attenuator at zero temperature: exactly at the CP boundary
        η = 0.5
        ch = GaussianChannel(sqrt(η)*Matrix(1.0I, 2, 2),
            (1-η)/2*Matrix(1.0I, 2, 2); nmodes = 1)
        @test abs(ch.cp_margin) < 1e-10
        @test ch.nmodes == 1
        @test JosephsonCircuits.nterminals(ch) == 2
        # quantum limited phase insensitive amplifier
        G = 4.0
        amp = GaussianChannel(sqrt(G)*Matrix(1.0I, 2, 2),
            (G-1)/2*Matrix(1.0I, 2, 2); nmodes = 1)
        @test abs(amp.cp_margin) < 1e-10
        # a noiseless amplifier is not completely positive
        @test_throws ArgumentError GaussianChannel(
            sqrt(G)*Matrix(1.0I, 2, 2), zeros(2, 2); nmodes = 1)
        # Y must be symmetric
        @test_throws ArgumentError GaussianChannel(Matrix(1.0I, 2, 2),
            [0.5 0.1; -0.1 0.5]; nmodes = 1)
        # odd dimension is rejected
        @test_throws DimensionMismatch GaussianChannel(zeros(3,3),
            zeros(3,3))
        # ideal squeezer: symplectic X, Y = 0 is completely positive
        r = 0.5
        sq = GaussianChannel([exp(r) 0.0; 0.0 exp(-r)], zeros(2,2);
            nmodes = 1)
        @test abs(sq.cp_margin) < 1e-10
        # Bogoliubov conversion: a beamsplitter swap
        X = quadraturetransform([0 1; 1 0], zeros(2, 2))
        @test X == [0 1 0 0; 1 0 0 0; 0 0 0 1; 0 0 1 0]
        # two mode channel from the Bogoliubov form of a two mode squeezer
        A = cosh(r)*Matrix(1.0I, 2, 2)
        B = sinh(r)*[0.0 1.0; 1.0 0.0]
        Xtms = quadraturetransform(A, B)
        tms = GaussianChannel(Xtms, zeros(4, 4); nmodes = 2)
        @test abs(tms.cp_margin) < 1e-10
        # channels embed in circuits and are rejected by the solver bridge
        cg = Circuit([:ch => GaussianChannel(sqrt(η)*Matrix(1.0I, 2, 2),
                (1-η)/2*Matrix(1.0I, 2, 2); nmodes = 1, grounded = true)],
            [((:ch, 1), Ground)])
        @test_throws ComponentNotSupportedError compile(cg)
    end

    @testset "nonlinear inductors and current phase relations" begin
        jj = JosephsonJunction(100e-12)
        @test jj isa NonlinearInductor
        @test JosephsonCircuits.issinusoidal(jj)
        @test JosephsonJunction(Ic = JosephsonCircuits.phi0/100e-12).L0 ≈ 100e-12

        p = PolynomialCPR([1.0, 0.0, -1/6])
        @test p(0.1) ≈ 0.1 - 0.1^3/6
        # inductors written from equal expansions are equal, as a set of
        # them sees, and ones from different expansions are not
        twice = Set([NonlinearInductor(1e-9, p), NonlinearInductor(1e-9, PolynomialCPR([1.0, 0.0, -1/6]))])
        @test length(twice) == 1
        @test NonlinearInductor(1e-9, p) != NonlinearInductor(1e-9, PolynomialCPR([1.0, 0.0, -1/5]))
        dp = JosephsonCircuits.cprderivative(p)
        @test dp(0.0) ≈ 1.0
        @test dp(0.2) ≈ 1.0 - 0.2^2/2
        # the linear coefficient must be one
        @test_throws ArgumentError PolynomialCPR([2.0, 0.0])
        @test_throws ArgumentError PolynomialCPR(Float64[])
        # unknown callables require an explicit derivative
        mycpr(x) = x - x^3/6
        @test_throws ArgumentError NonlinearInductor(1e-9, mycpr)
        nl = NonlinearInductor(1e-9, mycpr, x -> 1 - x^2/2)
        @test !JosephsonCircuits.issinusoidal(nl)
        # a snail written as its expansion compiles as a junction, with
        # the relation kept beside the table; see test/nonlinearinductor.jl
        c = Circuit([:snail => NonlinearInductor(1e-9, p),
                     :r => Resistor(50.0), :p1 => Port(1; termination = nothing)],
            [((:p1, 1), (:r, 1), (:snail, 1)),
             ((:p1, 2), (:r, 2), (:snail, 2), Ground)])
        @test JosephsonCircuits.ninstances(elaborate(c)) == 3
        psc = compile(c)
        i = findfirst(==("snail"), psc.componentnames)
        @test psc.componenttypes[i] == :Lj
        @test psc.junctioncprs[i].a == [0.0, 1.0, 0.0, -1/6]
        # a relation which is neither the Josephson one nor a polynomial
        # is refused, since no solver can write it down
        cn = Circuit([:nl => nl, :r => Resistor(50.0),
                      :p1 => Port(1; termination = nothing)],
            [((:p1, 1), (:r, 1), (:nl, 1)),
             ((:p1, 2), (:r, 2), (:nl, 2), Ground)])
        @test_throws ComponentNotSupportedError compile(cn)
        # a sinusoidal junction lowers to the legacy Lj component
        c2 = Circuit([:jj => JosephsonJunction(100e-12), :p1 => Port(1; termination = nothing),
                      :r => Resistor(50.0)],
            [((:p1, 1), (:r, 1), (:jj, 1)),
             ((:p1, 2), (:r, 2), (:jj, 2), Ground)])
        psc = compile(c2)
        @test psc.componenttypes[findfirst(==("jj"),
            psc.componentnames)] == :Lj
    end

    @testset "touchstone loading" begin
        path = joinpath(mktempdir(), "attenuator.s2p")
        open(path, "w") do io
            write(io, "# GHz S MA R 50\n")
            write(io, "1.0 0.0 0.0 0.5 0.0 0.5 0.0 0.0 0.0\n")
            write(io, "2.0 0.0 0.0 0.5 0.0 0.5 0.0 0.0 0.0\n")
        end
        @test_throws ArgumentError ScatteringParameters(path; noise = Lossless())
        # a line is lossless at every frequency, so it may say so
        @test TransmissionLine(50.0, 1e-3; noise = Lossless()).noise isa Lossless
        blk = ScatteringParameters(path)
        @test blk.nports == 2
        @test blk.zref == [50.0, 50.0]
        d = zeros(Complex{Float64}, 2, 2, 1)
        JosephsonCircuits.evaluatescattering!(d, blk, [2*pi*1.5e9])
        @test d[2,1,1] ≈ 0.5
        # a conflicting explicit zref is an error
        @test_throws ArgumentError ScatteringParameters(path; zref = 30.0)
        # an explicit 50 Ohms is a statement like any other, not an
        # omission: against a 75 Ohm file it conflicts, and omitted the
        # file's value is taken
        path75 = joinpath(mktempdir(), "attenuator75.s2p")
        open(path75, "w") do io
            write(io, "# GHz S MA R 75\n")
            write(io, "1.0 0.0 0.0 0.5 0.0 0.5 0.0 0.0 0.0\n")
            write(io, "2.0 0.0 0.0 0.5 0.0 0.5 0.0 0.0 0.0\n")
        end
        @test ScatteringParameters(path75).zref == [75.0, 75.0]
        @test ScatteringParameters(path75; zref = 75.0).zref == [75.0, 75.0]
        @test_throws ArgumentError ScatteringParameters(path75; zref = 50.0)
        # Y and Z parameters are converted to S at the file's reference
        # impedances: a 100 ohm load, normalized to 50 ohms in a version 1
        # file and in ohms and siemens at 25 ohms in a version 2 one, and
        # a series 100 ohm resistor as the admittances of a two port, each
        # reproduce the resistor in a circuit
        dir = mktempdir()
        file(name, text) = (p = joinpath(dir, name); write(p, text); p)
        v2(param, ref, rows) = "[Version] 2.0\n# GHz $(param) RI R $(ref)\n" *
            "[Number of Ports] $(isqrt(length(split(rows[1])) ÷ 2))\n" *
            (length(split(rows[1])) == 8 ? "[Two-Port Data Order] 12_21\n" : "") *
            "[Number of Frequencies] 2\n[Network Data]\n1.0 $(rows[1])\n2.0 $(rows[1])\n[End]\n"
        wt = 2pi*1.5e9
        loaded(entry) = hblinsolve([wt], Circuit([(:p1, 1, 0, Port(1)), entry]); keyedarrays = false).S
        for path in (file("z1.s1p", "# GHz Z RI R 50\n1.0 2.0 0.0\n2.0 2.0 0.0\n"),
                file("y1.s1p", "# GHz Y RI R 50\n1.0 0.5 0.0\n2.0 0.5 0.0\n"),
                file("z2.s1p", v2("Z", 25, ["100.0 0.0"])),
                file("y2.s1p", v2("Y", 25, ["0.01 0.0"])))
            @test loaded((:x, 1, ScatteringParameters(path))) ≈ loaded((:x, 1, 0, Resistor(100.0))) atol = 1e-12
        end
        series(x) = hblinsolve([wt], Circuit([(:p1, 1, 0, Port(1)), (:x, 1, 2, x),
            (:p2, 2, 0, Port(2))]); keyedarrays = false).S
        yseries = file("y2.s2p", v2("Y", 50, ["0.01 0.0 -0.01 0.0 -0.01 0.0 0.01 0.0"]))
        @test series(ScatteringParameters(yseries)) ≈ series(Resistor(100.0)) atol = 1e-12
        # an open circuit written as an admittance reflects everything, and
        # hybrid parameters, which mix units, are refused
        dopen = zeros(ComplexF64, 1, 1, 1)
        JosephsonCircuits.evaluatescattering!(dopen, ScatteringParameters(
            file("open.s1p", "# GHz Y RI R 50\n1.0 0.0 0.0\n2.0 0.0 0.0\n")), [wt])
        @test dopen[1, 1, 1] ≈ 1
        @test_throws ArgumentError ScatteringParameters(file("h.s2p",
            "# GHz H RI R 50\n1.0 0 0 1 0 -1 0 0 0\n2.0 0 0 1 0 -1 0 0 0\n"))
    end

    @testset "one contract for a block's data" begin
        JC = JosephsonCircuits
        # the same data gets the same verdict stored, where the block is
        # checked when it is built, and as a callable, which the solve
        # checks at its frequency: a 40 dB amplifier's covariance which is
        # Hermitian to a part in 1e9 of its size, and one which is not
        w0 = 2pi*5e9
        amp = ComplexF64[0 0; 100 0]
        V = ComplexF64[2e5 1e5; 1e5+2e-4 2e5]
        two(b) = Circuit([(:p1, 1, 0, Port(1)), (:b, 1, 2, b), (:p2, 2, 0, Port(2))])
        verdict(f) = try
            f(); true
        catch e
            e isa ArgumentError || rethrow()
            false
        end
        for (Vk, accepted) in ((V, true), (V .+ ComplexF64[0 0; 1 0], false))
            stored = verdict(() -> ScatteringParameters(amp; noise = NoiseCovariance(Vk)))
            called = verdict(() -> hblinsolve([w0], two(ScatteringParameters(amp;
                noise = NoiseCovariance(w -> Vk))); keyedarrays = false, returnCnoise = true))
            @test stored == called == accepted
        end
        # a nearly lossless two port declared lossless is held to the
        # block's atol, as a constant and as a table
        g = sqrt(1 - 1e-7)
        nearly = ComplexF64[0 g; g 0]
        @test_throws ArgumentError ScatteringParameters(nearly; noise = Lossless())
        @test ScatteringParameters(nearly; noise = Lossless(), atol = 1e-6).atol == 1e-6
        @test ScatteringParameters(([1.0, 2.0], repeat(nearly, 1, 1, 2)); noise = Lossless(),
            extrapolation = :constant, atol = 1e-6) isa ScatteringParameters
        # and with the default Passive model the same tolerance says
        # whether it needs the channels of its loss
        @test JC.provablylossless(ScatteringParameters(nearly; atol = 1e-6))
        @test !JC.provablylossless(ScatteringParameters(nearly))
        # every kind of stored data is held to the contract: a piecewise
        # table as scattering data and as a covariance
        fs = collect(range(1.0, 5.0; length = 5))
        bands(A) = JC.piecewisetable(vcat(fs, fs .+ 100), cat(A, A; dims = 3); step = 1.0)
        @test_throws ArgumentError ScatteringParameters(bands(repeat(ComplexF64[0 2; 2 0], 1, 1, 5)))
        @test_throws ArgumentError ScatteringParameters(bands(repeat(ComplexF64[0 1; 1 0], 1, 1, 5));
            noise = Lossless())
        skew = repeat(ComplexF64[2 0.5; 0.5 2], 1, 1, 5)
        skew[1, 2, 3] += 1e-3
        @test_throws ArgumentError ScatteringParameters(ComplexF64[0 0.5; 0.5 0];
            noise = NoiseCovariance(bands(skew)))
        # A rotation of a constant covariance against a constant scattering
        # matrix meets the commutation relations at some frequencies and
        # not at others, so the construction leaves them to the solve,
        # which checks every frequency it evaluates the block at: this
        # one is the least covariance the block can emit at wq, where the
        # delay turns a quarter, and less than it must emit at zero and at
        # 2 wq
        S = ComplexF64[0.3 0.5; 0.5 0.3]
        K = I - S*S'
        taus = [0.4e-9, 0.0]
        wq = pi/2/taus[1]
        turn(w) = [cis(w*(taus[p] - taus[q])) for p in 1:2, q in 1:2]
        rotated = ScatteringParameters(S; noise = NoiseCovariance(
            JC.RotatedMatrixProvider(JC.ConstantMatrixProvider(K .* conj.(turn(wq))), taus)))
        @test JC.quantumnoisemargin(K .* conj.(turn(wq)), S) < -0.1
        @test abs(hblinsolve([wq], two(rotated); keyedarrays = false, returnCM = true).CM[1, 1]) ≈ 1 atol = 1e-8
        @test_throws ArgumentError hblinsolve([2wq], two(rotated); keyedarrays = false, returnCnoise = true)
    end

    @testset "a derivative is read as the block's data is" begin
        JC = JosephsonCircuits
        # tabulated, it is interpolated and extrapolated as the block's
        # data: beyond the table the block holds its end value, and so
        # does its derivative
        fs = collect(range(1.0, 5.0; length = 9))
        S = reshape(ComplexF64.(0.1 .+ 0.01 .* fs), 1, 1, :)
        dS = reshape(fill(0.01 + 0im, 9), 1, 1, :)
        blk = ScatteringParameters((fs, S); extrapolation = :constant, interpolation = :linear,
            derivatives = (x = (fs, dS),))
        at(p, w) = JC.evaluateprovider!(zeros(ComplexF64, 1, 1, 1), p, [w])[1]
        @test at(blk.derivatives.x, 6.0) == at(blk.derivatives.x, 5.0)
        @test blk.derivatives.x.interpolation == :linear
        # and a block called one entry at a time may state its derivative
        # as data, which is not called at all
        entry = ScatteringParameters((p, q, w) -> 0.1 + 0.01w; nports = 1, form = :entry,
            derivatives = (x = (fs, dS), y = (p, q, w) -> 0.01 + 0im))
        @test at(entry.derivatives.x, 2.0) ≈ 0.01
        @test at(entry.derivatives.y, 2.0) ≈ 0.01
    end

    @testset "a rotation moves the reference planes of a covariance" begin
        JC = JosephsonCircuits
        # its two sides carry opposite phases, as a correlation's do, so it
        # is refused as the scattering data of a block, ordinary or pumped,
        # and its phase turns and does not scale
        through = JC.TabulatedMatrixProvider([1.0, 2.0], repeat(ComplexF64[0 1; 1 0], 1, 1, 2))
        turned = JC.RotatedMatrixProvider(through, [1e-10, 0.0])
        @test_throws ArgumentError ScatteringParameters(turned)
        @test_throws ArgumentError LinearizedScattering([turned], 2pi*1e9; harmonics = [0], nports = 2)
        @test_throws ArgumentError JC.RotatedMatrixProvider(through, [0.0, 0.0]; phase = 2.0)
        @test JC.RotatedMatrixProvider(through, [0.0, 0.0]; phase = cis(0.7)).phase ≈ cis(0.7)
    end
end

@testset "shared signed provider evaluation and rational workspace" begin
    JC = JosephsonCircuits
    S = [0.1 0.2im; -0.2im 0.1]
    V = Matrix(1.0I,2,2)
    for rule in (ConjugateSymmetry(),Native())
        block = ScatteringParameters(S;noise=NoiseCovariance(V),negative_frequency=rule)
        frequencies = [-1.,0.,2.]
        scattering = zeros(ComplexF64,2,2,3)
        covariance = similar(scattering)
        buffer = Float64[]
        JC.evaluatescattering!(scattering,block,frequencies,buffer)
        JC.evaluatecovariance!(covariance,block,frequencies,buffer)
        for (i,w) in enumerate(frequencies)
            @test scattering[:,:,i] == (rule isa ConjugateSymmetry && w<0 ? conj.(S) : S)
            @test covariance[:,:,i] == V
        end
    end
    # a defective realization, which has no usable eigenvectors, goes
    # through the Schur factors; two workspaces on one set of factors are
    # independent
    A = [-2. 1.;0. -2.]; B = [1. 0.;0. 1.]; C = [0.2 0.1;0. 0.3]; D=zeros(2,2)
    rf = JC.resolventfactors(A,B)
    work1 = JC.ResolventWorkspace(rf); work2 = JC.ResolventWorkspace(rf)
    out1 = zeros(ComplexF64,2,2); out2=similar(out1)
    for w in (-3.,0.,1.,4.)
        JC.rationaltransfer!(out1,rf,C*rf.Z,D,w,work1)
        saved = copy(out1)
        JC.rationaltransfer!(out2,rf,C*rf.Z,D,w+1,work2)
        @test out1 == saved
        @test out1 ≈ D+C*((im*w*I-A)\B)
        @test out2 ≈ D+C*((im*(w+1)*I-A)\B)
    end
    provider = JC.RationalScatteringProvider(A,B,C,D)
    ws = collect(range(-4.,4.;length=21))
    out = zeros(ComplexF64,2,2,length(ws))
    JC.evaluateprovider!(out,provider,ws)
    @test all(out[:,:,k] ≈ D+C*((im*w*I-A)\B) for (k,w) in enumerate(ws))
    # the provider holds a copy of its realization and the factors of it
    # it takes when it is built, so it evaluates the realization it was
    # given whatever becomes of the caller's matrices; another
    # realization is another provider
    held = copy(out)
    A[1,1] = -3.
    JC.evaluateprovider!(out,provider,ws)
    @test out == held
    JC.evaluateprovider!(out,JC.RationalScatteringProvider(A,B,C,D),ws)
    @test all(out[:,:,k] ≈ D+C*((im*w*I-A)\B) for (k,w) in enumerate(ws))
    # complex pairs are 2 by 2 blocks of the real Schur form: a
    # realization of three damped resonances under a similarity which is
    # not orthogonal, against a dense solve at each frequency, at the
    # resonances and away from them, one frequency at a time and many
    blocks = zeros(6,6)
    for (k,(a,w0)) in enumerate(((0.1,1.0),(1e-4,2.5),(0.3,4.0)))
        blocks[2k-1:2k,2k-1:2k] = [-a w0; -w0 -a]
    end
    similarity = Matrix(1.0I,6,6) + 0.4*triu(ones(6,6),1)
    A2 = similarity\blocks*similarity
    B2 = [1. 0.; 0. 1.; 1. 1.; 0.5 -1.; 2. 0.; 0. 1.]
    C2 = [0.2 0.1 0. 0.3 -0.1 0.2; 0.1 0. 0.4 0.1 0.2 -0.3]
    D2 = [0.1 0.; 0. -0.2]
    p2 = JC.RationalScatteringProvider(A2,B2,C2,D2)
    wr = vcat(collect(range(-5.,5.;length=41)), [1.0,2.5,-2.5,4.0,2.5+1e-4])
    dense(w) = D2+C2*((im*w*I-A2)\B2)
    out2 = zeros(ComplexF64,2,2,length(wr))
    JC.evaluateprovider!(out2,p2,wr)
    @test all(isapprox(out2[:,:,k],dense(w);rtol=1e-10) for (k,w) in enumerate(wr))
    # the same transfer function at a physical frequency scale, the state
    # matrix and the frequencies scaled together
    sc = 2pi*1e10
    outsc = similar(out2)
    JC.evaluateprovider!(outsc,JC.RationalScatteringProvider(sc*A2,sc*B2,C2,D2),sc .* wr)
    @test all(isapprox(outsc[:,:,k],dense(w);rtol=1e-10) for (k,w) in enumerate(wr))
    one2 = zeros(ComplexF64,2,2,1)
    @test all(isapprox(JC.evaluateprovider!(one2,p2,[w])[:,:,1],dense(w);rtol=1e-10) for w in wr)
    @test sort(JC.schureigenvalues(p2.factors.T);by=x->(imag(x),real(x))) ≈
        sort(eigvals(A2);by=x->(imag(x),real(x))) rtol=1e-10
end
