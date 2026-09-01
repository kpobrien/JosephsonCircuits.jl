using Symbolics
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
        # strictly increasing frequencies required
        @test_throws ArgumentError ScatteringParameters(([2.0, 1.0],
            vals[:,:,1:2]))

        # callable provider requires nports
        f(w) = [0.0 exp(-im*w*1e-12); exp(-im*w*1e-12) 0.0]
        @test_throws ArgumentError ScatteringParameters(f)
        cblk = ScatteringParameters(f; nports = 2)
        d2 = zeros(Complex{Float64}, 2, 2, 1)
        JosephsonCircuits.evaluatescattering!(d2, cblk, [2*pi*5e9])
        @test d2[2,1,1] ≈ exp(-im*2*pi*5e9*1e-12)

        # passivity validation
        @test_throws ArgumentError ScatteringParameters([0.0 2.0; 2.0 0.0])
        active = ScatteringParameters([0.0 2.0; 2.0 0.0];
            noise = NoiseCovariance([1.0 0.0; 0.0 1.0]))
        @test active.nports == 2
        # noise covariance must be Hermitian
        @test_throws ArgumentError ScatteringParameters([0.0 2.0; 2.0 0.0];
            noise = NoiseCovariance([1.0 1.0; 0.0 1.0]))
        # thermal equilibrium noise model carries the temperature
        blkT = ScatteringParameters(S; noise = ThermalEquilibrium(20e-3))
        @test blkT.noise.temperature == 20e-3

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
        # snail-style direct specification elaborates but does not yet lower
        c = Circuit([:snail => NonlinearInductor(1e-9, p),
                     :r => Resistor(50.0), :p1 => Port(1; termination = nothing)],
            [((:p1, 1), (:r, 1), (:snail, 1)),
             ((:p1, 2), (:r, 2), (:snail, 2), Ground)])
        @test JosephsonCircuits.ninstances(elaborate(c)) == 3
        @test_throws ComponentNotSupportedError compile(c)
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
        blk = ScatteringParameters(path)
        @test blk.nports == 2
        @test blk.zref == [50.0, 50.0]
        d = zeros(Complex{Float64}, 2, 2, 1)
        JosephsonCircuits.evaluatescattering!(d, blk, [2*pi*1.5e9])
        @test d[2,1,1] ≈ 0.5
        # a conflicting explicit zref is an error
        @test_throws ArgumentError ScatteringParameters(path; zref = 30.0)
    end
end
