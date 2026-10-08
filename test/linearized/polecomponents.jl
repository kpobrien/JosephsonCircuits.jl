using Test, JosephsonCircuits, LinearAlgebra, SparseArrays

@testset "Pole component coverage" begin
    JC = JosephsonCircuits
    distance(a,b) = isempty(a) ? 0.0 : isempty(b) ? Inf : maximum(x -> minimum(abs.(x .- b)), a)
    same(a,b; tol = 1e-8) = length(a) == length(b) && max(distance(a,b),distance(b,a)) < tol

    @testset "constant hybrid blocks, shorts, opens and differential terminals" begin
        for S in (-0.5, 0.0, 0.5, 1.0)
            block = ScatteringParameters(fill(S,1,1); zref = 2.0)
            c = Circuit([(:s,1,block),(:c,1,0,Capacitor(0.5))])
            result = hbstability(c)
            @test result.converged
            @test result.poles ≈ [-(1-S)/(1+S)] atol = 1e-12
        end
        short = Circuit([(:s,1,ScatteringParameters(fill(-1.0,1,1); zref = 2.0)),
            (:c,1,0,Capacitor(0.5))])
        @test isempty(hbstability(short).poles)
        for reversed in (false,true)
            block = ScatteringParameters(zeros(1,1); zref = 2.0, grounded = false)
            a,b = reversed ? (2,1) : (1,2)
            c = Circuit([(:s,a,b,block),(:c1,1,0,Capacitor(1.0)),(:c2,2,0,Capacitor(1.0))])
            @test same(hbstability(c).poles, [0,-1])
        end
        # An ideal through between unequal reference impedances is a
        # transformer, not a wire: both right-side impedance factors matter.
        block = ScatteringParameters([0.0 1; 1 0]; zref = [1.0,4.0])
        c = Circuit([(:s,1,2,block),(:c,1,0,Capacitor(1.0)),(:r,2,0,Resistor(4.0))])
        @test hbstability(c).poles ≈ [-1.0]
    end

    @testset "rational states, including hidden modes" begin
        for w0 in (1.0, 2pi*8e9), scale in (0.2w0,1.0w0,10.0w0)
            block = RationalScattering(w0*[-2.0 0;0 -7], reshape([1.0,0],2,1),
                w0*[0.4 0], zeros(1,1); zref = 1.0)
            c = Circuit([(:s,1,block),(:r,1,0,Resistor(3.0))])
            result = hbstability(c; frequencyscale = scale)
            @test result.converged
            @test same(result.poles/w0, [-1.8,-7])
            @test count(result.internalonly) == 1
            @test result.internalonly[argmin(abs.(result.poles .+ 7w0))]
            q = hbstability(c; method = ShiftInvert(-1.7w0; nev = 2), frequencyscale = scale)
            @test same(q.poles/w0, result.poles/w0)
            @test q.converged
            # the hidden mode is internal whichever method finds it, its
            # node coordinates roundoff, and has no profile
            r = hbstability(c; method = ContourIntegral(-7w0, 1.0w0), frequencyscale = scale)
            for x in (q, r)
                k = argmin(abs.(x.poles .+ 7w0))
                @test count(x.internalonly) == 1 && x.internalonly[k] && iszero(x.nodevoltage[:, :, k])
            end
            @test same(hbstability(c).poles/w0, [-1.8,-7])
        end
    end

    @testset "bounded contour checks and rank saturation" begin
        c = Circuit([(:r,1,0,Resistor(2.0)),(:c,1,0,Capacitor(0.5)),(:l,1,0,Inductor(2.0))])
        expected = hbstability(c)
        result = hbstability(c; method = ContourIntegral(-0.5, 1.2; moments = 3))
        @test result.converged
        @test same(result.poles,expected.poles)
        @test result.method isa ContourIntegral
        @test !only(result.searches).saturated
        limited = hbstability(c; method = ContourIntegral(-0.5, 1.2; moments = 2))
        @test only(limited.searches).saturated
        @test !limited.converged
        empty = hbstability(c; method = ContourIntegral(3.0, 0.5))
        @test isempty(empty.poles)
        @test empty.converged
        for kw in ((quadrature = 2,), (moments = 0,), (probes = 0,), (refinements = 0,), (ranktol = 0,))
            @test_throws ArgumentError ContourIntegral(0, 2; kw...)
        end
    end

    @testset "exact delayed line and its infinite pole ladder" begin
        # rho_left = rho_right = 1/2, round-trip condition
        # (1/4)*exp(-2s*tau) = 1: s = -log(2)/tau + i*k*pi/tau.
        for tau in (1.0, 1e-10)
            c = Circuit([(:r1,1,0,Resistor(3.0)),(:r2,2,0,Resistor(3.0)),
                (:line,1,2,TransmissionLine(1.0,tau; vp = 1.0))])
            @test_throws ArgumentError hbstability(c)
            @test_throws ArgumentError hbstability(c; method = Monodromy())
            result = hbstability(c; method = ContourIntegral(-log(2)/tau, 4.0/tau), frequencyscale = 1/tau)
            @test result.converged
            @test same(result.poles*tau, [-log(2)+im*k*pi for k in -1:1])
            @test maximum(result.residuals) < 1e-9
            # Matched lines have no returning wave and no discrete modes.
            matched = Circuit([(:r1,1,0,Resistor(1.0)),(:r2,2,0,Resistor(1.0)),
                (:line,1,2,TransmissionLine(1.0,tau; vp = 1.0))])
            q = hbstability(matched; method = ContourIntegral(-log(2)/tau, 4.0/tau), frequencyscale = 1/tau)
            @test isempty(q.poles)
            @test q.converged
        end
    end

    @testset "explicit Laplace R, C, L and scattering models" begin
        # Two equivalent series R=2,L=1 branches, in parallel with a C.
        # s^2 + 2s + 1 = 0 is avoided by taking the fixed C to be 0.5.
        models = [Resistor(FrequencyDependent(LaplaceResponse(s -> 2+s))),
            Capacitor(FrequencyDependent(LaplaceResponse(s -> 1/(s*(2+s)))))]
        for model in models
            c = Circuit([(:x,1,0,model),(:c,1,0,Capacitor(0.5))])
            result = hbstability(c; method = ContourIntegral(-1.0, 1.5; moments = 3), frequencyscale = 1.0)
            @test result.converged
            @test same(result.poles, [-1+im,-1-im])
        end
        # L(s)=2/(1+s) has inverse inductance (1+s)/2: exactly a
        # constant L=2 and parallel R=2, with finite static stiffness.
        c = Circuit([(:l,1,0,Inductor(FrequencyDependent(LaplaceResponse(s -> 2/(1+s))))),
            (:c,1,0,Capacitor(0.5))])
        result = hbstability(c; method = ContourIntegral(-0.5, 1.2; moments = 3))
        @test result.converged
        @test same(result.poles, [-0.5+im*sqrt(3)/2,-0.5-im*sqrt(3)/2])
        # Laplace callbacks use s=i*w on the HB axis, with the existing
        # negative-frequency reflection for lumped element values.
        f = LaplaceResponse(s -> 2+s)
        @test f(3.0) == 2+3im
        @test f(-3.0) == 2-3im
        @test JC.laplacevalue(2*FrequencyDependent(f)+1, 3+im) == 11+2im
        # Analytic scattering callback and an explicit rational realization
        # have the same closed-loop pole (but only the latter retains states).
        block = ScatteringParameters(LaplaceResponse(s -> fill(0.4/(s+2),1,1)); nports = 1, zref = 1.0)
        c = Circuit([(:s,1,block),(:r,1,0,Resistor(3.0))])
        result = hbstability(c; method = ContourIntegral(-1.8, 0.5))
        @test result.converged
        @test result.poles ≈ [-1.8]
        for (form,f) in ((:entry, LaplaceResponse((p,q,s) -> 0.4/(s+2))),
                (:inplace, LaplaceResponse((M,s) -> fill!(M,0.4/(s+2)))))
            block = ScatteringParameters(f; nports = 1, zref = 1.0, form)
            c = Circuit([(:s,1,block),(:r,1,0,Resistor(3.0))])
            result = hbstability(c; method = ContourIntegral(-1.8, 0.5))
            @test result.converged
            @test result.poles ≈ [-1.8]
        end
        # Independent source perturbations vanish, even with a complex drive.
        driven = Circuit([(:r,1,0,Resistor(2.0)),(:c,1,0,Capacitor(0.5)),(:i,1,0,CurrentSource(2im))])
        @test hbstability(driven).poles ≈ [-1.0]
        for bad in (FrequencyDependent(w -> 2+im*w), real(FrequencyDependent(f)))
            c = Circuit([(:r,1,0,Resistor(bad)),(:c,1,0,Capacitor(0.5))])
            @test_throws ArgumentError hbstability(c; method = ContourIntegral(-1.0, 2.0))
        end
    end

    @testset "periodic scattering and harmonic profile matching" begin
        a, b = 0.2, 0.05+0.025im
        block = LinearizedScattering([fill(a,1,1),fill(b,1,1)],1.0;
            harmonics = [0,1], nports = 1, zref = 1.0, phase = 0.37,
            noise = NoiseCovariance([fill(2.0,1,1),zeros(1,1)]))
        c = Circuit([(:s,1,block),(:c,1,0,Capacitor(1.0))])
        @test_throws ArgumentError hbstability(c)
        expected = 1-2/sqrt((1+a)^2-4abs2(b))
        p = hbstability(c; pumpfrequency = 1.0, Nmodulationharmonics = (6,), method = DenseSpectrum())
        q = hbstability(c; pumpfrequency = 1.0, Nmodulationharmonics = (8,), method = DenseSpectrum())
        i = argmin(abs.(imag.(p.poles)))
        j = argmin(abs.(imag.(q.poles)))
        @test p.poles[i] ≈ expected atol = 1e-10
        @test q.poles[j] ≈ expected atol = 1e-10
        matches = JC.matchpoles(p,q; indices = [i], aliases = false)
        @test only(matches).index == j
        @test only(matches).overlap > 1-1e-10
        @test only(matches).distance < 1e-10
        restricted = JC.matchpoles(p,q; indices = [i], candidateindices = [j])
        @test only(restricted).index == j
        @test only(JC.matchpoles(p,q; indices = [i], candidateindices = Int[])).index == 0
        @test_throws ArgumentError JC.matchpoles(p,q; candidateindices = [j,j])
        # Alias alignment keeps physical sideband frequencies unchanged.
        shifted = findfirst(k -> abs(q.poles[k]-q.poles[j]-im) < 1e-8, eachindex(q.poles))
        @test !isnothing(shifted)
        m = JC.matchpoles(q,p; indices = [shifted])
        @test only(m).overlap > 0.99999
        @test only(m).distance < 1e-8
        @test length(unique(x.index for x in JC.matchpoles(p,q; indices = [i,i]))) == 2
        @test only(JC.matchpoles(p,q; indices = [i], maxdistance = 0.0)).index == 0
        @test JC.poleassignment([-0.9 -0.8; -0.85 -0.1]) == [2,1]
        @test_throws ArgumentError JC.matchpoles(p,q; minoverlap = 2)
        @test_throws ArgumentError hbstability(c; pumpfrequency = 0.7)
        # a block pumped within harmonic balance's tolerance of a harmonic
        # is commensurate, as harmonic balance takes it
        @test hbstability(c; pumpfrequency = 1 + 1e-8, Nmodulationharmonics = (2,), method = DenseSpectrum()).converged
    end

    @testset "rational converted filters versus analytic harmonic transfers" begin
        # Independently assembled internal-state and transfer-function
        # formulations must agree, including the sine sign, pump phase,
        # input sideband argument and a commensurate block pump.
        p0 = JC.RationalScatteringProvider(fill(-3.0,1,1),ones(1,1),fill(0.3,1,1),fill(0.1,1,1))
        pc = JC.RationalScatteringProvider(fill(-4.0,1,1),ones(1,1),fill(0.08,1,1),zeros(1,1))
        ps = JC.RationalScatteringProvider(fill(-5.0,1,1),ones(1,1),fill(0.06,1,1),zeros(1,1))
        Hstate = [p0, JC.ModulatedRationalProvider(pc,ps)]
        Hlaplace = [LaplaceResponse(s -> fill(0.1+0.3/(s+3),1,1)),
            LaplaceResponse(s -> fill(0.08/(s+4)+im*0.06/(s+5),1,1))]
        function circuit(H)
            block = LinearizedScattering(H,2.0; harmonics = [0,1], nports = 1,
                phase = 0.63, zref = 1.0, noise = NoiseCovariance([fill(2.0,1,1),zeros(1,1)]))
            Circuit([(:s,1,block),(:c,1,0,Capacitor(1.0)),(:r,1,0,Resistor(5.0))])
        end
        c, d = circuit(Hstate), circuit(Hlaplace)
        p = hbstability(c; pumpfrequency = 1.0, Nmodulationharmonics = (4,), method = DenseSpectrum())
        region = (center = -1.0, radius = 0.45)
        q = hbstability(d; method = ContourIntegral(region.center, region.radius), pumpfrequency = 1.0,
            Nmodulationharmonics = (4,))
        inside = findall(z -> abs(z-region.center) < region.radius, p.poles)
        @test p.converged && q.converged
        @test !isempty(inside)
        @test same(p.poles[inside],q.poles)
        matches = JC.matchpoles(p,q; indices = inside, aliases = false)
        @test all(x -> x.index > 0 && x.overlap > 1-1e-9, matches)
        @test length(p.stateedgeweights) == length(p.poles)
        @test all(x -> 0 <= x <= 1, p.stateedgeweights)
    end

    @testset "the period map of a pumped block, its envelope aside" begin
        # a pumped block's converting filters realized in time: its
        # envelope gates a transient record, not the periodic device, so
        # the map's poles are the same with and without one, and those of
        # the polynomial modulo the pump frequency
        # the states realized at a scale b, each input times b and output
        # over b, the same response
        provider(a, c, d, b) = JC.RationalScatteringProvider(fill(a,1,1),fill(b,1,1),fill(c/b,1,1),fill(d,1,1))
        circuit(envelope; b = 1.0) = Circuit([(:s,1,LinearizedScattering([provider(-3.0, 0.3, 0.1, b),
            JC.ModulatedRationalProvider(provider(-4.0, 0.08, 0.0, b), provider(-5.0, 0.06, 0.0, b))],2.0;
            harmonics = [0,1], nports = 1, phase = 0.63, zref = 1.0, envelope)),
            (:c,1,0,Capacitor(1.0)),(:r,1,0,Resistor(5.0))])
        on, ramped = (hbstability(circuit(envelope); pumpfrequency = 1.0, method = Monodromy(nev = :all, steps = 128))
            for envelope in (nothing, t -> min(1.0, t/2pi)))
        dense = hbstability(circuit(nothing); pumpfrequency = 1.0, Nmodulationharmonics = (8,), method = DenseSpectrum())
        strip(z) = complex(real(z), mod(imag(z)+0.5, 1.0)-0.5)
        @test ramped.poles == on.poles
        @test maximum(s -> minimum(abs.(strip.(dense.poles .- s))), on.poles) < 1e-4
        # at the scale 1e24 the same modes, none internal to the block, with
        # the same profiles and harmonics, each profile to within ten times
        # its multiplier's relative error bound: a tolerance taken from this
        # test, whose profiles differ by two to five times the bound, each
        # multiplier lying about its own size from the nearest; it bounds
        # no profile's error in general
        scaled = hbstability(circuit(nothing; b = 1e24); pumpfrequency = 1.0, method = Monodromy(nev = :all, steps = 128))
        @test scaled.poles ≈ on.poles rtol = 1e-6
        @test scaled.internalonly == on.internalonly && !any(on.internalonly)
        @test all(eachindex(on.poles)) do j
            isapprox(abs.(scaled.nodevoltage[:, :, j]), abs.(on.nodevoltage[:, :, j]);
                rtol = max(1e-6, 10*max(on.residuals[j], scaled.residuals[j])))
        end
    end

    @testset "native delay in a converting harmonic" begin
        tau, omega = 0.2, 1.7
        native = TransmissionLine(1.0,tau; vp=1.0).provider
        analytic = LaplaceResponse(s -> [0 exp(-tau*s); exp(-tau*s) 0])
        function circuit(provider, phase)
            block = LinearizedScattering([zeros(2,2),provider],2omega;
                harmonics=[0,1], nports=2, zref=[1.0,2.0], phase,
                noise=NoiseCovariance([10Matrix{Float64}(I,2,2),zeros(2,2)]))
            return Circuit([(:s,1,2,block),(:c1,1,0,Capacitor(1.0)),
                (:c2,2,0,Capacitor(0.7)),(:r1,1,0,Resistor(3.0)),(:r2,2,0,Resistor(4.0))])
        end
        for phase in (0.0,0.63)
            c,d = circuit(native,phase),circuit(analytic,phase)
            a,b = (JC.hbpolesystem(compile(x),Dict(),nothing,(2,),omega;pumpfrequency=omega) for x in (c,d))
            for z in (-0.4+0.7im,0.3-0.8im)
                @test JC.polematrix(a,z) ≈ JC.polematrix(b,z) atol=1e-14 rtol=1e-13
            end
        end
        region=(center=-0.8,radius=0.7)
        p,q = (hbstability(circuit(provider,0.63); method=ContourIntegral(region.center, region.radius), pumpfrequency=omega,
            Nmodulationharmonics=(2,), frequencyscale=omega) for provider in (native,analytic))
        @test p.converged && q.converged
        @test !isempty(p.poles)
        @test same(p.poles,q.poles)
        @test all(x -> x.index > 0 && x.overlap > 1-1e-9,
            JC.matchpoles(p,q; aliases=false))
    end

    @testset "a pumped block beside a pumped junction against harmonic balance" begin
        # The block's conversion and the junction's modulation act on the
        # same harmonics, so the block's pump phase moves the poles. Those
        # near the pump against the poles of an AAA rational fit
        # (Nakatsukasa, Sete and Trefethen, SIAM J. Sci. Comput. 40, 2018)
        # to the reflection hblinsolve computes with its own block stamps.
        function aaapoles(z, f; tol = 1e-12, mmax = 20)
            J = collect(eachindex(z))
            support, values, C = ComplexF64[], ComplexF64[], zeros(ComplexF64, length(z), 0)
            r, w = fill(sum(f)/length(f), length(z)), ComplexF64[]
            for _ in 1:mmax
                j = J[argmax(abs.(f[J] .- r[J]))]
                push!(support, z[j]); push!(values, f[j])
                deleteat!(J, findfirst(==(j), J))
                C = hcat(C, 1 ./ (z .- z[j]))
                w = svd((Diagonal(f)*C .- C*Diagonal(values))[J, :]).V[:, end]
                r = copy(f)
                r[J] = (C*(w .* values))[J] ./ (C*w)[J]
                norm(f .- r, Inf) <= tol*norm(f, Inf) && break
            end
            m = length(w)
            B = Matrix{ComplexF64}(I, m + 1, m + 1)
            B[1, 1] = 0
            return filter(isfinite, eigvals([0 transpose(w); ones(m) Diagonal(support)], B))
        end
        wp = 2pi*4.75001e9
        block = LinearizedScattering([fill(-0.6 + 0.0im, 1, 1), fill(0.3 + 0.1im, 1, 1)], wp;
            harmonics = [0, 2], nports = 1, zref = 50.0, atol = 10.0, phase = 0.63)
        c = Circuit([(:p, 1, 0, Port(1; Z0 = 50.0)), (:cc, 1, 2, Capacitor(100e-15)),
            (:jj, 2, 0, JosephsonJunction(1e-9)), (:cj, 2, 0, Capacitor(1e-12)),
            (:cb, 2, 3, Capacitor(100e-15)), (:b, 3, block)])
        pump = hbnlsolve((wp,), (16,), [(mode = (1,), port = 1, current = 0.00565e-6)], c)
        ws = collect(range(wp - 2pi*300e6, wp + 2pi*300e6; length = 401))
        lin = hblinsolve(ws, c; nonlinear = pump, Nmodulationharmonics = (6,), keyedarrays = false)
        k = lin.signalindex
        p = hbstability(c; nonlinear = pump, Nmodulationharmonics = (6,),
            method = ShiftInvert(1e6 + im*wp; nev = 6))
        near(poles) = filter(s -> abs(imag(s) - wp) < 2pi*300e6 && abs(real(s)) < 5e9, poles)
        @test p.converged
        @test same(near(p.poles)/wp, near(aaapoles(im .* ws, lin.S[k, k, :]))/wp; tol = 1e-7)
    end

    @testset "the period map with a rational block and floating nodes" begin
        # the JPA fed through a rational series inductor: the block's state,
        # its port currents and the floating fluxes of the nodes it joins,
        # against the dense spectrum modulo the pump frequency
        R0, L = 50.0, 2e-9
        a = 2R0/L
        u = [1.0, -1.0]
        series = RationalScattering(fill(-a, 1, 1), reshape(u, 1, 2), reshape(-a .* u, 2, 1), Matrix(1.0I, 2, 2); zref = R0)
        c = Circuit([(:p1, 1, 0, Port(1; Z0 = R0)), (:l, 1, 2, series), (:cc, 2, 3, Capacitor(100e-15)),
            (:jj, 3, 0, JosephsonJunction(1000e-12)), (:cj, 3, 0, Capacitor(1000e-15))])
        wp = 2pi*4.75001e9
        pump = hbnlsolve((wp,), (16,), [(mode = (1,), port = 1, current = 0.00565e-6)], c)
        m = hbstability(c; nonlinear = pump, method = Monodromy(nev = :all, steps = 128))
        d = hbstability(c; nonlinear = pump, method = DenseSpectrum())
        strip(z) = complex(real(z), mod(imag(z) + wp/2, wp) - wp/2)
        @test m.converged && length(m.poles) == 4 && only(m.searches).gauge == 2
        @test maximum(s -> minimum(abs.(strip.(d.poles .- s))), m.poles) < 1e-5*maximum(abs, m.poles)
        # a junction beside a block with a state its port does not see: that
        # state's mode, internal to the block, which the map places by its
        # states, at no frequency, and leaves without a profile
        hidden = RationalScattering([-0.5 0.0; 0.0 -0.07], reshape([1.0, 1.0], 2, 1), reshape([0.3, 0.0], 1, 2),
            fill(0.2, 1, 1); zref = 20.0)
        beside = Circuit([(:p, 1, 0, Port(1; Z0 = 20.0)), (:j, 1, 0, JosephsonJunction(1.0)), (:c, 1, 0, Capacitor(1.0)),
            (:b, 1, hidden)])
        besidepump = hbnlsolve((0.8,), (12,), [(mode = (1,), port = 1, current = 0.04*JC.phi0)], beside;
            method = Newton(), atol = 1e-12, keyedarrays = false)
        h = hbstability(beside; nonlinear = besidepump, Nmodulationharmonics = (4,), method = Monodromy(nev = :all))
        k = findfirst(h.internalonly)
        @test count(h.internalonly) == 1 && isapprox(h.poles[k], -0.07; atol = 1e-9) && iszero(h.nodevoltage[:, :, k])
        # a junction biased behind an ideal through, whose own equations
        # leave the direct current through it to the circuit: the current
        # the solution determines, and the map's poles against those of the
        # biased junction's resonator, -1/(2RC) +- im sqrt(cos(phi)/(Lj C) - 1/(2RC)^2)
        through = Circuit([(:p1, 1, 0, Port(1; Z0 = 20.0)), (:t, 1, 2, ScatteringParameters([0.0 1.0; 1.0 0.0]; zref = 1.0)),
            (:jj, 2, 0, JosephsonJunction(1.0)), (:cj, 2, 0, Capacitor(1.0))])
        bias = hbnlsolve((1.0,), (0,), [(mode = (0,), port = 1, current = 0.2*JC.phi0)], through; dc = true)
        @test only(bias.blockcurrents) ≈ [0.2*JC.phi0 -0.2*JC.phi0]
        b = hbstability(through; nonlinear = bias, method = Monodromy(nev = :all))
        s = -1/40 + im*sqrt(sqrt(0.96) - 1/1600)
        @test length(b.poles) == 2 && maximum(abs, sort(b.poles; by = imag) .- [conj(s), s]) < 1e-6
    end

    @testset "the period map through a pumped block" begin
        # A junction beside a block pumped at twice its pump, whose
        # conversion, H_1(s) = g (s + i w)/((s + a)(s + b)), takes the
        # pump's odd harmonics into each other and none into a conjugate's,
        # so that harmonic balance's orbit through the block is the block's
        # in time: the orbit gives the block's states, harmonic by harmonic,
        # and the map's modes, the junction's and the block's, are those of
        # the dense spectrum, which takes the block's conversions in the
        # polynomial, modulo the pump frequency.
        omega = 0.8
        p0 = JC.RationalScatteringProvider(fill(-3.0, 1, 1), ones(1, 1), fill(0.3, 1, 1), fill(0.1, 1, 1))
        a, b, g = 2.0, 4.0, 0.2
        A, B = [-a 0.0; 0.0 -b], ones(2, 1)
        pc = JC.RationalScatteringProvider(A, B, g .* [a/(a - b) b/(b - a)], zeros(1, 1))
        ps = JC.RationalScatteringProvider(A, B, (g*omega/(b - a)) .* [1.0 -1.0], zeros(1, 1))
        block = LinearizedScattering([p0, JC.ModulatedRationalProvider(pc, ps)], 2omega; harmonics = [0, 1], nports = 1,
            zref = 1.0)
        beside = Circuit([(:p, 1, 0, Port(1; Z0 = 20.0)), (:j, 1, 0, JosephsonJunction(1.0)), (:c, 1, 0, Capacitor(1.0)),
            (:b, 1, block)])
        pump = hbnlsolve((omega,), (12,), [(mode = (1,), port = 1, current = 0.04*JC.phi0)], beside;
            method = Newton(), atol = 1e-12, keyedarrays = false)
        m = hbstability(beside; nonlinear = pump, Nmodulationharmonics = (4,), method = Monodromy(nev = :all, steps = 128))
        dense = hbstability(beside; nonlinear = pump, Nmodulationharmonics = (8,), method = DenseSpectrum())
        strip(z) = complex(real(z), mod(imag(z) + omega/2, omega) - omega/2)
        apart(z) = minimum(abs.(strip.(dense.poles) .- strip(z)))
        @test m.converged && length(m.poles) == 7
        @test apart(m.poles[1]) < 1e-7 && maximum(apart, m.poles) < 1e-4
    end

    @testset "the period map through a line" begin
        # A junction behind a mismatched line from its port: the line's
        # history over its delay, sampled at the step, is part of the map's
        # state, and the map's least damped pole against the contour's,
        # which takes the line's exact delay in the polynomial of the
        # harmonics: the contour converged in the harmonics, the map in its
        # steps, its error falling at the rule's order and its rate's
        # estimated error the error itself.
        omega = 0.8
        behind = Circuit([(:p, 1, 0, Port(1; Z0 = 20.0)), (:line, 1, 2, TransmissionLine(10.0, 1.0; vp = 1.0)),
            (:j, 2, 0, JosephsonJunction(1.0)), (:c, 2, 0, Capacitor(1.0))])
        pump = hbnlsolve((omega,), (12,), [(mode = (1,), port = 1, current = 0.04*JC.phi0)], behind;
            method = Newton(), atol = 1e-12, keyedarrays = false)
        coarse, fine = (hbstability(behind; nonlinear = pump, Nmodulationharmonics = (4,), method = Monodromy(; nev = 2, steps))
            for steps in (64, 128))
        s = fine.poles[1]
        contours = [hbstability(behind; nonlinear = pump, Nmodulationharmonics = (H,), method = ContourIntegral(s, 0.2))
            for H in (6, 8)]
        @test all(r -> r.converged && length(r.poles) == 1, contours)
        exact = only(contours[2].poles)
        @test abs(exact - only(contours[1].poles)) < 1e-12*abs(s)
        @test abs(s - exact) < 1e-7*abs(s) && abs(s - exact) < abs(coarse.poles[1] - exact)/12
        @test coarse.rateerrors[1] ≈ abs(real(exact) - real(coarse.poles[1])) rtol = 1e-2
    end

    @testset "the period map places a line's own modes by its waves" begin
        # A lossless line of unit delay shorted at both ends, beside a
        # ported node it does not touch: its current and its standing
        # waves, at k pi, reach no node, and the map places each by the
        # waves leaving the line's ports, internal to the line, of a
        # profile of zero, its rate, zero, within its error.
        shorted = Circuit([(:p, 1, 0, Port(1; Z0 = 1.0)), (:c, 1, 0, Capacitor(1.0)),
            (:line, 0, 0, TransmissionLine(1.0, 1.0; vp = 1.0))])
        r = hbstability(shorted; pumpfrequency = 1.0, Nmodulationharmonics = (4,), method = Monodromy(steps = 128, nev = 8))
        @test all(r.internalonly) && all(iszero, r.nodevoltage)
        @test sort(abs.(imag.(r.poles))) ≈ pi .* [0, 1, 1, 2, 2, 3, 3, 4] atol = 1e-3
        @test all(abs.(real.(r.poles)) .<= r.rateerrors)
    end

    @testset "the rate error of modes a coarse map mixes" begin
        # A resonator whose stiffness the pump modulates strongly, the
        # block's states z1' = -d z1 + z2, z2' = -w^2 (1 + q cos 2t) z1 -
        # d z2 beside a resistor its output does not load: in 16 steps of
        # the period its two modes mix, and the least damped rate,
        # positive, is no more resolved than its error says; in 64 the
        # error is the distance to the dense spectrum's rate, stable.
        w, q, d = 3.6779661016949152, 1.0882352941176472, 0.16461152882205513
        none = JC.RationalScatteringProvider(zeros(0, 0), zeros(0, 1), zeros(1, 0), zeros(1, 1))
        cosine = JC.RationalScatteringProvider([-d 1.0; -w^2 -d], reshape([0.0, 1.0], 2, 1),
            reshape([-w^2*q, 0.0], 1, 2), zeros(1, 1))
        block = LinearizedScattering([none, JC.ModulatedRationalProvider(cosine, none)], 2.0; harmonics = [0, 1],
            nports = 1, zref = 1.0)
        resonator = Circuit([(:block, 1, block), (:load, 1, 0, Resistor(3.0))])
        coarse, resolved = (hbstability(resonator; pumpfrequency = 1.0, Nmodulationharmonics = (0,),
            method = Monodromy(; steps)) for steps in (16, 64))
        dense = hbstability(resonator; pumpfrequency = 1.0, Nmodulationharmonics = (16,), method = DenseSpectrum())
        k = argmin(abs.(dense.poles .- resolved.poles[1]))
        rate = real(dense.poles[k])
        @test dense.edgeweights[k] < 1e-6 && rate < 0 && real(coarse.poles[1]) > 0
        @test coarse.rateerrors[1] > real(coarse.poles[1])
        @test real(resolved.poles[1]) + resolved.rateerrors[1] < 0
        @test resolved.rateerrors[1] ≈ abs(rate - real(resolved.poles[1])) rtol = 1e-2
    end

    @testset "no invented temporal model" begin
        table = ScatteringParameters(([0.0,1.0],zeros(1,1,2)); zref = 1.0)
        bad = Circuit([(:s,1,table),(:r,1,0,Resistor(1.0))])
        @test_throws ArgumentError hbstability(bad; method = ContourIntegral(-1.0, 0.5))
        mismatch = ScatteringParameters(zeros(1,1); zref = 1.0, dcmodel = OpenDC())
        bad = Circuit([(:s,1,mismatch),(:r,1,0,Resistor(1.0))])
        @test_throws ArgumentError hbstability(bad)
        consistent = ScatteringParameters(ones(1,1); zref = 1.0, dcmodel = OpenDC())
        good = Circuit([(:s,1,consistent),(:c,1,0,Capacitor(1.0))])
        @test hbstability(good).poles ≈ [0.0]
    end

    @testset "a realization's scale under the period map" begin
        # one state realizing S(s) = 0.4/(s + 2) at any scale b of its
        # coordinate, behind a resistor of 3 which reflects 1/2: the closed
        # loop pole, 1 = 0.2/(s + 2) at s = -1.8, observed at the node with
        # the same profile whatever b
        mapped(b) = hbstability(Circuit([(:s, 1, RationalScattering(fill(-2.0, 1, 1), fill(b, 1, 1),
            fill(0.4/b, 1, 1), zeros(1, 1); zref = 1.0)), (:r, 1, 0, Resistor(3.0))]);
            pumpfrequency = 1.0, Nmodulationharmonics = (0,), method = Monodromy(nev = :all))
        unit, scaled = mapped(1.0), mapped(1e24)
        @test only(unit.poles) ≈ -1.8 rtol = 1e-5
        @test scaled.poles ≈ unit.poles
        @test !any(unit.internalonly) && !any(scaled.internalonly)
        @test abs.(scaled.nodevoltage) ≈ abs.(unit.nodevoltage)
        @test norm(unit.nodevoltage) > 0
    end
end
