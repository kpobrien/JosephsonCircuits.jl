using JosephsonCircuits
using LinearAlgebra
using Test

# The quantum noise of a pumped transient against harmonic balance: a
# Josephson amplifier's gain and quantum efficiency on both rules and on
# a batch, a traveling wave amplifier's, and the convergence of a pulsed
# lossy line where no stationary reference exists. The passive and the
# warm circuits are in transientnoise.jl; this file is separate so that
# the two run on different workers.
@testset "the quantum noise of a pumped transient" begin
    JC = JosephsonCircuits

    @testset "a pumped amplifier against harmonic balance" begin
        circuit = [("P1", "1", "0", 1.0), ("R1", "1", "0", 50.0), ("C1", "1", "2", 100e-15),
            ("Lj1", "2", "0", 1e-9), ("C2", "2", "0", 1e-12)]
        fp, fs, ip = 4.75e9, 4.74e9, 0.00565e-6
        ramp(t) = t <= 0 ? 0.0 : t >= 2e-9 ? 1.0 : (1 - cospi(t/2e-9))/2
        # a physical cosine of peak 2ip has the positive frequency coefficient ip
        pump(t) = 2ip*ramp(t)*cospi(2fp*t)
        prob = transientproblem(circuit; sources = [TransientSource(1, pump)])
        # the trapezoidal rule warps the drive frequencies by (2 pi f dt)^2/12,
        # which the resonator's detuning and the amplifier's bifurcation
        # magnify: the gain converges to harmonic balance at second order,
        # 12% high at 5 ps and 2% at 0.3 ps, so the step is 2000 samples of
        # a pump period; the pumped resonator rings for the whole settle
        settle, record, dt = 200e-9, 100e-9, 0.15625e-12
        sol = transientsolve(prob, (0.0, settle + record - dt); dt, record = :phases)
        # the measurement on the settled window of the record; the baths are
        # driven from the equilibrium start through the settling
        first = round(Int, settle/dt) + 1
        times = sol.times[first:end]
        @test length(times) == round(Int, record/dt)
        measurement = transientquantumplan(times, [fs])
        # the stationary Floquet frequencies of the pumped response
        frequencies = sort!(abs.([fs + 2k*fp for k in -2:2]))
        noise = transientnoise(sol, measurement; frequencies, weights = fill(1/record, 5),
            inputs = measurement, commutationrtol = 3e-3)
        hb = hbsolve([2pi*fs], (2pi*fp,), [(mode = (1,), port = 1, current = ip)], (8,), (10,),
            circuit, Dict(); ftol = 5e-17)
        s = hb.linearized.S((0,), 1, (0,), 1, 1)
        qe = hb.linearized.QE((0,), 1, (0,), 1, 1)
        metrics = transientquantumefficiency(noise.gain, noise.covariance; rtol = 3e-3)
        @test abs2(s) > 1.1
        @test noise.diagnostics.passed
        @test metrics.gain ≈ abs2(s) rtol=1e-2
        @test metrics.QE ≈ qe rtol=3e-3
        @test metrics.normalizedQE ≈ 1 rtol=3e-3
        # the Gauss-Legendre rule at 84 samples per period, a sixteenth of
        # the trapezoidal step, on the same window: the gain to a part in
        # ten thousand and the forward method to roundoff of the adjoint
        dt = 2.5e-12
        gsol = transientsolve(prob, (0.0, settle + record - dt); dt, record = :phases, method = GaussLegendre())
        first = round(Int, settle/dt) + 1
        gmeasurement = transientquantumplan(gsol.times[first:end], [fs])
        gnoise = transientnoise(gsol, gmeasurement; frequencies, weights = fill(1/record, 5),
            inputs = gmeasurement, commutationrtol = 3e-3)
        gmetrics = transientquantumefficiency(gnoise.gain, gnoise.covariance; rtol = 3e-3)
        @test gnoise.diagnostics.passed
        @test gmetrics.gain ≈ abs2(s) rtol=1e-4
        @test gmetrics.QE ≈ qe rtol=1e-4
        gforward = transientnoise(gsol, gmeasurement; frequencies, weights = fill(1/record, 5),
            inputs = gmeasurement, commutationrtol = 3e-3, method = :forward)
        @test gforward.gain ≈ gnoise.gain rtol=1e-9
        @test gforward.covariance ≈ gnoise.covariance rtol=1e-9
        # the noise of three pump amplitudes on one pass, equal to each
        # member's, with the gain rising with the pump. The checks from
        # here are of equivalences, the batch to its members, the tiled
        # contraction to the whole, the forward method to the adjoint and
        # the checkpointed record to the kept one, so a short window with
        # the amplifier still ringing serves them: 20 ns settled and one
        # 20 ns beat recorded, the signal 50 MHz below the pump so that
        # the beat is on a bin, and the commutation tolerance loosened for
        # the ringing
        make(a) = transientproblem(prob; sources = [TransientSource(1, let a = a; t -> 2a*ramp(t)*cospi(2fp*t); end)])
        settle, record, fs = 20e-9, 20e-9, 4.70e9
        frequencies = sort!(abs.([fs + 2k*fp for k in -2:2]))
        batch = transientsolve(make.([0.004e-6, 0.005e-6, ip]), (0.0, settle + record - dt); dt, record = :phases)
        first = round(Int, settle/dt) + 1
        gmeasurement = transientquantumplan(batch.times[first:end], [fs])
        bnoise = transientnoise(batch, gmeasurement; frequencies, weights = fill(1/record, 5),
            inputs = gmeasurement, commutationrtol = 5e-2)
        @test size(bnoise.covariance) == (2, 2, 3) && size(bnoise.gain) == (2, 2, 3) && length(bnoise.diagnostics) == 3
        gains = Float64[]
        for j in 1:3
            single = transientnoise(batch[j], gmeasurement; frequencies, weights = fill(1/record, 5),
                inputs = gmeasurement, commutationrtol = 5e-2)
            @test bnoise.covariance[:, :, j] ≈ single.covariance rtol=1e-9
            @test bnoise.gain[:, :, j] ≈ single.gain rtol=1e-9
            @test bnoise.diagnostics[j].passed
            push!(gains, sum(abs2, bnoise.gain[:, :, j])/2)
        end
        @test issorted(gains)
        # the contraction tiled over the conditions and then over the
        # frequencies, within a budget too small for more than one of
        # either, is the same contraction in more passes
        JC.noisememorybudget[] = 1
        tnoise = transientnoise(batch, gmeasurement; frequencies, weights = fill(1/record, 5),
            inputs = gmeasurement, commutationrtol = 5e-2)
        JC.noisememorybudget[] = 0
        @test tnoise.covariance ≈ bnoise.covariance rtol=1e-9
        @test tnoise.commutator ≈ bnoise.commutator rtol=1e-9
        @test tnoise.gain ≈ bnoise.gain rtol=1e-9
        # a range of conditions is a batch of views
        pair = batch[2:3]
        @test pair isa TransientBatchSolution && length(pair) == 2
        @test transientnoise(pair, gmeasurement; frequencies, weights = fill(1/record, 5),
            inputs = gmeasurement, commutationrtol = 5e-2).covariance ≈ bnoise.covariance[:, :, 2:3] rtol=1e-9
        # the forward method on the batch is one tangent over every
        # condition, from each condition's stationary response
        fnoise = transientnoise(batch, gmeasurement; frequencies, weights = fill(1/record, 5),
            inputs = gmeasurement, commutationrtol = 5e-2, method = :forward)
        @test fnoise.covariance ≈ bnoise.covariance rtol=1e-9
        @test fnoise.gain ≈ bnoise.gain rtol=1e-9
        # the noise of a checkpointed record equals the phase record's
        ckpt = transientsolve(make.([0.004e-6, 0.005e-6, ip]), (0.0, settle + record - dt); dt, record = :checkpoints)
        cnoise = transientnoise(ckpt, gmeasurement; frequencies, weights = fill(1/record, 5),
            inputs = gmeasurement, commutationrtol = 5e-2)
        @test cnoise.covariance ≈ bnoise.covariance rtol=1e-6
        @test cnoise.gain ≈ bnoise.gain rtol=1e-6
        # the pulsed gain of the batch equals the members', and is the
        # periodic gain up to the amplifier's ring down inside the window
        bgain = transientgain(batch, gmeasurement, gmeasurement)
        @test size(bgain) == (2, 2, 3)
        for j in 1:3
            @test bgain[:, :, j] ≈ transientgain(batch[j], gmeasurement, gmeasurement) rtol=1e-9
        end
        @test bgain[:, :, 3] ≈ bnoise.gain[:, :, 3] rtol=0.2
    end

    @testset "a traveling wave amplifier against harmonic balance" begin
        # two hundred series junctions pumped to 2.6 radians of node flux,
        # a line rather than a resonator: the accumulated propagation phase
        # and the phase matching are the constraint, and the periodic gain
        # and the quantum efficiency agree with hbsolve at fourth order
        cells = 200
        circuit = Any[("P1", "1", "0", 1.0), ("R1", "1", "0", 50.0)]
        for k in 1:cells
            push!(circuit, ("Lj$k", "$k", "$(k+1)", 100e-12))
            push!(circuit, ("Cj$k", "$k", "$(k+1)", 300e-15))
            push!(circuit, ("C$k", "$(k+1)", "0", 50e-15))
        end
        push!(circuit, ("P2", "$(cells+1)", "0", 2.0))
        push!(circuit, ("R2", "$(cells+1)", "0", 50.0))
        fp, fs, ip = 7e9, 7.3e9, 1.5e-6
        hb = hbsolve([2pi*fs], (2pi*fp,), [(mode = (1,), port = 1, current = ip)], (8,), (10,), circuit, Dict(); ftol = 1e-14)
        s21 = hb.linearized.S((0,), 2, (0,), 1, 1)
        qe = hb.linearized.QE((0,), 2, (0,), 1, 1)
        @test maximum(abs.(hb.nonlinear.nodeflux)) > 2
        ramp(t) = t <= 0 ? 0.0 : t >= 2e-9 ? 1.0 : (1 - cospi(t/2e-9))/2
        prob = transientproblem(circuit; sources = [TransientSource(1, t -> 2ip*ramp(t)*cospi(2fp*t))])
        settle, record, dt = 20e-9, 20e-9, 2.5e-12
        sol = transientsolve(prob, (0.0, settle + record - dt); dt, method = GaussLegendre(), record = :checkpoints)
        first = round(Int, settle/dt) + 1
        measurement = transientquantumplan(sol.times[first:end], [fs]; ports = [2])
        inputs = transientquantumplan(sol.times[first:end], [fs]; ports = [1])
        frequencies = sort!(abs.([fs + 2k*fp for k in -3:3]))
        noise = transientnoise(sol, measurement; frequencies, weights = fill(1/record, 7), inputs, commutationrtol = 1e-2)
        @test noise.diagnostics.passed
        metrics = transientquantumefficiency(noise.gain, noise.covariance; rtol = 1e-2)
        @test metrics.gain ≈ abs2(s21) rtol=1e-4
        @test metrics.QE ≈ qe rtol=1e-4
        # the signed quadrature block of the transmission, phase included
        @test noise.gain ≈ [real(s21) imag(s21); -imag(s21) real(s21)] rtol=1e-3
    end

    @testset "a pulsed lossy line converges in the step, the bath spacing and the cutoff" begin
        # a hundred series junctions with a loss resistor in every cell,
        # pumped by a pulse that rises over 2 ns, holds 4 ns and falls
        # over 2 ns, measured in Hann windows of 2 ns on the rise, the
        # plateau and the fall, the output window delayed by the line's
        # transit; the bath is the whole band. There is no stationary
        # reference for this: the evidence is convergence, of the step at
        # fourth order, of the bath spacing, and of the cutoff, which must
        # pass the junction plasma frequency near 29 GHz where the line
        # responds most, beyond which nothing changes
        cells = 100
        circuit = Any[("P1", "1", "0", 1.0), ("R1", "1", "0", 50.0)]
        for k in 1:cells
            push!(circuit, ("Lj$k", "$k", "$(k+1)", 100e-12))
            push!(circuit, ("Cj$k", "$k", "$(k+1)", 300e-15))
            push!(circuit, ("C$k", "$(k+1)", "0", 50e-15))
            k < cells && push!(circuit, ("R$(k+2)", "$(k+1)", "0", 10e3))
        end
        push!(circuit, ("P2", "$(cells+1)", "0", 2.0))
        push!(circuit, ("R2", "$(cells+1)", "0", 50.0))
        fp, fs, ip = 7e9, 8e9, 1.5e-6
        pulse(t) = t <= 0 ? 0.0 : t < 2e-9 ? (1 - cospi(t/2e-9))/2 : t <= 6e-9 ? 1.0 : t < 8e-9 ? (1 + cospi((t - 6e-9)/2e-9))/2 : 0.0
        prob = transientproblem(circuit; sources = [TransientSource(1, t -> 2ip*pulse(t)*cospi(2fp*t))])
        T, delay = 2e-9, 0.25e-9
        solutions = Dict{Float64,Any}()
        function study(dt, spacing, fmax)
            sol = get!(() -> transientsolve(prob, (0.0, 10e-9 - dt); dt, method = GaussLegendre(), record = :checkpoints), solutions, dt)
            map((0.5e-9, 3.5e-9, 6e-9)) do start
                first = round(Int, start/dt) + 1
                nm, shift = round(Int, T/dt), round(Int, delay/dt)
                tin, tout = sol.times[first:first + nm - 1], sol.times[first + shift:first + shift + nm - 1]
                env = reshape(sinpi.((tin .- tin[1]) ./ T) .^ 2, :, 1)
                measurement = transientquantumplan(tout, [fs]; ports = [2], envelopes = env)
                inputs = transientquantumplan(tin, [fs]; ports = [1], envelopes = env)
                freqs = collect(spacing:spacing:fmax)
                noise = transientnoise(sol, measurement; frequencies = freqs, weights = fill(spacing, length(freqs)), commutationrtol = 5e-2)
                pulsed = transientgain(sol, measurement, inputs)
                @test noise.diagnostics.passed
                (gain = sum(abs2, pulsed)/2, covariance = noise.covariance, pulsed)
            end
        end
        coarse, fine, finest = study(5e-12, 0.5e9, 30e9), study(2.5e-12, 0.5e9, 30e9), study(1.25e-12, 0.5e9, 30e9)
        for w in 1:3
            # the step: fourth order in the gain and the covariance
            ratio = norm(coarse[w].covariance - fine[w].covariance)/norm(fine[w].covariance - finest[w].covariance)
            @test 8 < ratio < 32
            @test 8 < abs(coarse[w].gain - fine[w].gain)/abs(fine[w].gain - finest[w].gain) < 32
            @test norm(fine[w].covariance - finest[w].covariance)/norm(finest[w].covariance) < 1e-5
        end
        # the plateau is stationary enough for a scalar quantum efficiency
        plateau = transientquantumefficiency(finest[2].pulsed, finest[2].covariance; rtol = 5e-2)
        @test 0.5 < plateau.QE < 0.7 && 0.6 < plateau.gain < 0.7
        # the spacing: the change halves or better as the spacing halves
        half, quarter = study(5e-12, 0.25e9, 30e9), study(5e-12, 0.125e9, 30e9)
        for w in 1:3
            @test norm(half[w].covariance - quarter[w].covariance) < norm(coarse[w].covariance - half[w].covariance)
            @test norm(half[w].covariance - quarter[w].covariance)/norm(quarter[w].covariance) < 1e-4
        end
        # the cutoff: the band below the plasma frequency misses a part in
        # a thousand, and past it the covariance no longer moves
        low, wide, wider = study(5e-12, 0.5e9, 20e9), study(5e-12, 0.5e9, 40e9), study(5e-12, 0.5e9, 60e9)
        for w in 1:3
            @test norm(low[w].covariance - coarse[w].covariance)/norm(coarse[w].covariance) > 1e-4
            @test norm(wide[w].covariance - wider[w].covariance)/norm(wider[w].covariance) < 1e-7
        end
    end
end
