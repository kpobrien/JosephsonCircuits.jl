using JosephsonCircuits, Test, LinearAlgebra, Random

function testtransientiq(backend = JosephsonCircuits.CPU())
    device(x) = JosephsonCircuits.tobackend(backend, x)
    rng = Random.default_rng()
    times = collect(range(0.31e-9; step = 2e-12, length = 401))
    frequencies = [5e9, 8e9, 6.7e9]
    ports = [1, 2, 1]
    traces = randn(rng, 2, length(times))
    @testset "Finite-window I/Q: $backend" begin
        for window in (:hann, :rectangular)
            plan = transientiqplan(times, frequencies; duration = 0.2e-9,
                ports, window, stride = 7, phasereference = 0.17e-9, backend)
            measured = transientiq(plan, device(traces))
            expected = zeros(ComplexF64, size(measured))
            for c in eachindex(frequencies),
                (j, n) in enumerate(plan.ntaps:plan.stride:length(times))

                expected[c, j] = 2sum(plan.taps[k+1] * traces[ports[c], n-k] *
                                      exp(-2pi*im*frequencies[c]*(times[n-k]-plan.phasereference))
                for k in 0:(plan.ntaps-1))
            end
            @test Array(measured) ≈ expected rtol=3e-12 atol=1e-14
            @test plan.times == times[plan.ntaps:7:end]
            @test plan.centertimes ≈ plan.times .- 0.1e-9
            @test plan.groupdelay ≈ 0.1e-9
            @test plan.noisebandwidth ≈ sum(abs2, plan.taps)/(2plan.dt)
            at3db = sum(plan.taps[k+1]*exp(-2pi*im*plan.bandwidth3db*k*plan.dt)
            for k in 0:(plan.ntaps-1))
            @test abs2(at3db) ≈ 0.5 rtol=1e-12
            weights = randn(rng, ComplexF64, size(expected))
            gradient = device(zeros(size(traces)))
            transientiqvjp!(gradient, plan, device(weights))
            # exact as an identity, and checkable to the roundoff of the
            # two sums which meet in it, which cancel by an amount the
            # draw decides; the tolerance is measured from their terms
            cancellation = dot(abs.(Array(gradient)), abs.(traces)) +
                sum(abs, weights .* expected)
            @test dot(Array(gradient), traces) ≈ real(dot(weights, expected)) rtol=2e-12 atol=1e-12*cancellation
            direction = randn(rng, size(traces))
            epsilon = 1e-5
            plus = transientiq(plan, device(traces+epsilon*direction))
            minus = transientiq(plan, device(traces-epsilon*direction))
            @test dot(Array(gradient), direction) ≈ real(dot(weights, Array(plus-minus)))/(2epsilon) rtol=1e-8 atol=1e-9
            # Future samples cannot affect a measurement already available.
            changed = copy(traces)
            changed[:, (plan.ntaps+1):end] .+= 1e3
            @test Array(transientiq(plan, device(changed)))[:, 1] ≈ expected[:, 1] rtol=1e-9 atol=1e-12
            @test JosephsonCircuits.KernelAbstractions.get_backend(measured) == backend
        end

        # A rectangular window of 100 samples rejects the image exactly for
        # this tone (two full periods at twice the carrier frequency).
        t = collect((0:499) .* 2e-12)
        A, phi = 1.7, 0.37
        plan = transientiqplan(
            t, [5e9]; duration = 99*2e-12, window = :rectangular, backend)
        rf = reshape(A .* cospi.(2*5e9 .* t .+ phi/pi), 1, :)
        @test Array(transientiq(plan, device(rf))) ≈
              fill(A*exp(im*phi), 1, length(plan.times)) rtol=2e-13
        # Extra, unmeasured trace rows get an exactly zero pullback.
        grad = device(ones(3, length(t)))
        transientiqvjp!(grad, plan, device(ones(ComplexF64, 1, length(plan.times))))
        @test iszero(Array(grad)[2:3, :])
    end
end

# A CPU test set, as a function so that including this file in the GPU
# suite does not run it there.
function testiqcontract()
    @testset "I/Q input contract" begin
        t = collect((0:100) .* 1e-12)
        @test_throws ArgumentError transientiqplan(t, [5e9]; duration = 1e-9)
        @test_throws ArgumentError transientiqplan(t, [5e9]; duration = 20e-12, stride = 0)
        @test_throws ArgumentError transientiqplan(t, [5e9]; duration = 20e-12, window = :unknown)
        @test_throws ArgumentError transientiqplan(t, [0.0]; duration = 20e-12)
        @test_throws ArgumentError transientiqplan(t, [6e11]; duration = 20e-12)
        @test_throws ArgumentError transientiqplan(t .^ 2, [5e9]; duration = 20e-12)
        p = transientiqplan(t, [5e9]; duration = 20e-12)
        @test_throws DimensionMismatch transientiq(p, zeros(1, 100))
        @test_throws ArgumentError transientiq(p, zeros(ComplexF64, 1, 101))
        static = transientiqplan(t, [5e9]; duration = 20e-12, backend = JosephsonCircuits.CPU(; static = true))
        @test transientiq(static, ones(1, 101)) ≈ transientiq(p, ones(1, 101)) rtol=1e-13
    end
    return nothing
end

# the main suite runs the CPU backend by including this file; the GPU
# suite defines TRANSIENTBACKENDTESTS and calls the function with its device
if !@isdefined(TRANSIENTBACKENDTESTS)
    testtransientiq(JosephsonCircuits.CPU())
    testiqcontract()
end
