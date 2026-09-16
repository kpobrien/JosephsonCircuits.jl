using JosephsonCircuits, Test, LinearAlgebra, Random

# two ports declared out of numerical order, so that a port number is not
# the row of its trace: port 2 is row 1 and port 1 is row 2
quantumtestproblem() = transientproblem(Circuit([(:p2, 1, 0, Port(2)), (:c2, 1, 0, Capacitor(1e-12)),
    (:p1, 2, 0, Port(1)), (:c1, 2, 0, Capacitor(1e-12))]))

function testtransientquantum(backend = JosephsonCircuits.CPU())
    JC = JosephsonCircuits
    device(x) = JC.tobackend(backend, x)
    prob = quantumtestproblem()
    n, dt = 64, 1e-11
    times = collect(0:(n-1)) .* dt .+ 0.17e-9
    period = n*dt
    f = 3/period
    c = zeros(ComplexF64, fld(n-1, 2), 3)
    c[3, :] = [1, im, 1]
    plan = transientquantumplan(prob, times, c; ports = [1, 1, 2], backend)
    rng = Random.default_rng()
    @testset "Quantum temporal modes: $backend" begin
        @test plan.rows == [2, 2, 1]
        X, P = 1.7, -0.4
        rf = sqrt(JC.planck_constant*f/period) .* (X .* cospi.(2f .* (times .- times[1])) .+
              P .* sinpi.(2f .* (times .- times[1])))
        # port 1 carries the wave and port 2 twice it, on their rows
        traces = device(permutedims(hcat(2rf, rf)))
        @test Array(transientquantum(plan, traces)) ≈ [X, P, P, -X, 2X, 2P] rtol=1e-13
        @test plan.gram ≈ [1 im 0; -im 1 0; 0 0 1]
        # Independent canonical Fourier quadratures reconstruct both the
        # covariance and commutators, including complex overlapping modes.
        response = zeros(6, 4length(plan.frequencies))
        for port in 1:2, (k, fk) in enumerate(plan.frequencies), q in 1:2
            wave = zeros(2, n)
            phase = 2fk .* (times .- times[1])
            wave[port, :] .= sqrt(JC.planck_constant*fk/period) .*
                             (q==1 ? cospi.(phase) : sinpi.(phase))
            col = 2length(plan.frequencies)*(port-1)+2k-2+q
            response[:, col] .= Array(transientquantum(plan, device(wave)))
        end
        v, k = device(zeros(6, 6)), device(zeros(6, 6))
        JC.transientnoiseaccumulate!(v, k, device(response), device(fill(0.5, size(response, 2))))
        @test Array(v) ≈ plan.vacuum atol=5e-14
        @test Array(k) ≈ plan.commutator atol=5e-14
        @test transientquantumdiagnostics(v, k, plan.commutator).passed
        weights = randn(rng, 6)
        gradient = device(zeros(3, n))
        transientquantumvjp!(gradient, plan, device(weights))
        direction = randn(rng, 3, n) .* 1e-6
        @test dot(Array(gradient), direction) ≈
              dot(weights, Array(transientquantum(plan, device(direction)))) rtol=1e-13
        @test iszero(Array(gradient)[3, :])
        envelope = reshape(sinpi.((0:(n-1)) ./ n) .^ 2, :, 1)
        windowed = transientquantumplan(prob, times, [3.4/period]; envelopes = envelope, backend)
        @test norm(windowed.coefficients) ≈ 1
        @test windowed.vacuum ≈ Matrix(0.5I, 2, 2)
        @test JosephsonCircuits.KernelAbstractions.get_backend(transientquantum(plan, traces)) == backend

        # Analytic Gaussian channels fix the amplifier and attenuator factors.
        single = transientquantumplan(prob, times, [f]; backend)
        J = [0.0 1; -1 0]
        for g in (1.0, 2.0, 100.0)
            h = hcat(sqrt(g)*Matrix(1.0I, 2, 2), sqrt(g-1)*[1.0 0; 0 -1])
            v, k = device(zeros(2, 2)), device(zeros(2, 2))
            JC.transientnoiseaccumulate!(v, k, device(h), device(fill(0.5, 4)))
            @test transientquantumdiagnostics(v, k, J).passed
            metrics = transientquantumefficiency(h[:, 1:2], v)
            @test metrics.addednoise ≈ (1-1/g)/2 atol=1e-14
            @test metrics.QE ≈ g/(2g-1)
            @test metrics.normalizedQE ≈ 1
        end
        @test !transientquantumdiagnostics(2Matrix(1.0I, 2, 2), 4J, J).passed
        @test !transientquantumdiagnostics(0.4Matrix(1.0I, 2, 2), J, J).passed
        eta = 0.3
        @test transientquantumefficiency(sqrt(eta)*Matrix(1.0I, 2, 2), single.vacuum).QE ≈
              eta
        @test_throws ArgumentError transientquantumefficiency([2.0 0; 0 0.5], single.vacuum)
    end
end

# A CPU test set, as a function so that including this file in the GPU
# suite does not run it there.
function testquantumcontracts()
    @testset "Quantum measurement and bath contracts" begin
        JC = JosephsonCircuits
        ts = collect(0:63) .* 1e-11
        two = quantumtestproblem()
        @test_throws ArgumentError transientquantumplan(two, ts, [0.0])
        @test_throws ArgumentError transientquantumplan(two, ts, [0.5e11])
        @test_throws ArgumentError transientquantumplan(two, ts, [3.4/(64e-11)])
        @test_throws ArgumentError transientquantumplan(two, ts .^ 2, [1e9])
        @test_throws ArgumentError transientquantumplan(two, ts, zeros(31, 1))
        # a mode names a port by its number, which must exist
        @test_throws ArgumentError transientquantumplan(two, ts, [1/(64e-11)]; ports = [3])
        @test transientquantumplan(two, ts, [1/(64e-11)]; ports = [2]).rows == [1]
        circuit = Circuit(
            [:p => Port(1), :Rloss => Resistor(100.0; temperature = 0.1),
                :Ropen => Resistor(Inf), :c => Capacitor(1e-12)],
            [((:p, 1), (:Rloss, 1), (:Ropen, 1), (:c, 1)),
                ((:p, 2), (:Rloss, 2), (:Ropen, 2), (:c, 2), Ground)])
        prob = transientproblem(circuit)
        # the port takes the default temperature and the resistor its own
        baths = transientnoisebaths(prob; temperature = 0.05)
        @test length(baths.channels) == 2
        @test [b.port for b in baths.channels] == [1, 0]
        @test [b.temperature for b in baths.channels] == [0.05, 0.1]
        @test baths.channels[2].name == "Rloss"
        @test_throws ArgumentError transientnoisebaths(prob; temperature = -1)
        unmatched = Circuit(
            [:p=>Port(1; termination = nothing), :r=>Resistor(50.0),
                :c=>Capacitor(1e-12)],
            [((:p, 1), (:r, 1), (:c, 1)), ((:p, 2), (:r, 2), (:c, 2), Ground)])
        @test_throws ArgumentError transientnoisebaths(transientproblem(unmatched))
        for T in (0.0, 300.0)
            bath = JC.TransientNoiseBath("test", [1], [1.0], 1, 50.0, T)
            f, df = 5e9, 1e8
            variance = JC.thermaloccupation(2pi*f, T)/2
            psd = JC.bathamplitude(bath, f, df)^2*variance/(2df)
            expected = T == 0 ? JC.planck_constant*f/50 : 2JC.boltzmann_constant*T/50
            @test psd ≈ expected rtol=1e-7
        end
        # the baths as targets: a resistor bath injects the resistor's own
        # current direction, and its port is zero
        inj = transientinjection(prob, baths)
        @test size(inj) == (length(prob), 2)
        # a port's bath is injected as a port source is, into the port; an
        # internal resistor's into its first terminal, the reverse of a named
        # source's direction
        @test inj[:, 1] == transientinjection(prob, [1])[:, 1]
        @test inj[:, 2] == -transientinjection(prob, ["Rloss"])[:, 1]
        return nothing
    end
    return nothing
end

# the main suite runs the CPU backend by including this file; the GPU
# suite defines TRANSIENTBACKENDTESTS and calls the function with its device
if !@isdefined(TRANSIENTBACKENDTESTS)
    testtransientquantum(JosephsonCircuits.CPU())
    testquantumcontracts()
end
