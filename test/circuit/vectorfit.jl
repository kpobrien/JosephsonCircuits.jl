using JosephsonCircuits
using LinearAlgebra
using SparseArrays
using Random
using Test

# The fit of a rational scattering block to sampled scattering data: the
# pole relocation and its pruning, the residues at settled poles and the
# value stated at zero frequency, the order the samples support, the
# passivity enforcement, and the realization the block is built from.
# What the fitted block then does in time is in transient/system.jl.
@testset "a rational block fitted to scattering data" begin
    JC = JosephsonCircuits
    # The response of a realization by a dense solve, and its largest
    # singular value sampled densely, on a logarithmic grid past its
    # poles and across every resonance: a check of passivity independent
    # of the norm search the enforcement and the validation decide on.
    response(A, B, C, D, w) = D .+ C*((im*w*I - A) \ B)
    function densemax(A, B, C, D)
        poles = eigvals(A)
        mags = filter(>(0), abs.(poles))
        ws = vcat(0.0, exp.(range(log(minimum(mags)/100), log(maximum(mags)*100); length = 4000)))
        for l in poles
            imag(l) > 0 && append!(ws, imag(l) .+ abs(real(l)) .* range(-8, 8; length = 401))
        end
        return maximum(w -> opnorm(response(A, B, C, D, w)), filter(>=(0), ws))
    end
    densemax(p::JC.RationalScatteringProvider) = densemax(p.A, p.B, p.C, p.D)
    # the samples of a two port delay line of `delay` seconds
    function delaydata(gs, delay)
        S = zeros(ComplexF64, 2, 2, length(gs))
        for (k, g) in enumerate(gs)
            S[:, :, k] .= [0 cis(-2pi*g*delay); cis(-2pi*g*delay) 0]
        end
        return S
    end
    # the frequencies and samples of a lossy line two centimetres long at
    # `K` frequencies spread evenly in the logarithm from 1 kHz to 20 GHz:
    # diffusion through its resistance at the low end, skin effect and
    # dielectric loss above, a turn of delay at the top; no rational
    # function, so its fit needs poles over every decade
    function lossyline(K)
        local fl = exp.(range(log(1e3), log(20e9); length = K))
        local Sl = zeros(ComplexF64, 2, 2, K)
        for (k, f) in enumerate(fl)
            local z = 2e4 + (1 + im)*2e-3*sqrt(f) + im*2pi*f*4e-7
            local y = im*2pi*f*1.6e-10*(1 - 2e-3im)
            local zc, g = sqrt(z/y), 0.02*sqrt(z*y)
            local a, b, c = cosh(g), zc*sinh(g), sinh(g)/zc
            local den = 2a + b/50 + 50c
            Sl[:, :, k] .= [(b/50 - 50c)/den 2/den; 2/den (b/50 - 50c)/den]
        end
        return fl, Sl
    end
    # the scattering data of an RLC two-port, tabulated as the
    # linearized solver gives it, fitted by vector fitting at its
    # three poles and at more: the fit is exact, the poles the RLC's,
    # and the extra poles dropped; the raw fit is returned on request,
    # and the sampling is checked
    rlc(R) = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:l1, 1, 2, Inductor(1.5e-9)), (:c, 2, 0, Capacitor(0.6e-12)),
        (:r, 2, 0, R), (:l2, 2, 3, Inductor(1.0e-9)), (:p2, 3, 0, Port(2; Z0 = 50.0))])
    fs = collect(range(0.2e9, 12e9; length = 240))
    hb = hblinsolve(2pi .* fs, rlc(Resistor(120.0)); keyedarrays = false)
    data = ScatteringParameters((2pi .* fs, hb.S); nports = 2, zref = 50.0)
    fitted = RationalScattering(data, 3)
    # the realization is minimal: the RLC's residues have rank one, so
    # its three poles take three states
    @test fitted.provider isa JC.RationalScatteringProvider && size(fitted.provider.A) == (3, 3)
    # the RLC's poles are the finite eigenvalues of its descriptor
    # pencil in the node voltages and the inductor currents, each of
    # them a pole of the fit once per port
    L1, C1, R1, L2, Z = 1.5e-9, 0.6e-12, 120.0, 1.0e-9, 50.0
    E = zeros(5, 5); E[2, 2] = C1; E[4, 4] = L1; E[5, 5] = L2
    M = zeros(5, 5); M[1, 1] = 1/Z; M[1, 4] = 1; M[2, 2] = 1/R1; M[2, 4] = -1; M[2, 5] = 1; M[3, 3] = 1/Z; M[3, 5] = -1
    M[4, 1] = -1; M[4, 2] = 1; M[5, 2] = -1; M[5, 3] = 1
    order(x) = (imag(x), real(x))
    expected = sort(filter(x -> abs(x) < 1e12, eigvals(-M, E)); by = order)
    poles = sort(eigvals(fitted.provider.A); by = order)
    @test poles ≈ sort(expected; by = order) rtol=1e-6
    fit = zeros(ComplexF64, 2, 2, length(fs))
    JC.evaluateprovider!(fit, fitted.provider, 2pi .* fs)
    @test maximum(abs.(fit .- hb.S)) < 1e-10
    @test densemax(fitted.provider) <= 1 + 1e-10
    more = RationalScattering(data, 8; frequencies = 2pi .* fs)
    @test size(more.provider.A) == (3, 3)
    JC.evaluateprovider!(fit, more.provider, 2pi .* fs)
    @test maximum(abs.(fit .- hb.S)) < 1e-10
    raw = RationalScattering(data, 4; passivity = nothing)
    @test size(raw.provider.A) == (3, 3)
    @test densemax(raw.provider) <= 1 + 1e-10
    # its check runs the enforcement's sweep within the enforcement's time,
    # and refuses a fit the sweep has not settled by then: this one, which
    # takes intervals to settle, given none
    let ws = 2pi .* fs
        local p, Rr, Dr = JC.vectorfit(hb.S, ws, 4, VectorFitting())
        local sp = JC.residuespaces(p, Rr, 1e-12)
        @test isnothing(JC.checkrawpassive(p, Rr, Dr, ws, sp; atol = 1e-8))
        @test_throws ArgumentError JC.checkrawpassive(p, Rr, Dr, ws, sp; atol = 1e-8, maxtime = 0.0)
    end
    # and refused where the sweep of its residues finds it above one: a
    # lossy two port's samples raised by a quarter stay under one over the
    # band, and their exact fit stands above it below the band, as the
    # realization shows
    let w4 = 2pi*4e9, gs = collect(range(1e9, 8e9; length = 40)),
            lossy = RationalScattering(-w4 .* Matrix(1.0I, 2, 2), w4 .* Matrix(1.0I, 2, 2),
                0.8 .* [0.0 1.0; 1.0 0.0], zeros(2, 2); zref = 50.0),
            Sl = JC.evaluateprovider!(zeros(ComplexF64, 2, 2, length(gs)), lossy.provider, 2pi .* gs),
            lifted = ScatteringParameters((2pi .* gs, 1.26 .* Sl); nports = 2, zref = 50.0)
        @test densemax(RationalScattering(lifted, 2; passivity = nothing, atol = 0.1).provider) > 1.005
        @test_throws ArgumentError RationalScattering(lifted, 2; passivity = nothing)
        # Weighted, the correction's metric counts against the memory
        # `maxstates` stands for, and a budget of the fit's own two states
        # leaves it none: refused before the metric is formed
        local weighed = [1.0 + 0.1i + 0.2j for i in 1:2, j in 1:2, _ in gs]
        @test_throws ArgumentError RationalScattering(lifted, 2; weights = weighed, maxstates = 2)
        @test RationalScattering(lifted, 2; weights = weighed) isa ScatteringParameters
        # and weights whose squares underflow leave an entry's metric
        # singular, which is refused as a fit that cannot be made passive
        local tiny = copy(weighed)
        tiny[1, :, :] .= 1e-170
        @test_throws ArgumentError RationalScattering(lifted, 2; weights = tiny)
    end
    @test_throws ArgumentError RationalScattering(data, 4; frequencies = reverse(2pi .* fs))
    @test_throws ArgumentError RationalScattering(fitted, 4)
    @test_throws ArgumentError RationalScattering(data, 0)
    # The tolerances of the fit are the caller's, since a block with
    # little loss needs them tighter than the defaults and a hard fit
    # needs its relocation and its sweep for violations finer: the
    # parameters of the fit and of the enforcement are option objects.
    # They are checked, and the exact fit above is reached whatever they
    # are, since it needs no enforcement.
    @test_throws ArgumentError PassivityEnforcement(margin = -1.0)
    @test_throws ArgumentError PassivityEnforcement(margin = Inf)
    @test_throws ArgumentError PassivityEnforcement(rounds = 0)
    @test_throws ArgumentError PassivityEnforcement(scalelimit = 0.5)
    @test_throws ArgumentError PassivityEnforcement(regularization = 0.0)
    @test_throws ArgumentError PassivityEnforcement(maxtime = NaN)
    @test_throws ArgumentError VectorFitting(pruneslack = -0.5)
    @test_throws ArgumentError VectorFitting(iterations = 0)
    # a start is a spacing by name, or poles which can start a fit:
    # finite, stable, and each complex one with its conjugate
    for bad in (:geometric, 3, ComplexF64[], [-1.0, NaN], [0.0, -1.0], [complex(1.0, 1.0), complex(1.0, -1.0)],
            [complex(-1.0, 1.0)], [complex(-1.0, 1.0), complex(-1.0, -1.1)])
        @test_throws ArgumentError VectorFitting(start = bad)
    end
    @test_throws ArgumentError RationalScattering(data, 3; tol = -1.0)
    # A fixed order is held to its data as the search is, against a
    # tolerance of its own: two poles miss this RLC, which has three,
    # and are refused unless a tolerance that loose is asked for; the fit
    # returned is made passive by its corrections alone, with no
    # contraction to warn of
    @test_throws ArgumentError RationalScattering(data, 2)
    short = @test_logs RationalScattering(data, 2; tol = 1.0)
    @test 1e-2 < JC.relativefiterror(short, hb.S, 2pi .* fs) <= 1.0
    @test densemax(short.provider) <= 1 + 2e-8
    for kw in ((; passivity = PassivityEnforcement(margin = 1e-9, rounds = 40, scalelimit = 0.1)),
            (; passivity = PassivityEnforcement(regularization = 1e-10)),
            (; fitting = VectorFitting(pruneslack = 0.0, stallpatience = 10)),
            (; fitting = VectorFitting(pruneslack = 1.0, iterations = 60)))
        tuned = RationalScattering(data, 3; kw...)
        JC.evaluateprovider!(fit, tuned.provider, 2pi .* fs)
        @test maximum(abs.(fit .- hb.S)) < 1e-10
    end
    # a pole is kept unless the fit is as close without it: the
    # permissive threshold of a factor of two would let the error grow
    # by that much at every pole it drops
    @test size(RationalScattering(data, 8; fitting = VectorFitting(pruneslack = 0.0)).provider.A) == (3, 3)
    # A block which states what it does at zero frequency has the fit
    # meet it exactly, which nothing in the data can do when the
    # samples begin above zero. At zero this RLC is a shunt resistor
    # between two ports: its inductors are shorts and its capacitor an
    # open, so `S11 = -y/(2 + y)` and `S21 = 2/(2 + y)` for the
    # resistor's conductance `y` in units of the port's.
    y = 50.0/120.0
    S0 = [(-y/(2 + y)) (2/(2 + y)); (2/(2 + y)) (-y/(2 + y))]
    stated = ScatteringParameters((2pi .* fs, hb.S); nports = 2, zref = 50.0,
        dcmodel = JC.ScatteringDC(S0))
    told = RationalScattering(stated, 6)
    P = told.provider
    reached = real.(response(P.A, P.B, P.C, P.D, 0.0))
    # exactly, and through the passivity enforcement, not only out of
    # the residue solve: the condition is carried into the correction
    @test reached ≈ S0 atol=1e-12
    # and the fit is no worse for it, the statement being the truth
    JC.evaluateprovider!(fit, P, 2pi .* fs)
    @test maximum(abs.(fit .- hb.S)) < 1e-8
    # Contracting toward a statement is a different step from
    # contracting toward nothing, and the step is found on the path
    # rather than derived from a bound: the norm along the path is a
    # norm of something affine in the step, so it is convex in it,
    # and the steps which meet a target are an interval containing
    # zero whose end the search finds. The bounds do not stand in for
    # it. `1/M` answers only contraction toward nothing:
    # `0.5 + 0.501 s/(s + 1)` anchored at `S(0) = 0.5` plainly has a
    # step, while repeated `1/M` steps stall short of one. The
    # triangle bound `t <= (1 - a)/(M - a)`, `a` the statement's
    # norm, is sufficient and not necessary: at `a = 1` its
    # numerator vanishes and it says nothing, though a positive step
    # can exist.
    @test JC.contractionstep(τ -> 1.001*τ, 1.0) ≈ 1/1.001 rtol=1e-9
    @test JC.contractionstep(τ -> abs(1 - 2.001τ), 1.0) ≈ 2/2.001 rtol=1e-9
    @test JC.contractionstep(τ -> 0.5, 1.0) == 1        # nothing to do
    @test JC.contractionstep(τ -> 2.0, 1.0) == 0        # no step helps
    # the resolution of the step is the caller's; a coarse one still
    # returns a step which measures under the target
    @test JC.contractionstep(τ -> 1.001*τ, 1.0; tol = 1e-3) ≈ 1/1.001 atol=1e-3
    # and on a measure which curves, as a contraction toward a statement
    # of nearly unit norm does, the search reaches the boundary, the far
    # end of its bracket moving as well as the near one
    @test JC.contractionstep(τ -> 1 - 3e-8 + 0.0098τ^2, 1 + 5e-9) ≈ sqrt(3.5e-8/0.0098) rtol=1e-6
    # A statement of unit norm, which a through, an open and a short
    # all are, pins the norm of any fit meeting it at one: there is
    # nothing to contract away, and contracting toward anything else
    # would move the statement. Such a fit is accepted at one to the
    # tolerance it is validated against. Where the two cannot both
    # hold, the fit is refused and says so, rather than returning a
    # block which meets neither.
    @test_throws ArgumentError RationalScattering(ScatteringParameters(
        (2pi .* fs, hb.S); nports = 2, zref = 50.0, dcmodel = JC.ThroughDC()), 6)
    # a statement at odds with the poles the data gives is visible at
    # once, in the error it costs in band. That is not the data
    # judging the statement: nothing below the samples is measured,
    # and poles placed below the band would let the same fit meet a
    # value it refuses here. It is only that a fit cannot reach where
    # its poles do not go.
    wrong = ScatteringParameters((2pi .* fs, hb.S); nports = 2, zref = 50.0,
        dcmodel = JC.ThroughDC())
    pw, rw, Dw = JC.vectorfit(hb.S, 2pi .* fs, 6, VectorFitting();
        dc = JC.dcscatteringmatrix(JC.ThroughDC(), 2))
    Aw, Bw, Cw = JC.realization(pw, rw, 2)
    JC.evaluateprovider!(fit, JC.RationalScatteringProvider(Aw, Bw, Cw, Dw), 2pi .* fs)
    @test maximum(abs.(fit .- hb.S)) > 1e-3
    # A sample at zero frequency is data like any other, and the fit
    # is scaled by the geometric centre of the band, which zero has
    # no part in: the lowest positive frequency stands in for it, in
    # the scaling and in the spread the starting poles are given. A
    # pole started at zero would sit on top of the sample there,
    # where the basis it belongs to is infinite, and the failure
    # would surface as infinities inside an eigensolver.
    withdc = vcat(0.0, fs)
    Sdc = zeros(ComplexF64, 2, 2, length(withdc))
    Sdc[:, :, 2:end] .= hb.S
    Sdc[:, :, 1] .= hb.S[:, :, 1]
    pdc, rdc, Ddc = JC.vectorfit(Sdc, 2pi .* withdc, 6, VectorFitting())
    @test all(isfinite, pdc) && all(isfinite, rdc) && all(isfinite, Ddc)
    @test all(p -> real(p) < 0, pdc)
    # and a fit with no positive frequency has no band to be scaled by
    @test_throws ArgumentError JC.vectorfit(Sdc[:, :, 1:2], [0.0, 0.0], 4, VectorFitting())
    # Zero frequency is data, and the public fits take it: a sample
    # there is the block's value at DC, and the coordinates the fit
    # works in are normalized by the lowest positive frequency
    # rather than the lowest, so a band beginning at zero still has
    # a scale.
    Az = reshape([-1.0], 1, 1); Bz = reshape([1.0], 1, 1)
    Cz = reshape([0.5], 1, 1); Dz = reshape([0.0], 1, 1)
    wz = 2pi .* collect(range(0.0, 1.0; length = 20))
    Sz = zeros(ComplexF64, 1, 1, length(wz))
    for (k, w) in enumerate(wz)
        Sz[1, 1, k] = only(response(Az, Bz, Cz, Dz, w))
    end
    blkz = ScatteringParameters((wz, Sz); nports = 1, zref = 50.0)
    for got in (RationalScattering(blkz, 1), RationalScattering(blkz; tol = 1e-6))
        Sf1 = zeros(ComplexF64, 1, 1, length(wz))
        JC.evaluatescattering!(Sf1, got, wz)
        @test maximum(abs, Sf1 .- Sz) < 1e-12
        @test real(Sf1[1, 1, 1]) ≈ 0.5 atol=1e-12
    end
    # a scan from minpoles needs minpoles + 1 samples, and too few
    # are refused by name rather than reported as a scan of no orders
    @test_throws ArgumentError RationalScattering(blkz; tol = 1e-6, minpoles = 30)
    # The samples are checked once, by one function: duplicates are
    # refused, since two samples at one frequency make the divided
    # differences of the Loewner pencil singular and leave a seeded
    # pole pair no width.
    @test JC.checkfrequencies([0.0, 1.0, 2.0]) == [0.0, 1.0, 2.0]
    @test_throws ArgumentError JC.checkfrequencies([0.0, 0.0])
    @test_throws ArgumentError JC.checkfrequencies([1.0, 1.0, 2.0])
    @test_throws ArgumentError JC.checkfrequencies([0.0, Inf])
    for bad in ([0.0, 0.1, 0.1, 0.2], [0.0, 0.0, 0.0, 0.0],
                [-0.1, 0.0, 0.1, 0.2], [0.2, 0.1, 0.3, 0.4])
        @test_throws ArgumentError RationalScattering(blkz, 1; frequencies = 2pi .* bad)
    end
    # `vectorfit` is the fit and not the enforcement. Samples of
    # `1.5 - 1.2/(s + 1)` over a low band have an exact one pole fit
    # whose constant term is 1.5, and the fit returns it: bringing a
    # constant term under one is part of making the fit passive, done
    # for the caller who asks for it rather than to every fit.
    wsc = 2pi .* collect(range(0.01, 0.3; length = 40))
    Sc = zeros(ComplexF64, 1, 1, length(wsc))
    for (k, w) in enumerate(wsc)
        Sc[1, 1, k] = 1.5 - 1.2/(im*w + 1)
    end
    praw, rraw, Draw = JC.vectorfit(Sc, wsc, 1, VectorFitting())
    @test only(Draw) ≈ 1.5 rtol=1e-8
    @test only(praw) ≈ -1 rtol=1e-6
    @test only(rraw) ≈ -1.2 rtol=1e-6
    # and when it is asked for, the trigger and the target are the
    # caller's, not a threshold of the fit's own
    _, _, Dc = JC.vectorfit(Sc, wsc, 1, VectorFitting(); constanttol = 1e-8, constantmargin = 1e-6)
    @test only(Dc) ≈ 1 - 1e-6 rtol=1e-9
    _, _, Dw = JC.vectorfit(Sc, wsc, 1, VectorFitting(); constanttol = 1e-8, constantmargin = 1e-3)
    @test only(Dw) ≈ 1 - 1e-3 rtol=1e-9
    # a constant term under one is left alone whatever is asked
    @test JC.passiveconstant([0.5 0.0; 0.0 0.25]; tol = 1e-8)[2] == false
    # a value stated at zero frequency must be real, a real rational
    # function being real there
    @test_throws ArgumentError JC.vectorfit(Sc, wsc, 1, VectorFitting(); dc = fill(0.5 + 0.1im, 1, 1))
    # More iterations of the pole relocation never return a worse
    # fit, because what it returns is the iterate measured to fit
    # best and not the last one, nor the one whose poles moved least
    # between steps: a relocation which is circling can pass its
    # closest fit on a step where the poles happen to be moving
    # quickly. A delay, which no order fits exactly, shows it.
    gsr = 2pi .* collect(range(0.5e9, 12e9; length = 120))
    Sdr = delaydata(gsr ./ 2pi, 40e-12)
    xsr = gsr ./ sqrt(first(gsr)*last(gsr))
    startr = ComplexF64[]
    for w in range(first(xsr), last(xsr); length = 3)
        push!(startr, complex(-0.01w, w))
        push!(startr, complex(-0.01w, -w))
    end
    errsr = [JC.fiterror(Sdr, xsr, JC.converge(Sdr, xsr, copy(startr), VectorFitting(iterations = it)))
        for it in 1:12]
    @test all(k -> errsr[k + 1] <= errsr[k]*(1 + 1e-12), 1:length(errsr) - 1)
    @test minimum(errsr) == errsr[end]
    # The rounds go on while the fit improves, whatever its poles do: the
    # lossy line at 500 samples and 24 poles spread like them has its best
    # fit at the eighth round for four rounds, while the change from one
    # round's sorted poles to the next's stays at several times their
    # size, and then a third closer by the twentieth.
    let (fl, Sl) = lossyline(500)
        xl = fl ./ sqrt(fl[1]*fl[end])
        startl = JC.spreadpoles(:log, xl[1], xl[end], 24)
        lineerror(fitting) = JC.fiterror(Sl, xl, JC.converge(Sl, xl, copy(startl), fitting))
        @test lineerror(VectorFitting()) <= 0.8*lineerror(VectorFitting(iterations = 9))
    end
    # Weights weigh each entry at each sample. A weight the same at every
    # sample of an entry scales it: the fit so weighted, which relocates
    # on the entries with a least squares for each, is the fit of the
    # scaled samples, which relocates on their singular components.
    let (fl, Sl) = lossyline(200), a = [1.0 30.0; 0.5 2.0]
        local wl = 2pi .* fl
        local pw, Rw, Dw = JC.vectorfit(Sl, wl, 12, VectorFitting(); weights = repeat(a, 1, 1, length(fl)))
        local ps, Rs, Ds = JC.vectorfit(a .* Sl, wl, 12, VectorFitting())
        @test length(pw) == length(ps)
        local responses(p, R, D) = [D[i, j] + sum(R[i, j, q]/(im*w - p[q]) for q in eachindex(p)) for i in 1:2, j in 1:2, w in wl]
        @test a .* responses(pw, Rw, Dw) ≈ responses(ps, Rs, Ds) rtol = 1e-8
        # At fixed poles each residue solve is the least squares of its own
        # measure: the weighted residues fit closer in the weighted rms, the
        # plain ones in the plain rms.
        local W = 1 ./ abs.(Sl)
        local weighted = JC.residualerrors(Sl, wl, ps, JC.fitresidues(Sl, wl, ps; weights = W)...; weights = W)
        local plain = JC.residualerrors(Sl, wl, ps, JC.fitresidues(Sl, wl, ps)...; weights = W)
        @test weighted[2] < plain[2]
        @test JC.residualerrors(Sl, wl, ps, JC.fitresidues(Sl, wl, ps; weights = W)...)[2] >
            JC.residualerrors(Sl, wl, ps, JC.fitresidues(Sl, wl, ps)...)[2]
        # The weights reach the public fit, whose error is the weighted
        # one, and are refused unless one positive weight is given for
        # each entry at each sample.
        local line = ScatteringParameters((wl, Sl); nports = 2, zref = 50.0)
        local weightedfit = RationalScattering(line, 12; weights = W)
        @test JC.relativefiterror(weightedfit, Sl, wl; weights = W) <= 1e-2
        @test_throws ArgumentError RationalScattering(line, 12; weights = W[:, :, 2:end])
        @test_throws ArgumentError RationalScattering(line, 12; weights = -W)
        # Only the weights' ratios weigh: scaled by 1e-200 the fit is the
        # same, and weights the same everywhere weigh nothing, the fit
        # without weights.
        @test JC.relativefiterror(RationalScattering(line, 12; weights = 1e-200 .* W), Sl, wl; weights = W) ≈
            JC.relativefiterror(weightedfit, Sl, wl; weights = W) rtol = 1e-6
        @test RationalScattering(line, 12; weights = fill(3.0, size(W))).provider.A == RationalScattering(line, 12).provider.A
    end
    # The relocation runs on the singular components of the entries,
    # which have the entries' Gram matrix, all the relaxed weight reads,
    # so they relocate the poles as the entries do: the RLC's four
    # entries are three functions, its transmissions being equal.
    let entries = reshape(permutedims(hb.S, (3, 1, 2)), length(fs), 4),
        comps = JC.relocationcomponents(hb.S), xr = 2pi .* fs ./ sqrt(2pi*fs[1]*2pi*fs[end])
        @test size(comps, 2) == 3
        byorder(p) = sort(p; by = x -> (imag(x), real(x)))
        start = ComplexF64[complex(-0.01w, s*w) for w in range(xr[1], xr[end]; length = 3) for s in (1, -1)]
        @test byorder(JC.relocate(comps, xr, start)) ≈ byorder(JC.relocate(entries, xr, start)) rtol = 1e-10
    end
    # A round of the relocation then costs in proportion to the
    # components rather than to the entries: the delay tiled over an
    # eight port and scaled to the same Gram matrix has sixteen times the
    # entries and still one component, and a round allocates what the
    # two port's does, beyond the measurement of its iterate, which reads
    # every entry.
    function roundalloc(S)
        once, twice = VectorFitting(iterations = 1), VectorFitting(iterations = 2)
        JC.converge(S, xsr, copy(startr), once)
        JC.converge(S, xsr, copy(startr), twice)
        JC.fiterror(S, xsr, startr)
        return (@allocated JC.converge(S, xsr, copy(startr), twice)) -
            (@allocated JC.converge(S, xsr, copy(startr), once)) - (@allocated JC.fiterror(S, xsr, startr))
    end
    @test roundalloc(repeat(Sdr ./ 4, 4, 4, 1)) <= 1.1*roundalloc(Sdr)
    # Where the entries are many, the components are found by probes of
    # the stacked entries, at a cost in proportion to how many there are,
    # and relocate the poles as the whole decomposition's do; where every
    # component would be kept, the relocation runs on the entries, which
    # relocate them the same. `p` resonances and a constant over `n` ports
    # are `p + 1` components however many the ports, nine and, over more
    # than a block of probes and at a thousandth of the level, as a weakly
    # coupled block's entries are, 25, the rank being relative to the
    # largest response; with a relative noise of 1e-4 every entry is one.
    function lowrank(n, K; p = 8, noise = 0.0, level = 1.0)
        rng = Random.Xoshiro(5)
        xl = collect(range(0.5, 3.5; length = K))
        Ms = [randn(rng, n, n) for _ in 0:p]
        Sl = zeros(ComplexF64, n, n, K)
        for (k, x) in enumerate(xl), q in 0:p
            a = complex(-0.02q, 2.8q/p)
            t = q == 0 ? 1.0 : 0.05*(1/(im*x - a) + 1/(im*x - conj(a)))
            Sl[:, :, k] .+= t .* Ms[q + 1]
        end
        return xl, level .* (Sl .+ noise .* randn(rng, ComplexF64, size(Sl)))
    end
    let byorder = p -> sort(p; by = x -> (imag(x), real(x)))
        for (p, noise, level, columns) in ((8, 0.0, 1.0, 9), (24, 0.0, 1e-3, 25), (8, 1e-4, 1.0, 144))
            xl, Sl = lowrank(12, 200; p, noise, level)
            start = JC.spreadpoles(:linear, xl[1], xl[end], 12)
            C = JC.relocationcomponents(Sl)
            @test size(C, 2) == columns
            # every singular component of the stacked entries
            Xl = vcat(real.(transpose(reshape(Sl, :, 200))), imag.(transpose(reshape(Sl, :, 200))))
            F = svd(Xl)
            every = complex.(F.U[1:200, :], F.U[201:400, :]) .* transpose(F.S)
            @test byorder(JC.relocate(C, xl, start)) ≈ byorder(JC.relocate(every, xl, start)) rtol = 1e-10
        end
    end
    # and they allocate in proportion to the samples and the entries
    # rather than to their product: nine components of 16 and of 32 ports
    # at 1000 samples
    function componentbytes(n)
        _, Sl = lowrank(n, 1000)
        @test size(JC.relocationcomponents(Sl), 2) == 9
        return @allocated JC.relocationcomponents(Sl)
    end
    @test componentbytes(32) < 2*componentbytes(16)
    # The relocation starts where `start` puts it. Three resonances
    # decades apart, sampled evenly in the logarithm of the frequency:
    # from their own poles, given in rad/s in either order of each pair,
    # the fit is exact with no relocation, and from pairs spread evenly
    # in the logarithm it is nearly so after one round, where the linear
    # spread, one pair in the lowest four decades, is not yet there. A
    # start with a count other than the fit's is refused, and so is one
    # given to the order search, which fits many orders.
    let ws = 2pi .* exp.(range(log(1e2), log(1e8); length = 100)),
        a = [complex(-0.1w0, s*w0) for w0 in 2pi .* [1e3, 1e5, 1e7] for s in (1, -1)]
        Sw = reshape([sum(0.1*abs(imag(p))/(im*w - p) for p in a) for w in ws], 1, 1, :)
        wref = sqrt(ws[1]*ws[end])
        errorfrom(f) = JC.fiterror(Sw, ws ./ wref, JC.vectorfit(Sw, ws, 6, f)[1] ./ wref)
        @test errorfrom(VectorFitting(start = a, iterations = 1)) < 1e-13
        @test errorfrom(VectorFitting(start = reverse(a), iterations = 1)) < 1e-13
        @test errorfrom(VectorFitting(start = :log, iterations = 1)) < 1e-9 < errorfrom(VectorFitting(iterations = 1))
        # half the pairs spread each way, as VFdriver starts, relocate
        # to the resonances as well, and the order search starts every
        # order so
        @test errorfrom(VectorFitting(start = :linlog)) < 1e-13
        @test size(RationalScattering(data; tol = 1e-6, fitting = VectorFitting(start = :linlog)).provider.A) == (3, 3)
        @test_throws ArgumentError JC.vectorfit(Sw, ws, 4, VectorFitting(start = a))
        @test_throws ArgumentError RationalScattering(data; tol = 1e-6, fitting = VectorFitting(start = expected))
    end
    # The order can be searched for instead of given: the fewest poles
    # which hold the error over the samples under a tolerance, as a
    # fraction of the largest response. This RLC is exact at three
    # poles, so every tolerance it can meet it meets there, and the
    # search returns three however loose or tight the tolerance is.
    auto = RationalScattering(data; tol = 1e-6)
    @test size(auto.provider.A) == (3, 3)
    JC.evaluateprovider!(fit, auto.provider, 2pi .* fs)
    @test maximum(abs.(fit .- hb.S)) < 1e-10
    @test JC.relativefiterror(auto, hb.S, 2pi .* fs) <= 1e-6
    # the tolerance bounds what the caller receives: a looser one may
    # be met by fewer poles but never by a fit which misses it
    for tol in (1e-2, 1e-8)
        got = RationalScattering(data; tol = tol)
        @test JC.relativefiterror(got, hb.S, 2pi .* fs) <= tol
        @test size(got.provider.A, 1) <= 3
    end
    # One evaluator for the relocation's choice of iterate, the
    # pruning, and the acceptance: it measures in the spectral norm
    # the public acceptance uses, and under the same zero frequency
    # condition the final residue solve will impose.
    let ps = [complex(-0.3, 1.0), complex(-0.3, -1.0)]
        xq = collect(range(0.4, 2.5; length = 30))
        Sq = zeros(ComplexF64, 2, 2, length(xq))
        for (k, x) in enumerate(xq)
            Sq[:, :, k] .= [0.2 0.9; 0.9 0.2] ./ (1 + im*x)
        end
        res, Dq = JC.fitresidues(Sq, xq, ps)
        deviation(k) = Dq .+ sum(res[:, :, q] ./ (im*xq[k] - ps[q]) for q in eachindex(ps)) .-
            view(Sq, :, :, k)
        byhand = maximum(k -> opnorm(deviation(k)), eachindex(xq))
        @test JC.fiterror(Sq, xq, ps) ≈ byhand
        # and the rms over every entry and sample the pruning holds as
        # well comes from the same residue solve
        @test collect(JC.fiterrors(Sq, xq, ps)) ≈
            [byhand, sqrt(sum(k -> sum(abs2, deviation(k)), eachindex(xq))/length(Sq))]
        # and the condition is carried, so the error reported is the
        # error of the fit that will be built
        resd, Dd = JC.fitresidues(Sq, xq, ps; dc = [0.1 0.8; 0.8 0.1])
        byhandd = maximum(eachindex(xq)) do k
            opnorm(Dd .+ sum(resd[:, :, q] ./ (im*xq[k] - ps[q]) for q in eachindex(ps)) .-
                   view(Sq, :, :, k))
        end
        @test JC.fiterror(Sq, xq, ps; dc = [0.1 0.8; 0.8 0.1]) ≈ byhandd
        @test JC.fiterror(Sq, xq, ps; dc = [0.1 0.8; 0.8 0.1]) > JC.fiterror(Sq, xq, ps)
    end
    # a stated value at zero survives the pruning: the poles that are
    # left still reach it, because every candidate was judged with it
    let xr = 2pi .* fs ./ sqrt(2pi*fs[1]*2pi*fs[end]),
        thru = JC.dcscatteringmatrix(JC.ThroughDC(), 2)
        start = ComplexF64[]
        for w in range(xr[1], xr[end]; length = 5)
            push!(start, complex(-0.01w, w)); push!(start, complex(-0.01w, -w))
        end
        settled = JC.converge(hb.S, xr, start, VectorFitting(); dc = thru)
        kept = JC.prunepoles(hb.S, xr, copy(settled), VectorFitting(); dc = thru)
        @test !isempty(kept)
        resk, Dk = JC.fitresidues(hb.S, xr, kept; dc = thru)
        reached = real.(Dk .+ sum(resk[:, :, q] ./ (0.0 - kept[q]) for q in eachindex(kept)))
        # The value at zero is reached to the roundoff of reading it
        # back, which is the size of the terms which cancel in the sum
        # rather than an absolute figure: where the relocation leaves
        # two poles close together the residues are large and opposite.
        # The pole set depends on roundoff, so the tolerance is measured
        # from the set in hand.
        cancellation = maximum(abs, Dk) +
            sum(opnorm(view(resk, :, :, q))/abs(kept[q]) for q in eachindex(kept))
        @test reached ≈ thru atol = 1e-10 + 1e-12*cancellation
    end
    # Two poles which have coalesced give the residue basis two columns
    # the samples cannot tell apart, and the solve leaves that direction
    # out rather than meeting it with large cancelling residues: the
    # residues stay the size of the response and the value stated at
    # zero is reached.
    let xr = 2pi .* fs ./ sqrt(2pi*fs[1]*2pi*fs[end]),
        thru = JC.dcscatteringmatrix(JC.ThroughDC(), 2),
        coalesced = ComplexF64[-4.17712, -4.17712,
            complex(-2.93039, 5.12498), complex(-2.93039, -5.12498)]
        resc, Dc = JC.fitresidues(hb.S, xr, coalesced; dc = thru)
        @test maximum(abs, resc) < 100
        reached = real.(Dc .+ sum(resc[:, :, q] ./ (0.0 - coalesced[q])
            for q in eachindex(coalesced)))
        @test reached ≈ thru atol = 1e-12
        # and the fit is the one the distinct poles give, the
        # repeated pole adding nothing rather than corrupting it
        distinct = ComplexF64[-4.17712, complex(-2.93039, 5.12498),
            complex(-2.93039, -5.12498)]
        @test JC.fiterror(hb.S, xr, coalesced; dc = thru) ≈
            JC.fiterror(hb.S, xr, distinct; dc = thru) rtol = 1e-8
    end
    # Poles 1e-10 apart take coefficients of 1e10 which cancel to the
    # error, so the tiles the samples are screened by round differently
    # from each sample's own residual by as much as the error: the
    # screened largest norm is still the largest of every sample's.
    let xe = collect(range(0.01, 3.0; length = 129)), pe = ComplexF64[-1, -1 - 1e-10]
        Se = repeat(reshape(1 ./ (1 .+ im .* xe).^2, 1, 1, :), 2, 2, 1)
        X, M = JC.fitcoefficients(Se, xe, pe)
        K = length(xe)
        er, ei = zeros(4), zeros(4)
        every = maximum(1:K) do k
            mul!(er, X, view(M, k, :)); mul!(ei, X, view(M, K + k, :))
            opnorm(reshape(complex.(er, ei) .- vec(Se[:, :, k]), 2, 2))
        end
        @test JC.fiterror(Se, xe, pe) ≈ every rtol = 1e-12
    end
    # and the window is real, not hypothetical: a pure delay is
    # passive and irrational, so no order fits it exactly and the
    # error is not monotone in the order, with orders which fit
    # better than either of their neighbours. A tolerance only such
    # an order meets is found only by a scan. The orders and the
    # tolerances here are measured rather than written down, so this
    # does not depend on where roundoff puts the noise floor.
    gs = collect(range(0.5e9, 12e9; length = 200))
    Sdelay = delaydata(gs, 40e-12)
    delayed = ScatteringParameters((2pi .* gs, Sdelay); nports = 2, zref = 50.0)
    # `pruneslack` is a budget for the pruning and not for each
    # deletion: measured per deletion, a run of them could each
    # spend the whole slack and the fit drift as far from where it
    # started as the number of deletions allowed.
    let xd = 2pi .* gs ./ sqrt(2pi*gs[1]*2pi*gs[end])
        start = ComplexF64[]
        for w in range(xd[1], xd[end]; length = 8)
            push!(start, complex(-0.01w, w)); push!(start, complex(-0.01w, -w))
        end
        settled = JC.converge(Sdelay, xd, start, VectorFitting())
        base = JC.fiterror(Sdelay, xd, settled)
        for slack in (0.0, 0.05, 0.5)
            kept = JC.prunepoles(Sdelay, xd, copy(settled), VectorFitting(pruneslack = slack))
            @test length(kept) <= length(settled)
            @test JC.fiterror(Sdelay, xd, kept) <= (1 + slack)*base + JC.fitroundoff*JC.largestopnorm(Sdelay)
        end
    end
    # Each pass of the pruning solves the residues once and scores every
    # deletion from that solve, refitting only a candidate the score
    # cannot rule out: at the data's own order, where every pole is needed
    # and every one is tried, it allocates in proportion to the poles
    # rather than to their square.
    function prunebytes(npairs)
        xp = collect(range(0.5, 3.5; length = 200))
        pp = ComplexF64[]
        for w in range(1.0, 3.0; length = npairs)
            push!(pp, complex(-0.02w, w), complex(-0.02w, -w))
        end
        Sp = reshape([sum(0.04*abs(imag(a))/(im*x - a) for a in pp) for x in xp], 1, 1, :)
        @test length(JC.prunepoles(Sp, xp, pp, VectorFitting())) == 2npairs
        return @allocated JC.prunepoles(Sp, xp, pp, VectorFitting())
    end
    @test prunebytes(32) < 3*prunebytes(16)
    # The pruning holds the rms error as well as the largest: the largest
    # can be pinned by one feature no pole set follows, here a corrupted
    # sample, and would let a pole the data needs go unseen. A strong
    # resonance and a weak broad one, with one sample off by more than
    # the weak one's height: dropping the weak one leaves the largest
    # error at the corrupted sample and more than doubles the rms, and
    # it is kept.
    let xw = collect(range(0.5, 2.0; length = 100)), weak = complex(-0.2, 1.4)
        Sw = zeros(ComplexF64, 1, 1, length(xw))
        for (k, x) in enumerate(xw)
            Sw[1, 1, k] = 0.2 + 0.05/(im*x - complex(-0.05, 0.8)) + 0.05/(im*x - complex(-0.05, -0.8)) +
                0.006/(im*x - weak) + 0.006/(im*x - conj(weak))
        end
        Sw[1, 1, 30] += 0.05
        kept, _, _ = JC.vectorfit(Sw, xw, 4, VectorFitting())
        @test length(kept) == 4 && minimum(abs.(kept .- weak)) < 0.05*abs(weak)
    end
    reach = np -> try
        JC.relativefiterror(RationalScattering(delayed, np; tol = 1e3), Sdelay, 2pi .* gs)
    catch e
        e isa ArgumentError ? Inf : rethrow()
    end
    # as above, the warnings of the orders walked through are captured
    # rather than asserted
    errs = @test_logs match_mode = :any [reach(np) for np in 1:16]
    window = [np for np in 2:15 if errs[np] < min(errs[np-1], errs[np+1])]
    @test !isempty(window)
    # The search fits every order from `minpoles` up and returns the
    # first which meets the tolerance, so whatever error some order
    # achieves, the search asked for that error meets it, with no more
    # states than that order takes. This is the guarantee, and it needs
    # the scan: a tolerance between what the first order of the window
    # reaches and what the better of its neighbours reaches is met by
    # that order alone. Whether it warns of a contraction on the way to
    # its fit depends on roundoff, so the logs are captured.
    let np = first(window), tol = sqrt(errs[np]*min(errs[np-1], errs[np+1]))
        got = @test_logs match_mode = :any RationalScattering(delayed; tol = tol,
            minpoles = np - 1, maxpoles = np + 1)
        @test JC.relativefiterror(got, Sdelay, 2pi .* gs) <= tol
        @test size(got.provider.A, 1) <= 2np
    end
    # minpoles is a floor the search does not fit below: this delay is
    # met to 1e-3 at three poles, and from five the search returns five
    # of them, ten states
    @test errs[3] <= 1e-3 && errs[5] <= 1e-3
    @test size((@test_logs match_mode = :any RationalScattering(delayed; tol = 1e-3,
        minpoles = 5)).provider.A, 1) >= 10
    # The degree the samples determine budgets the search rather than
    # walling it. It is the numerical rank of a pencil built along
    # cycling directions from at most four hundred samples, with the
    # constant term contributing to it, so it can be short of what a
    # block needs; where the scan reaches it without meeting the
    # tolerance and the error is still falling, the search goes on.
    # This delay needs six poles and is estimated at four once the
    # noise floor is put high enough, and the search finds the six. The
    # orders it discards on the way are contracted to be made passive,
    # and only the fit returned would warn of its contraction; this one
    # has none.
    @test JC.supporteddegree(Sdelay, 2pi .* gs, 1e-2) == 4
    expanded = @test_logs RationalScattering(delayed; tol = 1e-8, noisefloor = 1e-2)
    @test size(expanded.provider.A, 1) ÷ 2 > 4
    @test JC.relativefiterror(expanded, Sdelay, 2pi .* gs) <= 1e-8
    # a `maxpoles` given by the caller is a wall, because the caller
    # made it one; a search which returns nothing warns of nothing
    @test_logs @test_throws ArgumentError RationalScattering(
        delayed; tol = 1e-8, noisefloor = 1e-2, maxpoles = 4)
    # and the expansion stops rather than running to the sample count:
    # a tolerance nothing reaches is still reported
    @test_logs @test_throws ArgumentError RationalScattering(
        delayed; tol = 1e-16, maxpoles = 8)
    # Without `maxpoles` the search stops once the order has doubled since
    # its errors last fell by a tenth, and a fall of the rms counts as well
    # as one of the largest deviation: sixteen resonances of equal strength
    # between two ports, their residues reflections turned by 1.3 rad from
    # one to the next, hold the largest deviation at one the fit has yet to
    # take from 9 poles to 23, while the rms falls a tenth by 19, and the
    # search reaches the 32 poles they need.
    let fr = range(1.5e9, 9e9; length = 16)
        a = [complex(-0.03*2pi*f, 2pi*f*sqrt(1 - 0.03^2)) for f in fr]
        Ms = [[cos(1.3k) sin(1.3k); sin(1.3k) -cos(1.3k)] for k in eachindex(fr)]
        equalresponse(w) = [0.05 0.02; 0.02 -0.04] .+
            sum((-2real(p)*im*w/((im*w - p)*(im*w - conj(p)))) .* M for (p, M) in zip(a, Ms))
        peak = maximum(w -> opnorm(equalresponse(w)), 2pi .* range(1e8, 2e10; length = 20000))
        we = 2pi .* collect(range(1e9, 10e9; length = 200))
        Se = stack(0.95/peak .* equalresponse(w) for w in we)
        found = RationalScattering(ScatteringParameters((we, Se); nports = 2, zref = 50.0); tol = 1e-6)
        @test size(found.provider.A, 1) == 64
        @test JC.relativefiterror(found, Se, we) <= 1e-6
    end
    # The degree the samples determine is a property of the data and
    # not of the unit its frequencies are written in. The two halves
    # of the Loewner pencil do not carry the same units, a divided
    # difference of the response being an inverse frequency and a
    # shifted one dimensionless, so a threshold on the pencil built
    # at the frequencies as they come counts a different number as
    # the unit changes; normalized to the band, the count is the
    # same in every unit.
    @test allequal(JC.supporteddegree(Sdelay, (2pi .* gs) .* c, 1e-10)
                   for c in (1e12, 1e6, 1.0, 1e-6, 1e-12))
    @test allequal(JC.supporteddegree(Sdelay, (2pi .* gs) .* c, 1e-12)
                   for c in (1e9, 1.0, 1e-9))
    @test_throws ArgumentError RationalScattering(data; tol = 1e-16, maxpoles = 8)
    @test_throws ArgumentError RationalScattering(data; tol = 0.0)
    @test_throws ArgumentError RationalScattering(data; tol = Inf)
    @test_throws ArgumentError RationalScattering(data; tol = 1e-3, minpoles = 0)
    @test_throws ArgumentError RationalScattering(data; tol = 1e-3, noisefloor = 0.0)
    @test_throws ArgumentError RationalScattering(data; tol = 1e-3, minpoles = 9, maxpoles = 4)
    # The ceiling on the search comes from the data when it is not
    # given: the rank of the Loewner pencil is the degree the samples
    # determine, which bounds the poles because each takes at least
    # one state. This RLC has three poles of rank one residues, so its
    # degree is at least three and the ceiling admits its fit. A
    # noisefloor so loose that it counts nothing still leaves a floor
    # of one rather than a ceiling of none.
    @test JC.supporteddegree(hb.S, 2pi .* fs, 1e-12) >= 3
    @test JC.supporteddegree(hb.S, 2pi .* fs, 0.5) >= 1
    # the rank is a property of the block, not of how finely it was
    # sampled, so thinning the samples does not change it
    @test JC.supporteddegree(hb.S[:, :, 1:2:end], 2pi .* fs[1:2:end], 1e-12) ==
          JC.supporteddegree(hb.S, 2pi .* fs, 1e-12)
    # The samples bound the order from below: a fit of N poles and a
    # constant combines N + 1 real functions of frequency, so its error is
    # at least what the samples leave beyond their best approximation of
    # that rank. Samples of a four port with six poles in three pairs and
    # residues of full rank have rank seven, and six is the fewest poles
    # any fit to roundoff can have.
    let n = 4, ps = ComplexF64[-0.1 + 1im, -0.1 - 1im, -0.05 + 2.5im, -0.05 - 2.5im, -0.2 + 4im, -0.2 - 4im]
        Rs = [reshape(sin.((1:n^2) .* (k + 1)) .+ im .* cos.((1:n^2) .* (2k + 3)), n, n) for k in 1:3]
        D0 = reshape(cos.(1:n^2), n, n)./10
        wx = collect(range(0.1, 6.0; length = 80))
        Sx = [D0[i, j] + sum(Rs[k][i, j]/(im*w - ps[2k - 1]) + conj(Rs[k][i, j])/(im*w - ps[2k]) for k in 1:3)
            for i in 1:n, j in 1:n, w in wx]
        @test JC.fewestpoles(JC.relocationcomponents(Sx), Sx, 1e-10) == 6
        # Four poles miss the third pair, and no fit with them comes
        # closer than the bound of the weighted least squares on the
        # components: not their least squares fit on the samples, nor
        # one with other residues.
        p4, r4, d4 = JC.vectorfit(Sx, wx, 4, VectorFitting())
        bound = JC.deviationfloor(JC.relocationcomponents(Sx), wx, p4, n, Inf)
        worst(R, D) = maximum(k -> opnorm(D .+ sum(R[:, :, p]./(im*wx[k] - p4[p]) for p in eachindex(p4)) .- Sx[:, :, k]), eachindex(wx))
        @test 0 < bound <= JC.fiterrors(Sx, wx, p4)[1]
        @test bound <= worst(r4 .+ 0.1 .* randn(ComplexF64, size(r4)), d4 .+ 0.1 .* randn(n, n))
    end
    # a tolerance the data cannot support names the degree it carries
    cannot = try
        RationalScattering(data; tol = 1e-16); ""
    catch e
        sprint(showerror, e)
    end
    @test occursin("determine a degree of about", cannot)
    # the refusal of a fit or a search as its message, empty where it fits
    message(args...; kw...) = try
        RationalScattering(args...; kw...); ""
    catch e
        e isa ArgumentError || rethrow()
        sprint(showerror, e)
    end
    # The search goes on past orders which miss the tolerance: six
    # coupled resonators between two ports, whose samples determine a
    # degree of about twelve and whose lower orders all miss, are fitted
    # at thirteen, with thirteen states; a maxstates of ten ends the scan
    # at the first order whose fit has more, and says so.
    let comps = Any[(:p1, 1, 0, Port(1))]
        for r in 1:6
            C0 = 1/((2pi*5e9*(1 + 0.02*(r - 3.5)))^2*1e-9)
            append!(comps, [(Symbol(:L, r), r, 0, Inductor(1e-9)), (Symbol(:C, r), r, 0, Capacitor(C0)),
                (Symbol(:R, r), r, 0, Resistor(2e4)), (Symbol(:Cc, r), r, r + 1, Capacitor(0.08e-12))])
        end
        push!(comps, (:p2, 7, 0, Port(2)))
        fc = collect(range(4e9, 6e9; length = 120))
        Sc = hblinsolve(2pi .* fc, Circuit(comps); keyedarrays = false).S
        chain = ScatteringParameters((2pi .* fc, Sc); nports = 2, zref = 50.0)
        @test JC.relativefiterror(RationalScattering(chain; tol = 1e-6), Sc, 2pi .* fc) <= 1e-6
        @test occursin(r"The scan ended at \d+ poles, whose fit has \d+ states, more than maxstates = 10",
            message(chain; tol = 1e-6, maxstates = 10))
        # The relocation ends when its error has settled, `stallpatience`
        # rounds in a row lowering it by less than a thousandth of
        # itself. The chain at 200 samples with a noise of 1e-3, at 18
        # poles, settles in three rounds, and the next five change its
        # error by 4e-4 of itself at most: the relocation returns the best
        # of those eight rounds, as it does with the rule off and the
        # rounds cut to seven. A patience longer than the rounds carries it
        # on through the slow descent after, to a closer fit.
        fn = collect(range(4e9, 6e9; length = 200))
        Sn = hblinsolve(2pi .* fn, Circuit(comps); keyedarrays = false).S
        k = reshape(1:length(Sn), size(Sn))
        Sn .+= 1e-3 .* complex.(sin.(1.3 .* k .^ 2), cos.(0.7 .* k .^ 2 .+ 1))
        xn = 2pi .* fn ./ sqrt(2pi*fn[1]*2pi*fn[end])
        startn = JC.spreadpoles(:linear, first(xn), last(xn), 18)
        settled = JC.converge(Sn, xn, copy(startn), VectorFitting())
        @test settled == JC.converge(Sn, xn, copy(startn), VectorFitting(stallpatience = 30, iterations = 7))
        @test JC.fiterror(Sn, xn, JC.converge(Sn, xn, copy(startn), VectorFitting(stallpatience = 30))) <=
            (1 - 5e-3)*JC.fiterror(Sn, xn, settled)
    end
    # A search which misses names how near it came, and why the orders
    # which could not be fitted failed, since that is a different
    # problem from a tolerance too tight: the RLC sampled at eight
    # frequencies fits at every order up to seven and at none above,
    # which would need more samples than there are. At a tolerance of
    # 1e-16 no order comes near enough to be made passive, and the
    # nearest is named by the bound which refused it. One in which no
    # order produced a fit says that, and has no closest fit to name.
    # Below the degree the samples determine, an order which produces no
    # fit has too few poles for the data and the search goes on past it;
    # past the degree four such orders in a row end it: a constant six
    # port, of degree six, which no order fits since a constant
    # reproduces it, is searched from four poles to ten.
    missed = message(data; tol = 1e-16, minpoles = 4, maxpoles = 9, frequencies = 2pi .* fs[1:30:end])
    @test occursin("misses the data by at least", missed) && occursin("could not be fitted at all", missed)
    # A ripple of 1e-4 on the RLC puts a tolerance of 1e-6 out of reach
    # and the degree the samples determine at the pencil's size, so the
    # scan to twice it would fit every order to 240; it ends instead once
    # the order has doubled since the fit's own error last fell by a
    # tenth, and says so, and with a maxpoles it scans to that. Every
    # order's fit keeps three poles, six states, so a maxstates of five
    # ends the scan where it starts, and refuses a fixed order.
    let ripple = [1e-4*cis(k^2 + 3i + 5j) for i in 1:2, j in 1:2, k in eachindex(fs)]
        noisy = ScatteringParameters((2pi .* fs, hb.S .+ ripple); nports = 2, zref = 50.0, atol = 1e-3)
        @test occursin("The scan ended at", message(noisy; tol = 1e-6))
        walled = message(noisy; tol = 1e-6, maxpoles = 12)
        @test occursin("between 4 and 12 poles", walled) && !occursin("The scan ended", walled)
        @test occursin("where the scan starts, has 6 states", message(noisy; tol = 1e-6, maxpoles = 12, maxstates = 5))
        @test occursin("more than maxstates = 5", message(noisy, 6; tol = 1.0, maxstates = 5))
        @test isempty(message(noisy, 6; tol = 1.0, maxstates = 6))
    end
    let fc = collect(range(0.1, 1.0; length = 21))
        @test JC.supporteddegree(repeat(Matrix{ComplexF64}(0.3I, 6, 6), 1, 1, 21), 2pi .* fc, 1e-12) == 6
        none = message(ScatteringParameters(Matrix(0.3I, 6, 6)); tol = 1e-6, frequencies = 2pi .* fc)
        @test occursin("between 4 and 10 poles could be fitted at all", none) && !occursin("closest", none)
    end
    # A block may be active, at zero frequency as at any other, so long
    # as it declares its own noise; without one an active statement is
    # refused by the block itself, before any fitting, and an active
    # realization by the block's constructor.
    @test_throws ArgumentError ScatteringParameters((2pi .* fs, hb.S); nports = 2,
        zref = 50.0, dcmodel = JC.ScatteringDC([1.4 0.0; 0.0 1.4]))
    @test_throws ArgumentError RationalScattering(zeros(0, 0), zeros(0, 2), zeros(2, 0), [0.0 0.0; 3.0 0.0]; zref = 50.0)
    @test RationalScattering(zeros(0, 0), zeros(0, 2), zeros(2, 0), [0.0 0.0; 3.0 0.0]; zref = 50.0,
        noise = JC.NoiseCovariance([0.5 0.0; 0.0 4.0])) isa ScatteringParameters
    # An amplifier whose gain rolls off, stating the noise a quantum
    # limited one has, |I - S S'|/2, and its active value at zero
    # frequency, is fitted as it is, without the passivity enforcement,
    # which could not meet the statement, and without the validation:
    # the fit is stable, meets the statement exactly, has the gain of
    # the samples, and carries the stated noise, whichever way it is
    # asked for.
    w0, g0 = 2pi*8e9, 10.0
    amp(w) = [0.0 0.0; g0*w0/(w0 + im*w) 0.0]
    ampnoise(w) = (K = I - amp(w)*amp(w)'; [abs(K[1, 1])/2 0.0; 0.0 abs(K[2, 2])/2])
    ampdata = ScatteringParameters(amp; nports = 2, zref = 50.0, noise = JC.NoiseCovariance(ampnoise),
        dcmodel = JC.ScatteringDC([0.0 0.0; g0 0.0]))
    ampws = 2pi .* collect(range(0.1e9, 40e9; length = 300))
    for ampfit in (RationalScattering(ampdata, 4; frequencies = ampws),
            RationalScattering(ampdata; tol = 1e-6, minpoles = 1, maxpoles = 6, frequencies = ampws),
            RationalScattering(ampdata, 4; frequencies = ampws, passivity = nothing))
        @test ampfit.noise isa JC.NoiseCovariance
        @test size(ampfit.provider.A) == (1, 1) && maximum(real.(eigvals(ampfit.provider.A))) < 0
        Sa = zeros(ComplexF64, 2, 2, 4)
        JC.evaluateprovider!(Sa, ampfit.provider, 2pi .* [0.0, 1e9, 5e9, 20e9])
        @test Sa ≈ cat([amp(2pi*f) for f in (0.0, 1e9, 5e9, 20e9)]...; dims = 3) rtol=1e-8
        @test abs(Sa[2, 1, 2]) > 1
        @test_throws ArgumentError JC.checkpassive(ampfit.provider)
    end
    # The constant term is the model at infinite frequency. One above one
    # is brought under one before the residues are fitted, and the
    # residues then fitted to what is left over the samples, where the
    # enforcement would contract the whole fit, or refuse it beyond its
    # `scalelimit`.
    @test JC.passiveconstant([0.5 0.0; 0.0 0.25])[2] == false
    @test JC.passiveconstant([0.5 0.0; 0.0 0.25])[1] == [0.5 0.0; 0.0 0.25]
    # a block which is lossless at infinite frequency has a unitary
    # constant term by right, and moving it would put an error into a
    # fit which was exact
    @test JC.passiveconstant([0.0 1.0; 1.0 0.0])[2] == false
    @test JC.passiveconstant(Matrix(1.0I, 3, 3))[2] == false
    # one genuinely above unity is brought down, its singular vectors
    # left alone
    let (Dc, hit) = JC.passiveconstant([4.0 0.0; 0.0 0.5])
        @test hit
        @test opnorm(Dc) <= 1
        @test svdvals(Dc) ≈ [1 - 1e-6, 0.5] rtol=1e-9
        @test svd(Dc).U ≈ svd([4.0 0.0; 0.0 0.5]).U rtol=1e-9
    end
    # the residue fit takes a fixed constant, fitting only the
    # strictly proper part to what is left
    let ws2 = 2pi .* collect(range(1e9, 5e9; length = 20)),
        pol = ComplexF64[-1e10 + 3e10im, -1e10 - 3e10im],
        Sd = zeros(ComplexF64, 1, 1, 20)
        for (k, w) in enumerate(ws2)
            Sd[1,1,k] = 0.3 + 2e10/(im*w - pol[1]) + 2e10/(im*w - pol[2])
        end
        r1, d1 = JC.fitresidues(Sd, ws2, pol)
        @test d1[1,1] ≈ 0.3 rtol=1e-8
        r2, d2 = JC.fitresidues(Sd, ws2, pol; constant = fill(0.1, 1, 1))
        @test d2[1,1] == 0.1
        # the strictly proper part absorbs the difference: refitted
        # around the held constant, it is closer to the data than the
        # residues fitted with the constant free, which miss by the
        # whole difference once the constant is swapped
        model(r, d, w) = d + sum(r[1,1,q]/(im*w - pol[q]) for q in 1:2)
        refit = maximum(abs(model(r2, 0.1, w) - Sd[1,1,k]) for (k, w) in enumerate(ws2))
        swapped = maximum(abs(model(r1, 0.1, w) - Sd[1,1,k]) for (k, w) in enumerate(ws2))
        @test refit < swapped
    end

    # Every fit which is returned is passive over the whole imaginary
    # axis and not merely at the samples: the enforcement accepts a fit
    # only where the certified sweep settles every frequency under the
    # level, or where the ceiling of the norm search, an upper bound on
    # its largest singular value over the axis, stands under it. A block
    # whose largest singular value exceeds one is active, and a transient
    # built on it grows without bound, so a nearly passive fit is
    # contracted the rest of the way rather than returned as it is.
    # The near lossless case is the one which needs it: there the
    # enforcement's own level of `sqrt(1 + atol)` sits above the block's
    # entire dissipation. Whether a given fit is contracted, which
    # warns, or refused, which throws into the catch below, turns on
    # the roundoff of the norm search, so the logs are captured
    # without requiring a warning.
    @test_logs match_mode = :any for (rr, nps) in (
        (Resistor(120.0), (2, 4, 8)), (Resistor(1e9), (2, 4, 8)))
        lossless = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:l1, 1, 2, Inductor(2e-9)),
            (:c, 2, 0, Capacitor(0.3e-12)), (:r, 2, 0, rr), (:l2, 2, 3, Inductor(2e-9)),
            (:p2, 3, 0, Port(2; Z0 = 50.0))])
        gs = collect(range(1e9, 12e9; length = 60))
        hbl = hblinsolve(2pi .* gs, lossless; keyedarrays = false)
        dat = ScatteringParameters((2pi .* gs, hbl.S); nports = 2, zref = 50.0)
        for np in nps
            # two poles miss this data, which is the fit that needs the
            # enforcement most, so any fit is accepted
            f = try
                RationalScattering(dat, np; tol = 1.0)
            catch e
                # a fit far from passive is refused rather than scaled
                # into a block which transmits nothing
                @test e isa ArgumentError
                continue
            end
            q = f.provider
            @test densemax(q) <= 1 + 2e-8
        end
    end
    # a network which is a perfect open at one port and a perfect
    # short at the other at infinite frequency fits with its
    # feedthrough exactly on the unit circle
    embed = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)),
        (:le, 1, 2, Inductor(1e-9)), (:re, 2, 3, Resistor(0.5)),
        (:ce, 3, 0, Capacitor(1e-12)), (:p2, 3, 0, Port(2; Z0 = 50.0))])
    ghz = collect(range(0.5e9, 10e9; length = 80))
    hbe = hblinsolve(2pi .* ghz, embed; keyedarrays = false)
    fite = @test_logs match_mode = :any RationalScattering(
        ScatteringParameters((2pi .* ghz, hbe.S); nports = 2, zref = 50.0), 2)
    @test maximum(abs.(abs.(diag(fite.provider.D)) .- 1)) < 1e-9
    # a fit needs its last pole unless a constant reproduces the data:
    # one real pole and one conjugate pair fitted at their order and
    # above keep it, and constant data is refused as a rational block
    onepole = RationalScattering(fill(-1.0, 1, 1), ones(1, 1), fill(0.5, 1, 1), zeros(1, 1))
    for np in (1, 2, 3, 4)
        f1 = RationalScattering(onepole, np; frequencies = 2pi .* collect(range(0.01, 1.0; length = 41)))
        @test size(f1.provider.A) == (1, 1)
        @test f1.provider.A[1, 1] ≈ -1.0 rtol=1e-8
    end
    pair = RationalScattering([-0.1 1.0; -1.0 -0.1], [0.0; 1.0;;], [0.05 0.0], zeros(1, 1))
    for np in (2, 3, 5)
        f2 = RationalScattering(pair, np; frequencies = 2pi .* collect(range(0.05, 3.0; length = 121)))
        @test size(f2.provider.A) == (2, 2)
        @test sort(eigvals(f2.provider.A); by = imag) ≈ [-0.1 - im, -0.1 + im] rtol=1e-8
    end
    @test_throws ArgumentError RationalScattering(ScatteringParameters(fill(0.3, 1, 1)), 2; frequencies = 2pi .* collect(range(0.1, 1.0; length = 21)))
    # A strictly proper fit, as the parts of a pumped block's conversions
    # are, holds its constant at zero in the relocation and the pruning as
    # well as in its residues, so that its poles are chosen for the fit it
    # returns: data which vanishes at infinite frequency, sampled over a
    # band that ends before it decays, with resonances in it, two poles
    # above it and a delay turning its phase by most of a radian over it,
    # fits at eight poles to within 1e-7 of its largest value.
    let xs = collect(range(0.3, 3.0; length = 200)),
            G = reshape([0.1*(0.2/(im*x - complex(-0.05, 1.3)) + 0.2/(im*x - complex(-0.05, -1.3)) +
                0.1/(im*x - complex(-0.1, 0.6)) + 0.1/(im*x - complex(-0.1, -0.6)) + 12/(im*x + 30) +
                36/(im*x + 90))*cis(-0.3x) for x in xs], 1, 1, :),
            (pp, Rp, Dp) = JC.vectorfit(G, xs, 8, VectorFitting(); proper = true, lastpole = :keep)
        @test iszero(Dp)
        @test maximum(k -> abs(G[1, 1, k] - sum(Rp[1, 1, p]/(im*xs[k] - pp[p]) for p in eachindex(pp))),
            eachindex(xs)) <= 1e-7*maximum(abs, G)
    end
    # A pole pair near the real axis, `-1 ± 1.1e-6 i`, whose two functions
    # of the real basis differ by the pair's small imaginary part: formed
    # without the difference of the pole's term and its conjugate's, as the
    # sweep forms them, the residue solve at the exact poles reproduces
    # `0.1 + 0.8/((1 + s)^2 + b^2)` to the roundoff, as its realization
    # evaluates it
    let b = 1.1e-6, xs = collect(range(0.01, 10.0; length = 400)),
            Sb = reshape([complex(0.1 + 0.8/((1 + im*x)^2 + b^2)) for x in xs], 1, 1, :),
            pb = ComplexF64[complex(-1, b), complex(-1, -b)], (Rb, Db) = JC.fitresidues(Sb, xs, pb),
            qb = JC.RationalScatteringProvider(JC.realization(pb, Rb, 1)..., Db)
        @test maximum(abs, JC.evaluateprovider!(zeros(ComplexF64, 1, 1, length(xs)), qb, xs) .- Sb) <= 1e-14
    end
    # a three port, whose nine entries are eliminated one by one in the
    # fit, so nothing the size of every entry's every sample is formed
    tee = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:r1, 1, 4, Resistor(20.0)), (:p2, 2, 0, Port(2; Z0 = 50.0)), (:r2, 2, 4, Resistor(20.0)),
        (:p3, 3, 0, Port(3; Z0 = 50.0)), (:l3, 3, 5, Inductor(2e-9)), (:c3, 5, 4, Capacitor(0.4e-12)), (:r4, 4, 0, Resistor(80.0)), (:c4, 4, 0, Capacitor(0.2e-12))])
    hb3 = hblinsolve(2pi .* fs, tee; keyedarrays = false)
    fit3 = RationalScattering(ScatteringParameters((2pi .* fs, hb3.S); nports = 3, zref = 50.0), 8)
    S3 = zeros(ComplexF64, 3, 3, length(fs))
    JC.evaluateprovider!(S3, fit3.provider, 2pi .* fs)
    @test maximum(abs.(S3 .- hb3.S)) < 1e-10 && size(fit3.provider.A, 1) <= 6
    @test_throws ArgumentError RationalScattering(data, 4; frequencies = 2pi .* fs[1:4])
end

# A delay is not a rational function, so a cable fitted whole spends its
# poles on the phase of the delay and still answers throughout it. Taken
# out before the fit and put back as a line in cascade, the delay is
# exact and only what is left of the cable is fitted.
@testset "a delay taken out before the fit" begin
    JC = JosephsonCircuits
    tau, wc, len = 1e-9, 2pi*8e9, 0.15
    vp = 2*len/tau                      # a line of tau/2 at each port
    cable(w) = (g = exp(-im*w*tau)/(1 + im*w/wc); ComplexF64[0 g; g 0])
    block = ScatteringParameters(cable; nports = 2, zref = 50.0)
    fs = collect(range(0.0, 10e9; length = 200))
    taus = (tau/2, tau/2)
    fitted = RationalScattering(block, 4; frequencies = 2pi .* fs, delays = taus)
    # what is left of the cable is the one pole it really is
    @test size(fitted.provider.A, 1) == 2
    # and putting the delay back reproduces the cable: the entry from
    # port q to port p carries the delay of both
    F = zeros(ComplexF64, 2, 2, length(fs))
    JC.evaluateprovider!(F, fitted.provider, 2pi .* fs)
    @test maximum(abs(F[p, q, i]*cis(-2pi*fs[i]*(taus[p] + taus[q])) - cable(2pi*fs[i])[p, q])
        for i in eachindex(fs), q in 1:2, p in 1:2) < 1e-12
    # the same cable fitted whole cannot be done at this order at all
    @test_throws ArgumentError RationalScattering(block, 4; frequencies = 2pi .* fs)
    # the line of the delay in cascade with the fit is the cable
    cascade = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:la, 1, 2, TransmissionLine(50.0, len; vp = vp)),
        (:b, 2, 3, fitted), (:lb, 3, 4, TransmissionLine(50.0, len; vp = vp)),
        (:p2, 4, 0, Port(2; Z0 = 50.0))])
    check = collect(range(1e9, 9e9; length = 12))
    S = hblinsolve(2pi .* check, cascade; keyedarrays = false).S
    @test maximum(abs(S[2, 1, i] - cable(2pi*check[i])[2, 1]) for i in eachindex(check)) < 1e-9
    # the order search takes them too, and reports the order of what is
    # left rather than of the delay
    searched = RationalScattering(block; tol = 1e-8, minpoles = 1, frequencies = 2pi .* fs, delays = taus)
    @test size(searched.provider.A, 1) == 2
    # no delay and every delay zero are the same fit, on the cable's
    # tail, which needs no delay taken out to be fitted at all
    tail(w) = (g = 1/(1 + im*w/wc); ComplexF64[0 g; g 0])
    tailblock = ScatteringParameters(tail; nports = 2, zref = 50.0)
    plain = RationalScattering(tailblock, 4; frequencies = 2pi .* fs, delays = (0.0, 0.0))
    bare = RationalScattering(tailblock, 4; frequencies = 2pi .* fs, delays = nothing)
    @test plain.provider.A == bare.provider.A && plain.provider.C == bare.provider.C
    # three ports at reference impedances of their own, with
    # reflections, a zero frequency statement and delays which are all
    # different: what is fitted is the core the delays came off, at
    # signed frequencies and at zero, and the delays put back are the
    # data
    let taus3 = [0.125e-9, 0.3e-9, 0.75e-9], M = [0.15 0.2 0.1; 0.2 0.1 0.05; 0.1 0.05 0.2],
            w3 = 2pi*3e9, fs3 = collect(range(0.0, 5e9; length = 100)),
            ws3 = 2pi .* [-4.71e9, -1.11e9, 0.0, 0.31e9, 2.72e9, 4.9e9]
        core(w) = M/(1 + im*w/w3)
        turn(w) = Diagonal(cis.(-w .* taus3))
        data(w) = turn(w)*core(w)*turn(w)
        multi = ScatteringParameters((2pi .* fs3, stack(data.(2pi .* fs3)));
            zref = [40.0, 50.0, 60.0], dcmodel = ScatteringDC(M))
        for fit in (RationalScattering(multi, 4; delays = taus3),
                RationalScattering(multi; tol = 1e-9, minpoles = 1, delays = taus3))
            F = zeros(ComplexF64, 3, 3, length(ws3))
            JC.evaluatescattering!(F, fit, ws3)
            @test all(isapprox(F[:, :, i], core(w); rtol = 1e-10) for (i, w) in enumerate(ws3))
            @test all(isapprox(turn(w)*F[:, :, i]*turn(w), data(w); rtol = 1e-10) for (i, w) in enumerate(ws3))
            @test fit.zref == [40.0, 50.0, 60.0] && isapprox(F[:, :, 3], M; atol = 1e-12)
        end
        # a covariance of uncorrelated ports is what it was, however
        # unequal the delays
        uncorrelated = RationalScattering(ScatteringParameters(data; nports = 3, zref = 50.0,
            noise = NoiseCovariance(Diagonal(fill(2.0, 3)))), 4; frequencies = 2pi .* fs3, delays = taus3)
        V3 = zeros(ComplexF64, 3, 3, 2)
        JC.evaluatecovariance!(V3, uncorrelated, ws3[end - 1:end])
        @test all(V3[:, :, i] ≈ 2I for i in 1:2)
    end
    @test_throws ArgumentError RationalScattering(block, 4; frequencies = 2pi .* fs, delays = (tau,))
    @test_throws ArgumentError RationalScattering(block, 4; frequencies = 2pi .* fs, delays = (tau, -1.0))
    # a covariance the delays rotate is rotated; one they do not is not.
    # The cable's gain reaches one at zero frequency, so a covariance
    # which is to be one the block can emit over the whole band needs
    # its diagonal above what the commutator asks there
    C = [1.0 0.25; 0.25 1.0]
    stated = ScatteringParameters(cable; nports = 2, zref = 50.0,
        noise = NoiseCovariance((2pi .* fs, repeat(complex(C), 1, 1, length(fs)))))
    even = RationalScattering(stated, 4; frequencies = 2pi .* fs, delays = taus, passivity = nothing)
    V = zeros(ComplexF64, 2, 2, 2)
    JC.evaluateprovider!(V, even.noise.provider, 2pi .* fs[2:3])
    @test maximum(abs.(V .- repeat(complex(C), 1, 1, 2))) < 1e-14
    uneven = RationalScattering(stated, 4; frequencies = 2pi .* fs, delays = (tau, 0.0), passivity = nothing)
    JC.evaluateprovider!(V, uneven.noise.provider, 2pi .* fs[2:3])
    @test all(abs(V[1, 2, i] - C[1, 2]*cis(2pi*fs[i + 1]*tau)) < 1e-14 for i in 1:2)
    @test all(abs(V[1, 1, i] - C[1, 1]) < 1e-14 for i in 1:2)
    # the covariance is turned at the frequency it is read at and not at
    # the samples it stores, so between them it is the stated one turned
    # and not the turned samples interpolated
    between = 2pi .* [(fs[2] + fs[3])/2, (fs[40] + fs[41])/2]
    JC.evaluateprovider!(V, uneven.noise.provider, between)
    @test all(abs(V[1, 2, i] - C[1, 2]*cis(between[i]*tau)) < 1e-14 for i in 1:2)
    @test all(abs(V[1, 1, i] - C[1, 1]) < 1e-14 for i in 1:2)
    # and the fitted block behind the line of its delay is the block it
    # was fitted to, in the noise it emits as well as in what it
    # scatters, between the covariance's samples as well as on them
    off = 2pi*(fs[40] + fs[41])/2
    whole = hblinsolve([off], Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:b, 1, 2, stated),
        (:p2, 2, 0, Port(2; Z0 = 50.0))]); keyedarrays = false, returnCnoise = true)
    behind = hblinsolve([off], Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)),
        (:la, 1, 2, TransmissionLine(50.0, len; vp = len/tau)), (:b, 2, 3, uneven),
        (:p2, 3, 0, Port(2; Z0 = 50.0))]); keyedarrays = false, returnCnoise = true)
    @test maximum(abs.(whole.S .- behind.S)) < 1e-9
    @test maximum(abs.(whole.Cnoise .- behind.Cnoise)) < 1e-9
    # the covariance a block states is turned whatever provider holds
    # it: a table of reals, which the turned entries are not, and a
    # callable, which has no samples to turn
    realtable = ScatteringParameters(cable; nports = 2, zref = 50.0,
        noise = NoiseCovariance(JC.TabulatedMatrixProvider(2pi .* fs, repeat(C, 1, 1, length(fs)))))
    called = ScatteringParameters(cable; nports = 2, zref = 50.0,
        noise = NoiseCovariance(w -> complex(C)))
    for held in (realtable, called)
        fit = RationalScattering(held, 4; frequencies = 2pi .* fs, delays = (tau, 0.0), passivity = nothing)
        JC.evaluateprovider!(V, fit.noise.provider, between)
        @test all(abs(V[1, 2, i] - C[1, 2]*cis(between[i]*tau)) < 1e-14 for i in 1:2)
    end
    # and beyond a table it keeps turning with the frequency, being a
    # phase on the covariance the extrapolation states there
    narrow = ScatteringParameters(cable; nports = 2, zref = 50.0,
        noise = NoiseCovariance((2pi .* fs[1:20], repeat(complex(C), 1, 1, 20)); extrapolation = :constant))
    beyond = RationalScattering(narrow, 4; frequencies = 2pi .* fs, delays = (tau, 0.0), passivity = nothing)
    out = 2pi .* [fs[60], fs[80]]
    JC.evaluateprovider!(V, beyond.noise.provider, out)
    @test all(abs(V[1, 2, i] - C[1, 2]*cis(out[i]*tau)) < 1e-14 for i in 1:2)
    # a turned covariance is still data the commutation relations can be
    # checked against: this block has a power gain of four and states the
    # turned covariance of the fitted block, which is far less than it
    # must emit
    @test_throws ArgumentError ScatteringParameters(ComplexF64[0 2; 2 0]; zref = 50.0,
        noise = NoiseCovariance(uneven.noise.provider))
    # and it is checked on the covariance the rotation returns, not on
    # the one it turns: a phase which is not real turns a Hermitian
    # covariance into one which is not, and one which is not into one
    # which is
    phased(A) = NoiseCovariance(JC.RotatedMatrixProvider(JC.ConstantMatrixProvider(A), [0.0, 0.0]; phase = cis(pi/3)))
    @test_throws ArgumentError ScatteringParameters([0.5 0.0; 0.0 0.5]; zref = 50.0, noise = phased(complex(C)))
    @test ScatteringParameters([0.5 0.0; 0.0 0.5]; zref = 50.0, noise = phased(cis(-pi/3)*C)) isa ScatteringParameters
    # the same provider carries the rotation a pumped fit makes, the
    # emitting side read a harmonic above the frequency asked for and a
    # constant factor for the pump phase, and it carries the conjugate
    # ladder relation of a pumped covariance exactly
    let taus4 = [0.4e-9, 0.1e-9], off = 2pi*3e9, ph = cis(0.7), nu = 2pi*1.3e9
        inner = JC.TabulatedMatrixProvider(2pi .* fs, repeat(complex(C), 1, 1, length(fs));
            interpolation = :linear)
        turned = JC.RotatedMatrixProvider(inner, taus4; offset = off, phase = ph)
        ws = 2pi .* [(fs[2] + fs[3])/2, fs[40]]
        A, B = zeros(ComplexF64, 2, 2, 2), zeros(ComplexF64, 2, 2, 2)
        JC.evaluateprovider!(A, inner, ws)
        JC.evaluateprovider!(B, turned, ws)
        @test all(abs(B[p, q, i] - ph*cis((ws[i] + off)*taus4[p] - ws[i]*taus4[q])*A[p, q, i]) < 1e-14
            for i in 1:2, q in 1:2, p in 1:2)
        symmetric = JC.RotatedMatrixProvider(JC.ConstantMatrixProvider(complex(C)), taus4;
            offset = off, phase = ph)
        at, image = zeros(ComplexF64, 2, 2, 1), zeros(ComplexF64, 2, 2, 1)
        JC.evaluateprovider!(at, symmetric, [nu])
        JC.evaluateprovider!(image, symmetric, [-nu - off])
        @test maximum(abs.(image[:, :, 1] .- transpose(at[:, :, 1]))) < 1e-14
    end
end

@testset "a pumped block fitted with its zero frequency statement" begin
    wp = 2pi*1e9
    zero1 = zeros(ComplexF64, 1, 1)
    # the fit of the unconverted response meets the stated open at zero
    # frequency, which the data's extrapolation would not; the statement
    # contradicts the data near zero, so the fit is far from the data
    # there and the noise it adds to obey the commutation relations is
    # large
    b = LinearizedScattering([w -> fill(im*w/(wp + im*w), 1, 1), zero1], wp; harmonics = [0, 1], nports = 1,
        noise = NoiseCovariance([fill(2.0, 1, 1), zero1]), dcmodel = OpenDC())
    f = RationalScattering(b, 2; frequencies = 2pi .* collect(range(0.1e9, 0.9e9; length = 10)), tol = 10.0, noisetol = 10.0)
    H = zeros(ComplexF64, 1, 1, 1)
    JosephsonCircuits.evaluateprovider!(H, f.providers[1], [0.0])
    @test H[1] ≈ 1
    @test f.dcmodel isa OpenDC
    @test f.atol == b.atol && f.noise.completed && f.noise.atol == b.noise.atol
    # a band which leaves the unconverted response no sample is refused
    lossless = LinearizedScattering([w -> fill((wp - im*w)/(wp + im*w), 1, 1), zero1], wp; harmonics = [0, 1], nports = 1, dcmodel = OpenDC())
    @test_throws ArgumentError RationalScattering(lossless, 2; frequencies = 2pi .* [0.1e9, 0.2e9, 0.3e9], band = (2pi*0.5e9, 2pi*0.9e9))
end

# The fit of a pumped block checks the declaration of its data, and the
# noise its fit adds, over the modes its frequencies reach, ladder by
# ladder of the pump, the modes of each formed from its least frequency.
# A JPA solved with sixteen modulation harmonics holds sidebands of its
# band many pumps up, where the roundoff of a frequency exceeds the
# tolerance at the edge of a table; formed from the least frequency of
# each ladder, its modes fall within its tables, and the block is
# fitted.
@testset "a pumped block solved over many sidebands is fitted" begin
    jpa = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:cc, 1, 2, Capacitor(100.0e-15)),
        (:jj, 2, 0, JosephsonJunction(1000.0e-12)), (:cj, 2, 0, Capacitor(1000.0e-15))])
    fp, ip, Nmod = 4.75e9, 0.00565e-6, 16
    sol = hbsolve(2pi*collect(range(4.0e9, 5.5e9; length = 11)), (2pi*fp,), [(mode = (1,), port = 1, current = ip)],
        (Nmod,), (2Nmod,), jpa; atol = 1e-14)
    @test RationalScattering(LinearizedScattering(sol.linearized, 2pi*fp), 10) isa LinearizedScattering
end

# A residue of a fit takes as many states as it has singular values above
# `ranktol` of its largest, and the part it loses outside them takes its
# share of the value at zero frequency with it, which the constant takes
# back: a value stated there is met by the raw fit, by a fit which states
# its noise and by the fit of a pumped block's unconverted response, as by
# the enforced one. A pole pair whose residue's second direction is 1/150
# of its first keeps one state per pole at a `ranktol` of 1e-2.
@testset "a stated zero frequency value through the residues' rank" begin
    JC = JosephsonCircuits
    Q = [cos(0.3) -sin(0.3); sin(0.3) cos(0.3)]
    a = complex(-0.4e10, 2.0e10)
    R = Q*Diagonal([0.3e10, 0.002e10])*Q' .* (1 + 0.7im)
    S(w) = [0.1 0.05; 0.05 0.1] .+ R ./ (im*w - a) .+ conj.(R) ./ (im*w - conj(a))
    S0 = real.(S(0.0))
    ws = 2pi .* collect(range(0.5e9, 8e9; length = 120))
    fitting = VectorFitting(ranktol = 1e-2)
    at0(p) = JC.evaluateprovider!(zeros(ComplexF64, 2, 2, 1), p, [0.0])[:, :, 1]
    table(; kw...) = ScatteringParameters((ws, stack(S.(ws))); nports = 2, zref = 50.0,
        dcmodel = ScatteringDC(S0), kw...)
    for (block, passivity) in ((table(), nothing), (table(), PassivityEnforcement()),
            (table(noise = NoiseCovariance(Matrix(1.0I, 2, 2))), nothing))
        fit = RationalScattering(block, 2; fitting, passivity)
        @test size(fit.provider.A, 1) == 2
        @test maximum(abs, at0(fit.provider) .- S0) < 1e-12
    end
    pumped = LinearizedScattering([S], 2pi*30e9; harmonics = [0], nports = 2,
        noise = NoiseCovariance([Matrix{ComplexF64}(I, 2, 2)]), dcmodel = ScatteringDC(S0))
    fit = RationalScattering(pumped, 2; frequencies = ws, fitting)
    @test size(fit.providers[1].A, 1) == 2
    @test maximum(abs, at0(fit.providers[1]) .- S0) < 1e-12
end

# A stated value at zero frequency which settles the only coefficient, the
# constant of a fit without poles or the residue of one real pole under a
# held constant, is met, weighted or not; so the pruning keeps a pole the
# samples barely show where it alone reaches the stated value, which the
# constant meeting the samples misses.
@testset "a stated zero frequency value which settles the only coefficient" begin
    JC = JosephsonCircuits
    ws = collect(range(1.0, 2.0; length = 21))
    S = fill(0.5 + 0.0im, 1, 1, length(ws))
    for weights in (JC.noweights, ones(1, 1, length(ws)))
        X, _ = JC.fitcoefficients(S, ws, ComplexF64[]; dc = fill(0.2, 1, 1), weights)
        @test abs(only(X) - 0.2) < 1e-12
        R, D = JC.fitresidues(S, ws, ComplexF64[-1]; constant = zeros(1, 1), dc = fill(0.2, 1, 1), weights)
        @test abs(only(D .+ R[:, :, 1]) - 0.2) < 1e-12
    end
    e = 1e-13
    block = ScatteringParameters((ws, reshape([0.5 - 0.3e/(im*w + e) for w in ws], 1, 1, :));
        dcmodel = ScatteringDC(fill(0.2, 1, 1)))
    fit = RationalScattering(block, 1; fitting = VectorFitting(start = [complex(-e)], iterations = 1), tol = 1e-8)
    @test abs(only(JC.evaluateprovider!(zeros(ComplexF64, 1, 1, 1), fit.provider, [0.0])) - 0.2) < 1e-12
end

# An amplifier stating the least noise a quantum limited one emits,
# |I - S S'|/2, whose gain carries a ripple of 2e-3 no low order follows:
# its fit misses the samples by its error, and the stated covariance is
# completed to the commutation relations of the fitted scattering matrix,
# which a solve then meets: the noise the block adds at port 2 is the
# stated V22, raised to the quantum limit of the fitted gain,
# (|S21|^2 - 1)/2, where that exceeds the data's. A covariance the data
# itself does not meet is refused before the fit, and the noise the
# completion adds is held to `noisetol`.
@testset "a fitted block which states its noise" begin
    JC = JosephsonCircuits
    w0, g0 = 2pi*8e9, 10.0
    gain(w) = g0*w0/(w0 + im*w)*(1 + 2e-3*sin(w/(0.7e9*2pi)))
    amp(w) = [0.0 0.0; gain(w) 0.0]
    limit(w) = Diagonal([0.5, (abs2(gain(w)) - 1)/2])
    ws = 2pi .* collect(range(0.1e9, 12e9; length = 300))
    data = ScatteringParameters((ws, stack(amp.(ws))); nports = 2, zref = 50.0,
        noise = NoiseCovariance((ws, stack(complex.(limit.(ws))))))
    fit = RationalScattering(data, 2)
    @test fit.noise.completed
    @test 1e-3 < JC.relativefiterror(fit, stack(amp.(ws)), ws) < 1e-2
    check = 2pi .* [1e9, 3e9, 5e9]
    out = hblinsolve(check, Circuit([(:p1, 1, 0, Port(1)), (:a, 1, 2, fit), (:p2, 2, 0, Port(2))]);
        keyedarrays = false, returnCnoise = true)
    S21 = out.S[2, 1, :]
    @test maximum(abs.(S21 .- gain.(check))) < 0.03
    Vtab = zeros(ComplexF64, 2, 2, length(check))
    JC.evaluatecovariance!(Vtab, data, check)
    @test real.(out.Cnoise[2, 2, :]) ≈ max.(real.(Vtab[2, 2, :]), (abs2.(S21) .- 1)./2) rtol = 1e-9
    @test_throws ArgumentError RationalScattering(data, 2; noisetol = 1e-5)
    @test_throws ArgumentError RationalScattering(data; tol = 1e-2, minpoles = 1, maxpoles = 2, noisetol = 1e-5)
    @test RationalScattering(data; tol = 1e-2, minpoles = 1).noise.completed
    below = ScatteringParameters(amp; nports = 2, zref = 50.0, noise = NoiseCovariance(w -> 0.9 .* limit(w)))
    @test_throws ArgumentError RationalScattering(below, 2; frequencies = ws)
    # A covariance stated completed carries the noise the block's own
    # scattering requires, which the fit carries as well and is not
    # charged for (review of 2026-10-07, finding 1): a lossy one-port and
    # an amplifying one, stated as zero completed, fit as with that noise
    # written out, and in a solve both fits add at the port the least noise
    # a lossy or an amplifying one-port can, |1 - |S|^2|/2
    wc = 2pi*1e9
    band = wc .* collect(range(0.1, 2.0; length = 40))
    for s in (0.5, 2.0)
        response(w) = fill(s/(1 + im*w/wc), 1, 1)
        completed = ScatteringParameters(response; nports = 1, zref = 50.0,
            noise = NoiseCovariance(zeros(1, 1); completed = true))
        written = ScatteringParameters(response; nports = 1, zref = 50.0,
            noise = NoiseCovariance(w -> fill(abs(1 - abs2(s/(1 + im*w/wc)))/2, 1, 1)))
        @test RationalScattering(completed; tol = 1e-2, minpoles = 1, maxpoles = 4, frequencies = band).noise.completed
        for block in (completed, written)
            fitted = RationalScattering(block, 2; frequencies = band)
            out = hblinsolve(wc .* [0.5, 1.5], Circuit([(:p1, 1, 0, Port(1)), (:a, 1, fitted)]);
                keyedarrays = false, returnCnoise = true)
            @test real.(out.Cnoise[1, 1, :]) ≈ abs.(1 .- abs2.(out.S[1, 1, :]))./2 rtol = 1e-9
        end
    end
    # and so for a pumped block, over the ladder its fit completes on: the
    # lossy one-port unconverted fits stated either way, the two fits
    # stating the same noise, and converting, stated completed, it fits
    h0(w) = fill(0.5/(1 + im*w/wc), 1, 1)
    pumped(providers, harmonics, noise) = RationalScattering(LinearizedScattering(providers, 2pi*20e9;
        harmonics, nports = 1, noise), 2; frequencies = band)
    a = pumped([h0], [0], NoiseCovariance([zeros(ComplexF64, 1, 1)]; completed = true))
    b = pumped([h0], [0], NoiseCovariance([w -> fill(complex((1 - abs2(only(h0(w))))/2), 1, 1)]))
    rows, cols, K = JC.pumpedfamily(a, (wc,))
    @test JC.pumpednoisematrices(a, rows, cols, K)[3] ≈ JC.pumpednoisematrices(b, rows, cols, K)[3] rtol = 1e-9
    @test pumped([h0, w -> fill(0.1/(1 + im*w/wc), 1, 1)], [0, 1],
        NoiseCovariance([zeros(ComplexF64, 1, 1), zeros(ComplexF64, 1, 1)]; completed = true)) isa LinearizedScattering
end

# A fit refuses what it cannot fit by what is wrong with it and where: a
# block which is not finite at a sample, in the ordinary fit and the order
# search and in the pumped fit's harmonics, and a pumped block whose data
# does not meet what it declares, which is checked before any harmonic is
# fitted, here a lossy delay declared lossless which one pole cannot fit.
# The refusals name the frequency, which those of the fits they forestall
# do not.
@testset "a fit refuses its data by where it fails" begin
    ws = 2pi .* collect(range(0.5e9, 5e9; length = 10))
    h(w) = fill(w == ws[3] ? NaN + 0im : 0.5/(1 + im*w/(2pi*3e9)), 1, 1)
    refusal(f) = try
        f(); ""
    catch e
        e isa ArgumentError || rethrow()
        sprint(showerror, e)
    end
    bad = ScatteringParameters(h; nports = 1, zref = 50.0)
    stated = NoiseCovariance([fill(1.0 + 0im, 1, 1), zeros(ComplexF64, 1, 1)])
    for f in (() -> RationalScattering(bad, 2; frequencies = ws),
            () -> RationalScattering(bad; tol = 1e-3, frequencies = ws),
            () -> RationalScattering(LinearizedScattering([h, w -> zeros(ComplexF64, 1, 1)], 2pi*20e9;
                harmonics = [0, 1], nports = 1, noise = stated), 2; frequencies = ws),
            () -> RationalScattering(LinearizedScattering([w -> fill(0.1 + 0im, 1, 1), h], 2pi*20e9;
                harmonics = [0, 1], nports = 1, noise = stated), 2; frequencies = ws))
        @test occursin(string(ws[3]), refusal(f))
    end
    lossy = LinearizedScattering([w -> fill(0.5*cis(-w*1e-9), 1, 1)], 2pi*20e9; harmonics = [0], nports = 1)
    @test occursin(string(ws[1]), refusal(() -> RationalScattering(lossy, 1; frequencies = ws)))
end

@testset "a nearly defective realization is evaluated accurately" begin
    # an all pass whose state matrix is within 3e-10 of a Jordan block:
    # its eigenvectors are nearly parallel, and the resolvent over many
    # frequencies must be as accurate as a solve at each one, and as
    # lossless
    delta = 3e-10
    A = [-1.0 1.0; 0.0 -1.0 - delta]
    B = reshape([0.0, 1.0], 2, 1)
    C = reshape([4 + 2delta, -4 - 2delta], 1, 2)
    D = ones(1, 1)
    b = RationalScattering(A, B, C, D; noise = Lossless())
    ws = collect(range(0.01, 2.0; length = 1000))
    a = zeros(ComplexF64, 1, 1, length(ws))
    JosephsonCircuits.evaluateprovider!(a, b.provider, ws)
    direct = [only(D + C*((im*w*I - A) \ B)) for w in ws]
    @test maximum(abs, vec(a) .- direct) < 1e-12
    @test maximum(abs, abs2.(a) .- 1) < 1e-12
    c = Circuit([(:p, 1, 0, Port(1)), (:b, 1, b)])
    for count in (8, 9)
        o = hblinsolve(collect(range(0.01, 2.0; length = count)), c; keyedarrays = false)
        @test maximum(abs, abs.(o.S) .- 1) < 1e-10
    end
end

# The passivity correction is the least change of the residues and the
# constant over the samples, solved in the residues' own coordinates at
# the frequency scaled to the band, where the problem carries no unit.
@testset "the passivity correction is the least change of the residues" begin
    JC = JosephsonCircuits
    response(A, B, C, D, w) = D .+ C*((im*w*I - A) \ B)
    # the largest singular value sampled densely past the poles and
    # across every resonance, independent of the norm search
    function densemax(A, B, C, D)
        poles = eigvals(A)
        mags = filter(>(0), abs.(poles))
        ws = vcat(0.0, exp.(range(log(minimum(mags)/100), log(maximum(mags)*100); length = 4000)))
        for l in unique(round.(poles; sigdigits = 10))
            imag(l) > 0 && append!(ws, imag(l) .+ abs(real(l)) .* range(-8, 8; length = 401))
        end
        return maximum(w -> opnorm(response(A, B, C, D, w)), ws)
    end
    densemax(p::JC.RationalScatteringProvider) = densemax(p.A, p.B, p.C, p.D)
    # the poles and residues of a small realization, the form the fit
    # hands the enforcement
    function poleresidues(A, B, C)
        λ, V = eigen(A)
        W = inv(V)
        return ComplexF64.(λ), ComplexF64.(cat([C*V[:, k]*transpose(W[k, :])*B for k in eachindex(λ)]...; dims = 3))
    end
    realized(poles, R, D) = (JC.realization(poles, R, size(D, 1))..., D)
    # `nl` coupled lines of ten lumped sections, ports at both ends,
    # nearly lossless (a dissipation of about 1e-6)
    function lumpedlines(nl)
        node(l, k) = (l - 1)*11 + k + 1
        L, Cs = 50.0*10e-12, 10e-12/50.0
        comps = Any[]
        for l in 1:nl
            push!(comps, (Symbol(:pa, l), node(l, 0), 0, Port(l; Z0 = 50.0)))
            push!(comps, (Symbol(:pb, l), node(l, 10), 0, Port(nl + l; Z0 = 50.0)))
            for k in 1:10
                push!(comps, (Symbol(:L, l, :_, k), node(l, k - 1), node(l, k), Inductor(L*(1 + 0.03l))))
                push!(comps, (Symbol(:G, l, :_, k), node(l, k), 0, Resistor(5e4)))
                push!(comps, (Symbol(:C, l, :_, k), node(l, k), 0, Capacitor(Cs*(k == 10 ? 0.5 : 1.0))))
                l < nl && push!(comps, (Symbol(:Cm, l, :_, k), node(l, k), node(l + 1, k), Capacitor(0.15Cs)))
            end
            push!(comps, (Symbol(:C0, l), node(l, 0), 0, Capacitor(0.5Cs)))
        end
        return Circuit(comps)
    end
    fs = collect(range(20e9/400, 20e9; length = 400))
    S2 = hblinsolve(2pi .* fs, lumpedlines(1); keyedarrays = false).S
    # At twenty poles the raw fit of the line is good to 2e-10 and its
    # only violations lie above the band, where its norm reaches 1.075:
    # the least change removing them costs next to nothing in band, and
    # carries no unit of frequency.
    f2 = RationalScattering(ScatteringParameters((2pi .* fs, S2); nports = 2, zref = 50.0), 20; tol = 1.0)
    e2 = JC.relativefiterror(f2, S2, 2pi .* fs)
    @test e2 < 1e-5
    # the largest spectral norm over the samples, found with a singular
    # value decomposition only where the Frobenius norm allows it to be
    # the largest, is the largest of them all
    M = randn(ComplexF64, 5, 5, 200) .* reshape(exp.(range(-3, 3; length = 200)), 1, 1, :)
    @test JC.largestopnorm(M) == maximum(k -> opnorm(M[:, :, k]), axes(M, 3))
    @test densemax(f2.provider) <= 1 + 2e-8
    # the fit's state matrix is block diagonal in real poles and
    # rotations, its residues of rank two, and a frequency is evaluated
    # pole by pole
    @test isempty(f2.provider.states.poles)
    # the same responses a billion times lower in frequency: the same fit
    f2s = RationalScattering(ScatteringParameters((2pi .* fs .* 1e-9, S2); nports = 2, zref = 50.0), 20; tol = 1.0)
    @test JC.relativefiterror(f2s, S2, 2pi .* fs .* 1e-9) ≈ e2 rtol = 1e-2
    # The data is reciprocal and so is the raw fit; the correction keeps
    # it so to within the fit's own error, beyond the band as well, where
    # nothing in the samples holds it.
    wo = 2pi .* exp.(range(log(fs[1]/10), log(fs[end]*10); length = 200))
    F = zeros(ComplexF64, 2, 2, length(wo))
    JC.evaluateprovider!(F, f2.provider, wo)
    @test maximum(k -> opnorm(F[:, :, k] - transpose(F[:, :, k])), axes(F, 3)) <= e2
    # Two coupled lines at twenty-four poles, whose first correction
    # would leave the constant term with a singular value of 1.013 were
    # it left alone, which no sweep over the band sees: the constant is
    # constrained in every round.
    S4 = hblinsolve(2pi .* fs, lumpedlines(2); keyedarrays = false).S
    f4 = RationalScattering(ScatteringParameters((2pi .* fs, S4); nports = 4, zref = 50.0), 24; tol = 1e-2)
    @test densemax(f4.provider) <= 1 + 2e-8
    # At forty poles the same lines' fit is active only above the band, where
    # a change costs little over the samples and moves the response freely:
    # each round's correction pushes points above one that an earlier round
    # brought under, unless those stay constrained. They do, so the rounds
    # end on a passive fit, near the raw one, rather than on a contraction.
    f40 = RationalScattering(ScatteringParameters((2pi .* fs, S4); nports = 4, zref = 50.0), 40; tol = 1.0)
    @test JC.relativefiterror(f40, S4, 2pi .* fs) <= 2e-6
    # Eight lines at forty poles meet new violations above the band round
    # after round, fewer and smaller, and are passive after about forty:
    # the rounds the default allows, where fewer would end on a contraction.
    S8 = hblinsolve(2pi .* fs, lumpedlines(4); keyedarrays = false).S
    f8 = RationalScattering(ScatteringParameters((2pi .* fs, S8); nports = 8, zref = 50.0), 40; tol = 1.0)
    @test JC.relativefiterror(f8, S8, 2pi .* fs) <= 3e-4
    # The certificate is the sweep that settles every interval of the axis
    # by a bound across it (see passivitysweep). Its bands are the crossing
    # test's, on the raw fit of the two lines at twenty-four poles and on
    # the same fit turned by an orthogonal matrix on the left, which keeps
    # its singular values and makes it nonreciprocal: every band between
    # two crossings where the response stands above the level meets one of
    # the sweep's, and every band of the sweep meets one of those.
    level = 1 + 5e-9
    ws4 = 2pi .* fs
    p4, R4, D4 = JC.vectorfit(S4, ws4, 24, VectorFitting(); constanttol = 1e-8, constantmargin = 1e-6)
    Q4 = [cos(0.4) -sin(0.4) 0 0; sin(0.4) cos(0.4) 0 0; 0 0 cos(1.1) sin(1.1); 0 0 -sin(1.1) cos(1.1)]
    for (Rk, Dk) in ((R4, D4), (stack(Q4*R4[:, :, k] for k in axes(R4, 3)), Q4*D4))
        form = JC.ResidueForm(p4, Rk, Dk, sqrt(ws4[1]*ws4[end]), JC.residuespaces(p4, Rk, 1e-12))
        sweep = JC.passivitysweep(form, level, Inf)
        @test sweep.verdict === :active
        A4, B4, C4 = JC.realization(form)
        An4, Bn4, Cn4, w4 = JC.balancedrealization(A4, B4, C4 ./ level)
        crossed = JC.unitcrossings(An4, Bn4, Cn4, form.D ./ level) .* w4
        F4, phi4 = zeros(ComplexF64, 4, 4), zeros(ComplexF64, length(form.poles))
        edges = vcat(0.0, crossed)
        above = [(edges[k], edges[k + 1]) for k in 1:length(edges) - 1
            if opnorm(JC.responseat!(F4, form, (edges[k] + edges[k + 1])/2, phi4)) > level]
        meets(a, b) = a[1] <= b[2] && b[1] <= a[2]
        @test !isempty(above)
        @test all(b -> any(s -> meets(s, b), sweep.bands), above)
        @test all(s -> any(b -> meets(s, b), above), sweep.bands)
    end
    # A band narrower than any grid: a resonance of half width 1e-9 whose
    # peak stands 1e-7 above one, found at its level and settled above it.
    let d = 1e-9, R = zeros(ComplexF64, 2, 2, 2)
        R[1, 1, 1] = R[1, 1, 2] = (0.5 + 1e-7)*d
        narrow = JC.ResidueForm(ComplexF64[complex(-d, 1), complex(-d, -1)], R, [0.5 0.0; 0.0 0.2], 1.0,
            [Matrix{ComplexF64}(I, 2, 2) for _ in 1:2])
        sweep = JC.passivitysweep(narrow, level, Inf)
        @test sweep.verdict === :active && any(b -> b[1] <= 1 <= b[2], sweep.bands)
        @test JC.passivitysweep(narrow, 1 + 2e-7, Inf).verdict === :passive
    end
    # A band reaching past the poles into the variable `1/x`: a constant at
    # `1 - 1e-6` and `0.5/(s + 1)` on the first port, whose magnitude
    # squared, `(1 - 1e-6 + 0.5/(1 + x^2))^2 + (0.5 x/(1 + x^2))^2`, stands
    # above the level until about `x = 790`.
    let tail = JC.ResidueForm(ComplexF64[-1], reshape(ComplexF64[0.5 0; 0 0], 2, 2, 1), [1 - 1e-6 0.0; 0.0 0.3],
            1.0, [Matrix{ComplexF64}(I, 2, 2)])
        sweep = JC.passivitysweep(tail, level, Inf)
        @test sweep.verdict === :active
        @test all(x -> any(b -> b[1] <= x <= b[2], sweep.bands), (10.0, 100.0, 700.0))
    end
    # The roundoff allowance counts a pole pair's terms before they cancel.
    # The pair -1 +- 1.1e-6 i with coefficients [0, 1/2.2e-6] has as its
    # second basis function i (u - v), whose terms stand six orders of
    # magnitude above it at the frequency one and past the poles: at an
    # interval there, in x and in the tail's 1/x, the dissipation's ends
    # are within the allowance of a 256-bit evaluation.
    let pp = ComplexF64[complex(-1, 1.1e-6), complex(-1, -1.1e-6)], Rp = zeros(ComplexF64, 1, 1, 2)
        Rp[1, 1, 1], Rp[1, 1, 2] = complex(0, 1/2.2e-6), complex(0, -1/2.2e-6)
        pair = JC.ResidueForm(pp, Rp, fill(0.8660254, 1, 1), 1.0, [ones(ComplexF64, 1, 1) for _ in 1:2])
        work = JC.SweepWork(pair)
        setprecision(BigFloat, 256) do
            for (c, tail) in ((1.0, false), (0.02, true))
                h = 1e-9*c
                M, roundoff, _ = JC.expand!(work, pair, c, h, tail)
                q = Complex{BigFloat}(pair.poles[1])
                term(a) = tail ? (big(c)/(im - a*big(c)), im/(im - a*big(c))^2) : (1/(im*big(c) - a), -im/(im*big(c) - a)^2)
                (u, du), (v, dv) = term(q), term(conj(q))
                X1, X2 = big(pair.X[1, 1, 1]), big(pair.X[1, 1, 2])
                S0 = big(pair.D[1, 1]) + X1*(u + v) + X2*im*(u - v)
                S1 = X1*(du + dv) + X2*im*(du - dv)
                @test abs(real(work.P0[1, 1]) - (1 - abs2(S0))) + h*abs(2real(conj(S0)*S1) - 2real(work.G[1, 1])) <= roundoff
            end
            # and the pair's functions are formed without the difference of
            # its terms, so the response there is accurate to its own
            # roundoff, where the difference erred by 3.6e-11
            u, v = 1/(im*big(1.0) - Complex{BigFloat}(pp[1])), 1/(im*big(1.0) - Complex{BigFloat}(pp[2]))
            exact = big(pair.D[1, 1]) + big(pair.X[1, 1, 1])*(u + v) + big(pair.X[1, 1, 2])*im*(u - v)
            @test abs(JC.responseat!(zeros(ComplexF64, 1, 1), pair, 1.0, zeros(ComplexF64, 2))[1, 1] - exact) <= 1e-14
        end
    end
    # An all pass, `(1 - s)/(1 + s)`, dissipates nothing at any frequency:
    # it is settled at a level above one, and at one it is not, where it has
    # nothing to spare and the sweep halves until its time runs out.
    allpass = JC.ResidueForm(ComplexF64[-1], fill(2.0 + 0im, 1, 1, 1), fill(-1.0, 1, 1), 1.0, [ones(ComplexF64, 1, 1)])
    @test JC.passivitysweep(allpass, 1 + 1e-6, Inf).verdict === :passive
    @test JC.passivitysweep(allpass, 1.0, time_ns() + 1e7).verdict !== :passive
    # Three ports lossless in some directions at every frequency and lossy
    # in the others, `Q diag(...) Q'`, whose first-order bound dips along
    # the lossless directions until the intervals are narrow: the bound to
    # second order settles them in a tenth of the intervals, with one
    # lossless direction, `(1 - s)/(1 + s)` beside `0.5/(s + 1)` and
    # `0.25/(s + 1)`, where its test takes the basis of the lossless span,
    # and with two, `(2 - s)/(2 + s)` in place of the second, where it takes
    # that of the lossy one and tests along both. A resonance of half width
    # 1e-4 at 0.7 along the first direction, lifting it 2e-8 above one, is
    # found, its band about the resonance, and one lifting it 2e-9, under
    # the level, is settled; the bracket of the norm holds the peak of a
    # fine grid.
    let Q = Matrix(qr(reshape(sin.((1:9).^2), 3, 3)).Q), u = Q[:, 1], level = 1 + 5e-9
        turned(d) = Q*Diagonal(d)*Q'
        for (poles, X, D) in ((ComplexF64[-1], [turned([2.0, 0.5, 0.25])], turned([-1.0, 0, 0])),
                (ComplexF64[-1, -2], [turned([2.0, 0, 0.5]), turned([0, 4.0, 0])], turned([-1.0, -1, 0])))
            form(p, extra...) = JC.ResidueForm(p, ComplexF64.(cat(X..., extra...; dims = 3)), D, 1.0,
                [Matrix{ComplexF64}(I, 3, 3) for _ in p])
            @test JC.passivitysweep(form(poles), level, Inf).evaluated < 2000
            for (lift, verdict) in ((2e-8, :active), (2e-9, :passive))
                R = lift*1e-4*(1 - 0.7im)/(1 + 0.7im)
                bumped = form(vcat(poles, complex(-1e-4, 0.7), complex(-1e-4, -0.7)), R .* (u*u'), conj(R) .* (u*u'))
                sweep = JC.passivitysweep(bumped, level, Inf)
                @test sweep.verdict === verdict
                @test all(b -> b[1] <= 0.7 <= b[2], sweep.bands)
                F, phi = zeros(ComplexF64, 3, 3), zeros(ComplexF64, length(bumped.poles))
                peak = maximum(x -> opnorm(JC.responseat!(F, bumped, x, phi)), range(0.698, 0.702; length = 4001))
                lower, _, ceiling = JC.residuenorm(bumped; rtol = 1e-8)
                @test lower <= peak*(1 + 1e-12) && peak <= ceiling <= lower*(1 + 2e-8)
            end
        end
    end
    # The test to second order works in the arrays of its work: on sixteen
    # ports, eight two ports lossless in one direction,
    # `Q diag((1 - s)/(1 + s), 0.5/(s + 1)) Q'`, turned by an orthogonal
    # matrix, an attempt which settles an interval the first order leaves
    # allocates less than one of the form's matrices, the LAPACK routines'
    # workspace.
    let Qr = [cos(0.3) -sin(0.3); sin(0.3) cos(0.3)], Q16 = Matrix(qr(reshape(sin.((1:256).^2), 16, 16)).Q)
        part = JC.ResidueForm(ComplexF64[-1], reshape(ComplexF64.(Q16*kron(I(8), Qr*[2.0 0; 0 0.5]*Qr')*Q16'), 16, 16, 1),
            Q16*kron(I(8), Qr*[-1.0 0; 0 0]*Qr')*Q16', 1.0, [Matrix{ComplexF64}(I, 16, 16)])
        work, level, c, h = JC.SweepWork(part), 1 + 5e-9, 0.15, 0.003
        function attempt()
            M, roundoff, magnitude, e = JC.expand!(work, part, c, h, false)
            shift, allowance = JC.testshift(work, h, M, roundoff, magnitude, level)
            return JC.settled(work, h, shift),
                JC.secondorder!(work, part, c, h, false, (level - 1)*(level + 1), allowance, e, h^2*M/64)
        end
        @test attempt() == (false, true)
        attemptbytes() = @allocated attempt()
        @test attemptbytes() < sizeof(ComplexF64)*16^2
    end
    # A lossless form is settled by the bound of its dissipation over the
    # whole axis, from its residues and constant, before any interval, and
    # its norm with it: an all pass of two poles turned by an orthogonal
    # matrix, and two lossless sections whose projectors do not commute,
    # `(I - 2P/(s + 1))(I - 6Q/(s + 3))`, at six ports. A form 1e-7 above one
    # at zero frequency with a unitary constant is not settled by it, and
    # the sweep finds it active.
    let n = 6, Q = Matrix(qr(reshape(sin.((1:36).^2), 6, 6)).Q)
        U, V = Matrix(qr(reshape(cos.((1:18).^2), 6, 3)).Q), Matrix(qr(reshape(sin.((2:13).^2), 6, 2)).Q)
        P, Pq = U*U', V*V'
        spaces = [Matrix{ComplexF64}(I, n, n) for _ in 1:2]
        cascade = JC.ResidueForm(ComplexF64[-1, -3], ComplexF64.(cat(4Q, -12Q; dims = 3)), Q, 1.0, spaces)
        sections = JC.ResidueForm(ComplexF64[-1, -3], ComplexF64.(cat(-2P + 6P*Pq, -6Pq - 6P*Pq; dims = 3)),
            Matrix(1.0I, n, n), 1.0, spaces)
        for form in (cascade, sections)
            sweep = JC.passivitysweep(form, 1 + 5e-9, Inf)
            @test sweep.verdict === :passive && sweep.evaluated == 0
            lower, _, ceiling = JC.residuenorm(form; rtol = 1e-8)
            @test 1 - 1e-12 <= lower <= ceiling <= 1 + 2e-8
        end
        active = JC.ResidueForm(ComplexF64[-1], fill(-(2 + 1e-7) + 0im, 1, 1, 1), fill(1.0, 1, 1), 1.0, [ones(ComplexF64, 1, 1)])
        @test JC.passivitysweep(active, 1 + 5e-9, Inf).verdict === :active
    end
    # A lossless fit, whose largest singular value stands at one throughout,
    # is settled before any interval: three coupled resonators between
    # matched ports, fitted at their six poles, are certified with no time
    # for intervals at all, and are passive on a dense grid.
    let comps = Any[(:p1, 1, 0, Port(1; Z0 = 50.0)), (:p2, 3, 0, Port(2; Z0 = 50.0))]
        for r in 1:3
            push!(comps, (Symbol(:L, r), r, 0, Inductor(1e-9)),
                (Symbol(:C, r), r, 0, Capacitor(1/((2pi*5e9*(1 + 0.02*(r - 2)))^2*1e-9))))
            r < 3 && push!(comps, (Symbol(:Cc, r), r, r + 1, Capacitor(0.08e-12)))
        end
        fc = collect(range(4e9, 6e9; length = 100))
        Sc = hblinsolve(2pi .* fc, Circuit(comps); keyedarrays = false).S
        chain = RationalScattering(ScatteringParameters((2pi .* fc, Sc); nports = 2, zref = 50.0), 6;
            tol = 1e-6, passivity = PassivityEnforcement(maxtime = 0.0))
        @test densemax(chain.provider) <= 1 + 5e-9
    end
    # The norm the enforcement decides on is found by the same bounds, the
    # response at each interval's centre raising a lower bound until every
    # interval stands under it times `1 + rtol`; its bracket agrees with
    # the crossing test of passivityassessment on the realization, which
    # finds the form active a part in a million under the lower bound and
    # passive as far over the ceiling, on the raw fit, turned or not, and
    # on an all pass, which touches one at every frequency.
    for (Rk, Dk) in ((R4, D4), (stack(Q4*R4[:, :, k] for k in axes(R4, 3)), Q4*D4))
        form = JC.ResidueForm(p4, Rk, Dk, sqrt(ws4[1]*ws4[end]), JC.residuespaces(p4, Rk, 1e-12))
        lower, _, ceiling = JC.residuenorm(form; rtol = 1e-8)
        A, B, C = JC.realization(form)
        verdictat(level) = JC.passivityassessment(A, B, C ./ level, form.D ./ level; atol = 0.0).verdict
        @test verdictat(lower*(1 - 1e-6)) === :active && verdictat(ceiling*(1 + 1e-6)) === :passive
        @test ceiling <= lower*(1 + 1e-8)
    end
    lower, _, ceiling = JC.residuenorm(allpass; rtol = 1e-8)
    @test lower <= 1 + 1e-15 && 1 <= ceiling <= 1 + 2e-8
    # A level, or a tolerance, within the roundoff of the norm leaves
    # centres which halving cannot settle, standing within the roundoff of
    # the level where the remainder over the interval is within it too:
    # they belong to a band, or count their bound toward the ceiling,
    # rather than halving on until the time runs out. The all pass, one at
    # every frequency, swept at 1 + 1e-15 and its norm searched to 1e-15 of
    # itself, whose ceiling is then looser than the roundoff, but found.
    tight = JC.passivitysweep(allpass, 1 + 1e-15, time_ns() + 1e10)
    @test tight.verdict === :indeterminate && tight.evaluated < 100
    tightnorm = JC.residuenorm(allpass; rtol = 1e-15, deadline = time_ns() + 1e10)
    @test tightnorm[1] ≈ 1 atol = 1e-15
    @test tightnorm[1] <= tightnorm[3] < Inf
    # and a fit the sweep cannot settle in the time allowed is refused:
    # a resonance of half width 0.01 peaking at 0.9, which the bound over the
    # whole axis leaves to the intervals; the norm search, out of time as
    # well, leaves it no ceiling, so that a contraction it measures is
    # refused rather than accepted on the intervals it reached
    @test_throws ArgumentError JC.enforcepassivity(ComplexF64[complex(-0.01, 1), complex(-0.01, -1)],
        fill(0.009 + 0im, 1, 1, 2), zeros(1, 1), [0.1, 1.0, 10.0], PassivityEnforcement(maxtime = 0.0))
    @test isinf(last(JC.residuenorm(JC.ResidueForm(ComplexF64[complex(-0.01, 1), complex(-0.01, -1)],
        fill(0.009 + 0im, 1, 1, 2), zeros(1, 1), 1.0, [ones(ComplexF64, 1, 1) for _ in 1:2]); rtol = 1e-8,
        deadline = 0)))
    # A resonance at 1e-16 of the scaled frequency, far under the samples,
    # peaking at 1.006: the norm search halves an interval down to the
    # floating-point resolution of its centre, so that its ceiling closes
    # on the peak, and a fit is accepted only where its ceiling stands
    # under the level, so that the fit returned stays under one there as
    # well. The norm alone, of 1.0001e-16/(s + 1e-16) beside a pole at -1
    # with no residue, is bracketed to rtol about its value.
    let p = ComplexF64[complex(-2.643763489844016e-17, 1e-16), complex(-2.643763489844016e-17, -1e-16), -1.0],
            R = reshape(ComplexF64[complex(2.556635034502508e-17, 1.1121039570931664e-17),
                complex(2.556635034502508e-17, -1.1121039570931664e-17), 1e-3], 1, 1, 3)
        Rk, Dk, _ = JC.enforcepassivity(p, R, fill(0.025214912183540387, 1, 1), [0.1, 1.0, 10.0])
        @test maximum(w -> abs(Dk[1, 1] + sum(Rk[1, 1, k]/(im*w - p[k]) for k in 1:3)),
            range(0.0, 3e-16; length = 3001)) <= 1 + 1e-8
        tiny = JC.ResidueForm(ComplexF64[-1e-16, -1.0], reshape(ComplexF64[1.0001e-16, 0.0], 1, 1, 2), zeros(1, 1), 1.0,
            [ones(ComplexF64, 1, 1) for _ in 1:2])
        lower, _, ceiling = JC.residuenorm(tiny; rtol = 1e-8)
        @test lower <= 1.0001 <= ceiling <= lower*(1 + 2e-8)
    end
    # A value stated at zero frequency is held exactly by the correction:
    # the change there is linear in the unknowns and is eliminated before
    # the inequalities are solved, so the answer is the least change which
    # holds it. This resonance needs correcting, and correcting it without
    # the statement moves its value at zero by nine parts in a hundred.
    Ae, Be, Ce = [-0.25 1.0; -1.0 -0.25], [1.0 0.0; 0.0 1.0], [0.35 0.0; 0.0 0.35]
    pe, Re = poleresidues(Ae, Be, Ce)
    wse = collect(range(0.05, 4.0; length = 160))
    S0e = real.(response(Ae, Be, Ce, zeros(2, 2), 0.0))
    @test JC.passivityassessment(Ae, Be, Ce, zeros(2, 2)).lower > 1.3
    Rh, Dh, _ = JC.enforcepassivity(pe, Re, zeros(2, 2), wse; dc = S0e)
    @test densemax(realized(pe, Rh, Dh)...) <= 1 + 2e-8
    @test real.(response(realized(pe, Rh, Dh)..., 0.0)) ≈ S0e atol = 1e-12
    Rf, Df, _ = JC.enforcepassivity(pe, Re, zeros(2, 2), wse)
    @test densemax(realized(pe, Rf, Df)...) <= 1 + 2e-8
    @test maximum(abs, real.(response(realized(pe, Rf, Df)..., 0.0)) .- S0e) > 0.05
    # A constant term within the tolerance above one is left to the rounds,
    # which correct it where it stands above their level, and to nothing
    # else: a series resonator to ground, an open at both ends of the axis
    # whose samples stand above one by half the tolerance, with a leakage
    # stating its value at zero frequency 2e-7 under one, fits at three
    # poles as closely as its raw fit, its constant at one to the fit's
    # roundoff and its statement held.
    let gs = collect(range(0.5e9, 10e9; length = 300)), Rl = 5e8,
            Sr = reshape([(1 + 5e-9)*(Z - 50)/(Z + 50) for Z in (5.0 .+ im .* 2pi .* gs .* 2e-9 .+
                1 ./ (im .* 2pi .* gs .* 0.5e-12))], 1, 1, :),
            open = ScatteringParameters((2pi .* gs, Sr); nports = 1, zref = 50.0,
                dcmodel = ScatteringDC(fill((Rl - 50)/(Rl + 50), 1, 1)))
        @test JC.relativefiterror(RationalScattering(open, 3), Sr, 2pi .* gs) <= 1e-6
    end
    # The rounds hold a fit to the dissipation its samples and constant term
    # are held to, `I - S'S` no lower than `-atol`: a constant term above
    # `sqrt(1 + atol)` but under one plus half of `atol`, whose samples meet
    # their tolerance, is corrected rather than refused by the constructor
    let atol = 1e-3, xs = collect(range(0.1, 10.0; length = 80)),
            Sd = reshape(ComplexF64[1 + atol/2 - atol^2/10 - 0.2/(1 + im*x) for x in xs], 1, 1, :)
        local fitd = RationalScattering(ScatteringParameters((xs, Sd); nports = 1, zref = 50.0, atol), 1; atol, tol = 0.01)
        @test JC.passivitymargin(fitd.provider.D) >= -atol
        @test JC.relativefiterror(fitd, Sd, xs) <= 1e-3
    end
    # A fit held at a value of unit norm at zero frequency, where the rounds
    # end uncertified: the all pass (s - 1)/(s + 1), held at -1 there, with a
    # bump of 1e-6 or 1e-7 at the frequency 3, where it stands at 0.8 + 0.6i,
    # enforced in one round. The norm search resolves under the margin to
    # the level, and the contraction toward the held value aims at the
    # level, since no point of its path has a ceiling under one: the fit
    # returned holds -1 at zero frequency and stays under 1 + atol.
    for delta in (1e-6, 1e-7)
        q = complex(-0.1, 3.0)
        r = delta*(0.8 + 0.6im)*0.1
        p = ComplexF64[-1.0, q, conj(q)]
        R = reshape(ComplexF64[-2.0, r, conj(r)], 1, 1, 3)
        at(Rk, Dk, w) = Dk[1, 1] + sum(Rk[1, 1, k]/(im*w - p[k]) for k in 1:3)
        D = fill(-1 - real(at(R, zeros(1, 1), 0.0)), 1, 1)
        Rk, Dk, _ = JC.enforcepassivity(p, R, D, collect(range(0.1, 10.0; length = 100)),
            PassivityEnforcement(rounds = 1); dc = fill(-1.0, 1, 1))
        @test real(at(Rk, Dk, 0.0)) ≈ -1 atol = 1e-12
        @test maximum(w -> abs(at(Rk, Dk, w)), range(0.0, 12.0; length = 12001)) <= 1 + 1e-8
    end
    # The rank of each residue is decided once, and the enforcement and the
    # realization take the residues within the same spaces, so that the
    # block realized is the one made passive. Two nearly equal poles whose
    # residues nearly cancel on the first port, `diag(100, -0.5)` and
    # `diag(-100, 0)`, and `diag(0, 1.2)` at the first: the `-0.5` is under
    # the `ranktol` of 0.01 of its residue, and without it the second port
    # is `1.2/(s + 1)`, active, where with it the port is `0.7/(s + 1)`.
    let pc = ComplexF64[-1, -1 - 1e-4, -1], Rc = zeros(ComplexF64, 2, 2, 3)
        Rc[:, :, 1] = Diagonal([100, -0.5])
        Rc[:, :, 2] = Diagonal([-100, 0])
        Rc[:, :, 3] = Diagonal([0, 1.2])
        spaces = JC.residuespaces(pc, Rc, 0.01)
        Rk, Dk, _ = JC.enforcepassivity(pc, Rc, zeros(2, 2), wse; spaces)
        A, B, C = JC.realization(pc, Rk, 2; spaces)
        @test size(A, 1) == 3
        @test densemax(A, B, C, Dk) <= 1 + 2e-8
    end
    # Where every residue has full rank the metric of the least change is
    # the basis's Gram matrix repeated for each unknown of a port's row, and
    # is factored as that; its solves agree with the whole normal matrix's,
    # a value held at zero or not.
    let pm = ComplexF64[-0.4, -0.2 + 2im, -0.2 - 2im], Rm = zeros(ComplexF64, 3, 3, 3)
        Rm[:, :, 1] = [1.0 0.3 0.1; 0.3 -0.7 0.2; 0.1 0.2 0.5]
        Rm[:, :, 2] = [0.5+0.1im 0.2im 0.1; 0.2im -0.3 0.4im; 0.1 0.4im 0.6-0.2im]
        Rm[:, :, 3] = conj.(Rm[:, :, 2])
        form = JC.ResidueForm(pm, Rm, zeros(3, 3), 1.0, JC.residuespaces(pm, Rm, 1e-12))
        blocks, nu = JC.correctioncoordinates(form)
        xm = exp.(range(log(0.01), log(100.0); length = 60))
        G = [sin(3i + j) for i in 1:nu, j in 1:5]
        for dcm in (JC.nodc, ones(3, 3))
            metric = JC.correctionmetric(form, blocks, nu, xm, 1e-8, dcm)
            @test metric.stride == 3
            along(f, X) = JC.alongbasis(f, X, metric.stride)
            Zf = isnothing(metric.Z) ? identity : X -> along(Y -> metric.Z*Y, X)
            Ztf = isnothing(metric.Z) ? identity : X -> along(Y -> transpose(metric.Z)*Y, X)
            factored = Zf(along(Y -> transpose(only(metric.L)) \ Y, along(Y -> only(metric.L) \ Y, Ztf(G))))
            H = JC.correctiongram(blocks, nu, pm, xm)
            H[diagind(H)] .+= 1e-8*tr(H)/nu
            Z = Matrix(1.0I, nu, nu)
            if !isempty(dcm)
                phi0 = JC.basisat!(zeros(ComplexF64, 3), pm, 0.0)
                E = zeros(3, nu)
                for (c, cols, M) in blocks
                    E[:, cols] .+= real(c == 0 ? 1.0 : phi0[c]) .* M
                end
                Z = nullspace(E)
            end
            @test factored ≈ Z*((transpose(Z)*H*Z) \ (transpose(Z)*G)) rtol = 1e-9
        end
        # Weighted, each entry has a Gram matrix of its own and each port a
        # normal matrix: the factored solves agree with the normal matrix
        # formed from the unknowns' responses at every sample, each squared
        # and weighed by its entry's weight there, for residues of full
        # rank, factored entry by entry, and of rank one, port by port.
        Wm = [1 + 0.5*sin(i + 2j + 0.3k) for i in 1:3, j in 1:3, k in eachindex(xm)]
        Rone = zeros(ComplexF64, 3, 3, 3)
        Rone[:, :, 1] = [1.0, 0.3, -0.5]*transpose([1.0, 0.3, -0.5])
        Rone[:, :, 2] = [0.5 + 0.2im, 0.1, 0.4im]*transpose([0.5 + 0.2im, 0.1, 0.4im])
        Rone[:, :, 3] = conj.(Rone[:, :, 2])
        for (Rw, stride) in ((Rm, 3), (Rone, 1)), dcm in (JC.nodc, ones(3, 3))
            formw = JC.ResidueForm(pm, Rw, zeros(3, 3), 1.0, JC.residuespaces(pm, Rw, 1e-12))
            blocksw, nuw = JC.correctioncoordinates(formw)
            metric = JC.correctionmetric(formw, blocksw, nuw, xm, 0.0, dcm, Wm)
            @test metric.stride == stride
            along(f, X) = JC.alongbasis(f, X, stride)
            Gw = [sin(3i + j) for i in 1:nuw, j in 1:5]
            Z = Matrix(1.0I, nuw, nuw)
            if !isempty(dcm)
                phi0 = JC.basisat!(zeros(ComplexF64, 3), pm, 0.0)
                E = zeros(3, nuw)
                for (c, cols, M) in blocksw
                    E[:, cols] .+= real(c == 0 ? 1.0 : phi0[c]) .* M
                end
                Z = nullspace(E)
            end
            Ztw = isnothing(metric.Z) ? Gw : along(Y -> transpose(metric.Z)*Y, Gw)
            g, phi = zeros(ComplexF64, nuw), zeros(ComplexF64, 3)
            for i in 1:3
                H = zeros(nuw, nuw)
                for (k, x) in enumerate(xm), j in 1:3
                    JC.coordinateresponse!(g, blocksw, JC.basisat!(phi, pm, x), Matrix(1.0I, 3, 3)[:, j])
                    H .+= Wm[i, j, k]^2 .* real.(conj.(g) .* transpose(g))
                end
                solved = JC.portsolve(metric, i, JC.portsolve(metric, i, Ztw); adjoint = true)
                isnothing(metric.Z) || (solved = along(Y -> metric.Z*Y, solved))
                @test solved ≈ Z*((transpose(Z)*H*Z) \ (transpose(Z)*Gw)) rtol = 1e-8
            end
        end
    end
    # and a value stated at zero frequency is held through the rank decided,
    # the constant taking back what a residue loses there, for a real pole
    # and a pair: `diag(0.5, 0.1)` at a `ranktol` of 0.3 keeps its first
    # direction alone
    for pd in (ComplexF64[-1], ComplexF64[-1 + 2im, -1 - 2im])
        Rd = zeros(ComplexF64, 2, 2, length(pd))
        Rd[:, :, 1] = Diagonal(length(pd) == 1 ? [0.5, 0.1] : [0.5 + 0.2im, 0.1 + 0.1im])
        length(pd) == 2 && (Rd[:, :, 2] = conj.(Rd[:, :, 1]))
        S0d = -real.(sum(Rd[:, :, k] ./ pd[k] for k in eachindex(pd)))
        spaces = JC.residuespaces(pd, Rd, 0.3)
        Dd = JC.dcthroughspaces!(zeros(2, 2), pd, Rd, spaces, S0d)
        Rk, Dk, _ = JC.enforcepassivity(pd, Rd, Dd, [0.1, 1.0, 10.0]; dc = S0d, spaces)
        A, B, C = JC.realization(pd, Rk, 2; spaces)
        @test size(A, 1) == length(pd)
        @test real.(response(A, B, C, Dk, 0.0)) ≈ S0d atol = 1e-14
    end
    # A feedthrough above one is brought under first, toward the value
    # stated at zero where there is one, which holds it: `S -> S0 + t (S -
    # S0)`. `-0.501/(s + 1) + 1.001` is anchored at 0.5; and an anchor of
    # unit norm is contracted toward, not refused: `-1.001 + 2.001/(s + 1)`
    # anchored at 1 goes at `t = 2/2.001` to the all pass
    # `(1 - s)/(1 + s)`, holding it.
    pa, Ra = poleresidues(fill(-1.0, 1, 1), ones(1, 1), fill(-0.501, 1, 1))
    Rg, Dg, _ = JC.enforcepassivity(pa, Ra, fill(1.001, 1, 1), [0.001, 0.01, 0.1]; dc = fill(0.5, 1, 1))
    @test densemax(realized(pa, Rg, Dg)...) <= 1 + 2e-8
    @test only(real.(response(realized(pa, Rg, Dg)..., 0.0))) ≈ 0.5 atol = 1e-12
    p1, R1 = poleresidues(fill(-1.0, 1, 1), ones(1, 1), fill(2.001, 1, 1))
    Ru, Du, _ = JC.enforcepassivity(p1, R1, fill(-1.001, 1, 1), [0.1, 1.0, 10.0]; dc = ones(1, 1))
    @test densemax(realized(p1, Ru, Du)...) <= 1 + 2e-8
    @test only(real.(response(realized(p1, Ru, Du)..., 0.0))) ≈ 1 atol = 1e-10
    # A feedthrough over one which the stated value cannot be contracted
    # toward: both are positive, so every step along the path moves away
    # from one, and there is no step. The same feedthrough with the
    # opposite sign is repaired exactly, as above.
    pz, Rz = poleresidues(fill(-1.0, 1, 1), ones(1, 1), fill(0.1, 1, 1))
    @test_throws ArgumentError JC.enforcepassivity(pz, Rz, fill(1.001, 1, 1), [0.1, 1.0]; dc = ones(1, 1))
    # A resonance whose peak stands 2e-8 above one, `k/(s^2 + 2 z s + 1)`
    # with `z = 3e-4`, is brought under where it peaks, with no
    # contraction to warn of; and a peak decades from both of its real
    # poles, `k s/((s + a)(s + b))` peaking at `sqrt(a b)`, is found and
    # brought under there.
    let z = 3e-4, kk = (1 + 2e-8)*2z*sqrt(1 - z^2), wpk = sqrt(1 - 2z^2)
        Ah, Bh, Ch = [0.0 1.0; -1.0 -2z], reshape([0.0, kk], 2, 1), reshape([1.0, 0.0], 1, 2)
        ph, Rh = poleresidues(Ah, Bh, Ch)
        Rp, Dp, contraction = @test_logs JC.enforcepassivity(ph, Rh, zeros(1, 1),
            2pi .* collect(range(0.01, 1.0; length = 200)))
        @test isnothing(contraction)
        @test abs(only(response(realized(ph, Rp, Dp)..., wpk))) <= 1
    end
    let aa = 2.7592661119815163e-5, bb = 1.0
        Ar2, Br2 = [0.0 1.0; -aa*bb -(aa + bb)], reshape([0.0, 1.0], 2, 1)
        Cr2 = reshape([0.0, (1 + 2e-8)*(aa + bb)], 1, 2)
        p2, R2 = poleresidues(Ar2, Br2, Cr2)
        Rq, Dq, _ = JC.enforcepassivity(p2, R2, zeros(1, 1), exp.(range(log(1e-4), log(1e2); length = 101)))
        @test abs(only(response(realized(p2, Rq, Dq)..., sqrt(aa*bb)))) <= 1
        @test JC.passivityassessment(realized(p2, Rq, Dq)...; atol = 0.0).lower <= 1
    end
    # and a violation narrower than the sweep's grid, a high Q resonance
    # above one, is corrected where it is, without a contraction, and
    # leaves the response away from it where it was
    for peak in (1.002, 1.05)
        Ar, Br = [0.0 1.0; -1.0 -2e-4], reshape([0.0, 1.0], 2, 1)
        Cr, Dr = reshape([0.0, (peak - 0.5)*2e-4], 1, 2), fill(0.5, 1, 1)
        pr, Rr = poleresidues(Ar, Br, Cr)
        R2, D2, contraction = JC.enforcepassivity(pr, Rr, Dr, collect(range(0.5, 1.5; length = 401)))
        @test isnothing(contraction)
        @test densemax(realized(pr, R2, D2)...) <= 1 + 2e-8
        @test maximum(abs(only(response(realized(pr, R2, D2)..., w)) - only(response(Ar, Br, Cr, Dr, w)))
            for w in (0.6, 0.9, 1.1, 1.4)) < 1e-4
    end
    # The dual of the least change factors its working set as the set
    # changes rather than afresh at every step: on `m` constraints of `2m`
    # unknowns, an identity beside a block coupling them, with bounds of
    # both signs, where the working set grows to most of the constraints,
    # it allocates in proportion to the Gram matrix, the square of the
    # constraints, rather than to the cube of the working set.
    function dualbytes(m)
        Ad = hcat(Matrix(1.0I, m, m), [sinpi(0.37*(i + 3j)) + 0.4*cospi(0.11*(2i + j)) for i in 1:m, j in 1:m])
        bd = [sinpi(0.29*k) - 0.7 for k in 1:m]
        Gd = Symmetric(Ad*transpose(Ad))
        @test last(JC.dualactiveset(Gd, bd))
        return @allocated JC.dualactiveset(Gd, bd)
    end
    @test dualbytes(200) < 6*dualbytes(100)
    # The multiplier which sets a step's length retires at zero, whatever
    # rounding error the step leaves it: a multiplier left a rounding error
    # above zero blocks the next step at its own size, and the steps after
    # it, retiring nothing, spend the sweep's budget and end it off the
    # optimum, which the enforcement takes for a round whose constraints
    # cannot be met. Two duals of `min x'x/2` subject to `A x <= b` on which
    # the sweep did so, found by a search over random feasible duals with
    # pairs of nearly parallel rows, as reciprocal constraints give, and
    # rows and bounds scaled over several decades: three constraints, where
    # the sweep alone ended 7e-2 of the bounds off the optimum, and nine,
    # where the coordinate descent then left the round unsolved.
    let A3 = [0.006936863441321314 0.00149191056612436 -0.013686011167687383;
              -0.0016790194863737934 -6.614127098434518e-5 -0.011355245454936443;
              0.0006620272682071061 -0.00023572090196369994 -0.00040351231639398536],
        b3 = [-1754.5968607294744, -104.89822891961165, -241.4954662993951],
        A9 = [1.0493191801220174 1.4429812564867417 -18.74657985188262;
              1.143224797253575 15.666616473491201 2.1309602039272337;
              -22.34197093717954 23.723821551264876 2.6515014958160736;
              -12.441486899904568 41.60713025718348 12.126429964726412;
              -7.434111908596136 -5.356332841825573 -2.8063830835154224;
              1.8574997681483516 -3.28808655425062 -1.2851001091614866;
              -22.346276395071932 23.728388123401217 2.6520213291784205;
              1.14391941418418 15.676240099294937 2.1322682193288207;
              -7.4389187257161975 -5.35979533906799 -2.808193964378653],
        b9 = [2.7466495000486013, -10.038576690349537, -11.059019365023284, -27.669235053121998,
              6.996117346830494, 2.0271661408220005, -10.991944514048773, -9.838564069621517, 6.112565963370436]
        # the conditions of the optimum, relative to the largest bound
        optimality(G, c, μ) = (g = G*μ .+ c; maximum(k -> μ[k] > 0 ? abs(g[k]) : max(0.0, -g[k]), eachindex(c))/maximum(abs, c))
        G3 = Symmetric(A3*transpose(A3))
        @test optimality(G3, b3, JC.lawsonhanson(G3, b3)) <= 1e-12
        μ9, ok9 = JC.dualactiveset(Symmetric(A9*transpose(A9)), b9)
        @test ok9 && maximum(A9*(-(transpose(A9)*μ9)) .- b9) <= 1e-12*maximum(abs, b9)
    end
end
