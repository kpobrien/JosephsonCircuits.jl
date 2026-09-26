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
    more = RationalScattering(data, 8; frequencies = fs)
    @test size(more.provider.A) == (3, 3)
    JC.evaluateprovider!(fit, more.provider, 2pi .* fs)
    @test maximum(abs.(fit .- hb.S)) < 1e-10
    raw = RationalScattering(data, 4; passivity = nothing)
    @test size(raw.provider.A) == (3, 3)
    @test densemax(raw.provider) <= 1 + 1e-10
    @test_throws ArgumentError RationalScattering(data, 4; frequencies = reverse(fs))
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
    @test_throws ArgumentError PassivityEnforcement(rtol = 0.0)
    @test_throws ArgumentError VectorFitting(pruneslack = -0.5)
    @test_throws ArgumentError VectorFitting(iterations = 0)
    @test_throws ArgumentError RationalScattering(data, 3; tol = -1.0)
    # A fixed order is held to its data as the search is, against a
    # tolerance of its own: two poles miss this RLC, which has three,
    # and are refused unless a tolerance that loose is asked for, the
    # fit returned then warning of the contraction it was made passive by
    @test_throws ArgumentError RationalScattering(data, 2)
    short = @test_logs (:warn,) RationalScattering(data, 2; tol = 1.0)
    @test 1e-2 < JC.relativefiterror(short, hb.S, fs) <= 1.0
    for kw in ((; passivity = PassivityEnforcement(margin = 1e-9, rounds = 40, scalelimit = 0.1)),
            (; passivity = PassivityEnforcement(density = 40, rtol = 1e-7)),
            (; fitting = VectorFitting(pruneslack = 0.0, settletol = 1e-12, stallpatience = 10)),
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
    # The passivity correction holds a stated value at zero frequency,
    # rather than imposing it in the residue solve and then walking
    # away from it. `S(0) = D + C X0` with `X0 = (-A)^-1 B`, so a
    # perturbation holds it exactly when `dC X0 + dD = 0`, which is
    # linear and is eliminated before the inequalities are solved; the
    # answer is then the least change which holds it and not the least
    # change which nearly does. This resonance needs correcting, and
    # correcting it without the condition moves its value at zero by
    # nine parts in a hundred.
    Ae = [-0.25 1.0; -1.0 -0.25]
    Be = [1.0 0.0; 0.0 1.0]
    Ce = [0.35 0.0; 0.0 0.35]
    De = [0.0 0.0; 0.0 0.0]
    wse = collect(range(0.05, 4.0; length = 160))
    S0e = real.(response(Ae, Be, Ce, De, 0.0))
    @test first(JC.hinfnorm(Ae, Be, Ce, De)) > 1.3
    @test opnorm(S0e) < 1
    # the enforcement returns the contraction its end applied, which the
    # caller warns of for a fit it returns
    held = JC.enforcepassivity(Ae, Be, Ce, De, wse; dc = S0e)
    @test !isnothing(held[5])
    @test real.(response(held[1:4]..., 0.0)) ≈ S0e atol=1e-12
    # passive to the norm search's own tolerance, twice its rtol, at
    # which the enforcement decides
    @test densemax(held[1:4]...) <= 1 + 2e-8
    free = JC.enforcepassivity(Ae, Be, Ce, De, wse)
    @test !isnothing(free[5])
    @test densemax(free[1:4]...) <= 1 + 2e-8
    @test maximum(abs, real.(response(free[1:4]..., 0.0)) .- S0e) > 0.05
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
    let Aa = fill(-1.0, 1, 1), Ba = ones(1, 1), Ca = fill(-0.501, 1, 1),
        Da = fill(1.001, 1, 1), S0a = fill(0.5, 1, 1)
        @test first(JC.hinfnorm(Aa, Ba, Ca, Da)) > 1
        got = JC.enforcepassivity(Aa, Ba, Ca, Da, [0.001, 0.01, 0.1]; dc = S0a)
        @test densemax(got[1:4]...) <= 1 + 2e-8
        @test only(real.(response(got[1:4]..., 0.0))) ≈ 0.5 atol=1e-12
        # and an anchor of unit norm is contracted toward, not
        # refused: `-1.001 + 2.001/(s + 1)` is anchored at
        # `S(0) = 1`, and `t = 2/2.001` takes it to
        # `(1 - s)/(1 + s)`, exactly all pass, holding the anchor
        got1 = JC.enforcepassivity(fill(-1.0, 1, 1), ones(1, 1), fill(2.001, 1, 1),
            fill(-1.001, 1, 1), [0.1, 1.0, 10.0]; dc = ones(1, 1))
        @test densemax(got1[1:4]...) <= 1 + 2e-8
        @test only(real.(response(got1[1:4]..., 0.0))) ≈ 1 atol=1e-10
    end
    # A violation narrower than the grid the enforcement sweeps, a high Q
    # resonance standing above one, is found where its pole is and
    # corrected there: the fit comes back passive with no contraction,
    # and away from the resonance its response is where it was, where a
    # sweep which saw nothing left a contraction of the whole block, or
    # a refusal
    Sat(A, B, C, D, w) = only(response(A, B, C, D, w))
    for peak in (1.002, 1.05)
        Ar, Br = [0.0 1.0; -1.0 -2e-4], reshape([0.0, 1.0], 2, 1)
        Cr, Dr = reshape([0.0, (peak - 0.5)*2e-4], 1, 2), fill(0.5, 1, 1)
        fixed = JC.enforcepassivity(Ar, Br, Cr, Dr, collect(range(0.5, 1.5; length = 401)))
        @test isnothing(fixed[5])
        @test densemax(fixed[1:4]...) <= 1 + 2e-8
        @test maximum(abs(Sat(fixed[1:4]..., w) - Sat(Ar, Br, Cr, Dr, w)) for w in (0.6, 0.9, 1.1, 1.4)) < 1e-4
    end
    # The correction is assembled one output port block at a time,
    # which rests on the normal matrix being the same block for
    # every port: restricted to port `i`'s own unknowns the
    # perturbation reads `[X(w)' I]`, and that does not depend on
    # `i`. The property itself is tested.
    let Ab = [-0.25 1.0; -1.0 -0.25], Bb = [1.0 0.0; 0.0 1.0], nzb = 2, nb = 2
        blocks = JC.portblocks(nb, nzb)
        @test length(blocks) == nb
        @test sort(vcat(blocks...)) == 1:(nb*nzb + nb*nb)   # a partition
        for w in (0.3, 1.0, 2.5)
            X = (im*w*I - Ab) \ Bb
            M = zeros(ComplexF64, nb*nb, nb*nzb + nb*nb)
            for i in 1:nb, j in 1:nb
                row = (i - 1)*nb + j
                for k in 1:nzb
                    M[row, (k - 1)*nb + i] = X[k, j]
                end
                M[row, nb*nzb + (i - 1)*nb + j] = 1.0
            end
            first = real.(transpose(M[:, blocks[1]])*M[:, blocks[1]])
            for i in 2:nb
                @test real.(transpose(M[:, blocks[i]])*M[:, blocks[i]]) ≈ first atol=1e-12
            end
        end
    end
    # and the correction it produces is feasible and makes the block
    # passive, on a resonance which needs one
    let Ad = [-0.25 1.0; -1.0 -0.25], Bd = [1.0 0.0; 0.0 1.0],
        Cd2 = [0.35 0.0; 0.0 0.35], Dd = [0.0 0.0; 0.0 0.0],
        wsd = collect(range(0.05, 4.0; length = 60))
        @test first(JC.hinfnorm(Ad, Bd, Cd2, Dd)) > 1
        got = JC.enforcepassivity(Ad, Bd, Cd2, Dd, wsd)
        @test densemax(got[1:4]...) <= 1 + 2e-8
        S0d = real.(response(Ad, Bd, Cd2, Dd, 0.0))
        held = JC.enforcepassivity(Ad, Bd, Cd2, Dd, wsd; dc = S0d)
        @test densemax(held[1:4]...) <= 1 + 2e-8
        @test real.(response(held[1:4]..., 0.0)) ≈ S0d atol=1e-12
    end
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
    # A tabulated block is evaluated at the angular frequencies it
    # stores, not at those divided by `2 pi` and multiplied by it
    # again: that round trip puts this table's first sample at
    # 0.19999999999999998, one unit in the last place outside the
    # table's own range.
    let wt = [0.2, 1.0]
        @test 2pi*(wt[1]/(2pi)) != wt[1]
        St = reshape(ComplexF64[0.3/(1 + im*w) for w in wt], 1, 1, :)
        bt = ScatteringParameters((wt, St); nports = 1, zref = 50.0)
        ft = RationalScattering(bt, 1)
        @test size(ft.provider.A, 1) == 1
    end
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
        @test_throws ArgumentError RationalScattering(blkz, 1; frequencies = bad)
    end
    # `vectorfit` is the fit and not the enforcement. Samples of
    # `1.5 - 1.2/(s + 1)` over a low band have an exact one pole fit
    # whose constant term is 1.5, and the fit returns it: bringing a
    # constant term under one is the enforcement's work, since no
    # perturbation over a band reaches infinite frequency, and it is
    # done for the caller who asks for it rather than to every fit.
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
    # The order can be searched for instead of given: the fewest poles
    # which hold the error over the samples under a tolerance, as a
    # fraction of the largest response. This RLC is exact at three
    # poles, so every tolerance it can meet it meets there, and the
    # search returns three however loose or tight the tolerance is.
    auto = RationalScattering(data; tol = 1e-6)
    @test size(auto.provider.A) == (3, 3)
    JC.evaluateprovider!(fit, auto.provider, 2pi .* fs)
    @test maximum(abs.(fit .- hb.S)) < 1e-10
    @test JC.relativefiterror(auto, hb.S, fs) <= 1e-6
    # the tolerance bounds what the caller receives: a looser one may
    # be met by fewer poles but never by a fit which misses it
    for tol in (1e-2, 1e-8)
        got = RationalScattering(data; tol = tol)
        @test JC.relativefiterror(got, hb.S, fs) <= tol
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
        byhand = maximum(eachindex(xq)) do k
            opnorm(Dq .+ sum(res[:, :, q] ./ (im*xq[k] - ps[q]) for q in eachindex(ps)) .-
                   view(Sq, :, :, k))
        end
        @test JC.fiterror(Sq, xq, ps) ≈ byhand
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
            @test JC.fiterror(Sdelay, xd, kept) <= (1 + slack)*base + JC.roundoff(Sdelay)
        end
    end
    reach = np -> try
        JC.relativefiterror(RationalScattering(delayed, np; tol = 1e3), Sdelay, gs)
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
        @test JC.relativefiterror(got, Sdelay, gs) <= tol
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
    @test JC.relativefiterror(expanded, Sdelay, gs) <= 1e-8
    # a `maxpoles` given by the caller is a wall, because the caller
    # made it one; a search which returns nothing warns of nothing
    @test_logs @test_throws ArgumentError RationalScattering(
        delayed; tol = 1e-8, noisefloor = 1e-2, maxpoles = 4)
    # and the expansion stops rather than running to the sample count:
    # a tolerance nothing reaches is still reported
    @test_logs @test_throws ArgumentError RationalScattering(
        delayed; tol = 1e-16, maxpoles = 8)
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
    # a tolerance the data cannot support names the degree it carries
    cannot = try
        RationalScattering(data; tol = 1e-16); ""
    catch e
        sprint(showerror, e)
    end
    @test occursin("determine a degree of about", cannot)
    # Below the degree the samples determine, an order which produces
    # no fit has too few poles, and the search goes on past it: six
    # coupled resonators between two ports, whose samples determine a
    # degree of about twelve and whose orders seven to eleven are all
    # refused, are fitted at thirteen.
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
        @test JC.relativefiterror(RationalScattering(chain; tol = 1e-6), Sc, fc) <= 1e-6
        # A search which misses names the closest fit and why the orders
        # which could not be fitted failed, since that is a different
        # problem from a tolerance too tight; one in which no order
        # produced a fit says that, and has no closest fit to name.
        message(; kw...) = try
            RationalScattering(chain; tol = 1e-6, kw...); ""
        catch e
            e isa ArgumentError || rethrow()
            sprint(showerror, e)
        end
        missed = message(minpoles = 5, maxpoles = 12)
        @test occursin("closest was", missed) && occursin("could not be fitted at all", missed)
        none = message(minpoles = 7, maxpoles = 11)
        @test occursin("could be fitted at all", none) && !occursin("closest", none)
    end
    # A block may be active, at zero frequency as at any other, so long
    # as it declares its own noise; without one an active statement is
    # refused by the block itself, before any fitting, and an active
    # realization by the block's constructor.
    @test_throws ArgumentError ScatteringParameters((2pi .* fs, hb.S); nports = 2,
        zref = 50.0, dcmodel = JC.ScatteringDC([1.4 0.0; 0.0 1.4]))
    @test_throws ArgumentError RationalScattering(zeros(0, 0), zeros(0, 2), zeros(2, 0), [0.0 0.0; 3.0 0.0]; zref = 50.0)
    @test RationalScattering(zeros(0, 0), zeros(0, 2), zeros(2, 0), [0.0 0.0; 3.0 0.0]; zref = 50.0,
        noise = JC.NoiseCovariance([1.0 0.0; 0.0 8.0])) isa ScatteringParameters
    # An amplifier whose gain rolls off, stating the noise a quantum
    # limited one has, |I - S S'|, and its active value at zero
    # frequency, is fitted as it is, without the passivity enforcement,
    # which could not meet the statement, and without the validation:
    # the fit is stable, meets the statement exactly, has the gain of
    # the samples, and carries the stated noise, whichever way it is
    # asked for.
    w0, g0 = 2pi*8e9, 10.0
    amp(w) = [0.0 0.0; g0*w0/(w0 + im*w) 0.0]
    ampnoise(w) = (K = I - amp(w)*amp(w)'; [abs(K[1, 1]) 0.0; 0.0 abs(K[2, 2])])
    ampdata = ScatteringParameters(amp; nports = 2, zref = 50.0, noise = JC.NoiseCovariance(ampnoise),
        dcmodel = JC.ScatteringDC([0.0 0.0; g0 0.0]))
    ampfs = collect(range(0.1e9, 40e9; length = 300))
    for ampfit in (RationalScattering(ampdata, 4; frequencies = ampfs),
            RationalScattering(ampdata; tol = 1e-6, minpoles = 1, maxpoles = 6, frequencies = ampfs),
            RationalScattering(ampdata, 4; frequencies = ampfs, passivity = nothing))
        @test ampfit.noise isa JC.NoiseCovariance
        @test size(ampfit.provider.A) == (1, 1) && maximum(real.(eigvals(ampfit.provider.A))) < 0
        Sa = zeros(ComplexF64, 2, 2, 4)
        JC.evaluateprovider!(Sa, ampfit.provider, 2pi .* [0.0, 1e9, 5e9, 20e9])
        @test Sa ≈ cat([amp(2pi*f) for f in (0.0, 1e9, 5e9, 20e9)]...; dims = 3) rtol=1e-8
        @test abs(Sa[2, 1, 2]) > 1
        @test_throws ArgumentError JC.checkpassive(ampfit.provider)
    end
    # A feedthrough over one which the stated value cannot be
    # contracted toward: both are positive here, so every step along
    # the path moves away from one rather than toward it, and there is
    # no step. The same feedthrough with the opposite sign is repaired
    # exactly, which is the case above.
    @test_throws ArgumentError JC.enforcepassivity(fill(-1.0, 1, 1), ones(1, 1),
        fill(0.1, 1, 1), fill(1.001, 1, 1), [0.1, 1.0]; dc = ones(1, 1))
    # The constant term is the model at infinite frequency, which no
    # perturbation over a band of frequencies can reach, so a fit whose
    # constant term is active cannot be enforced passive at all. It is
    # brought under one before the residues are fitted, and the
    # residues then fitted to what is left.
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
    # axis and not merely at the samples, though to the accuracy of
    # the norm search and not exactly: it decides on the largest
    # singular value that search evaluated, which is a lower bound,
    # and the level the search ends at stands `2 rtol` above it. The
    # analytic case among the norm tests above pins that gap. A
    # block whose largest singular value exceeds one is active, and a
    # transient built on it grows without bound, so a nearly passive
    # fit is scaled the rest of the way rather than returned as it is.
    # The near lossless case is the one which needs it: there the
    # enforcement's own level of `1 + atol/2` sits above the block's
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
            # and the level, which is what termination establishes, is
            # over one by no more than the search's own tolerance
            @test JC.hinfnorm(q.A, q.B, q.C, q.D)[3] <= 1 + 4e-8
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
        f1 = RationalScattering(onepole, np; frequencies = range(0.01, 1.0; length = 41))
        @test size(f1.provider.A) == (1, 1)
        @test f1.provider.A[1, 1] ≈ -1.0 rtol=1e-8
    end
    pair = RationalScattering([-0.1 1.0; -1.0 -0.1], [0.0; 1.0;;], [0.05 0.0], zeros(1, 1))
    for np in (2, 3, 5)
        f2 = RationalScattering(pair, np; frequencies = range(0.05, 3.0; length = 121))
        @test size(f2.provider.A) == (2, 2)
        @test sort(eigvals(f2.provider.A); by = imag) ≈ [-0.1 - im, -0.1 + im] rtol=1e-8
    end
    @test_throws ArgumentError RationalScattering(ScatteringParameters(fill(0.3, 1, 1)), 2; frequencies = range(0.1, 1.0; length = 21))
    # a three port, whose nine entries are eliminated one by one in the
    # fit, so nothing the size of every entry's every sample is formed
    tee = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:r1, 1, 4, Resistor(20.0)), (:p2, 2, 0, Port(2; Z0 = 50.0)), (:r2, 2, 4, Resistor(20.0)),
        (:p3, 3, 0, Port(3; Z0 = 50.0)), (:l3, 3, 5, Inductor(2e-9)), (:c3, 5, 4, Capacitor(0.4e-12)), (:r4, 4, 0, Resistor(80.0)), (:c4, 4, 0, Capacitor(0.2e-12))])
    hb3 = hblinsolve(2pi .* fs, tee; keyedarrays = false)
    fit3 = RationalScattering(ScatteringParameters((2pi .* fs, hb3.S); nports = 3, zref = 50.0), 8)
    S3 = zeros(ComplexF64, 3, 3, length(fs))
    JC.evaluateprovider!(S3, fit3.provider, 2pi .* fs)
    @test maximum(abs.(S3 .- hb3.S)) < 1e-10 && size(fit3.provider.A, 1) <= 6
    @test_throws ArgumentError RationalScattering(data, 4; frequencies = fs[1:4])
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
    fitted = RationalScattering(block, 4; frequencies = fs, delays = taus)
    # what is left of the cable is the one pole it really is
    @test size(fitted.provider.A, 1) == 2
    # and putting the delay back reproduces the cable: the entry from
    # port q to port p carries the delay of both
    F = zeros(ComplexF64, 2, 2, length(fs))
    JC.evaluateprovider!(F, fitted.provider, 2pi .* fs)
    @test maximum(abs(F[p, q, i]*cis(-2pi*fs[i]*(taus[p] + taus[q])) - cable(2pi*fs[i])[p, q])
        for i in eachindex(fs), q in 1:2, p in 1:2) < 1e-12
    # the same cable fitted whole cannot be done at this order at all
    @test_throws ArgumentError RationalScattering(block, 4; frequencies = fs)
    # the line of the delay in cascade with the fit is the cable
    cascade = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:la, 1, 2, TransmissionLine(50.0, len; vp = vp)),
        (:b, 2, 3, fitted), (:lb, 3, 4, TransmissionLine(50.0, len; vp = vp)),
        (:p2, 4, 0, Port(2; Z0 = 50.0))])
    check = collect(range(1e9, 9e9; length = 12))
    S = hblinsolve(2pi .* check, cascade; keyedarrays = false).S
    @test maximum(abs(S[2, 1, i] - cable(2pi*check[i])[2, 1]) for i in eachindex(check)) < 1e-9
    # the order search takes them too, and reports the order of what is
    # left rather than of the delay
    searched = RationalScattering(block; tol = 1e-8, minpoles = 1, frequencies = fs, delays = taus)
    @test size(searched.provider.A, 1) == 2
    # no delay and every delay zero are the same fit, on the cable's
    # tail, which needs no delay taken out to be fitted at all
    tail(w) = (g = 1/(1 + im*w/wc); ComplexF64[0 g; g 0])
    tailblock = ScatteringParameters(tail; nports = 2, zref = 50.0)
    plain = RationalScattering(tailblock, 4; frequencies = fs, delays = (0.0, 0.0))
    bare = RationalScattering(tailblock, 4; frequencies = fs, delays = nothing)
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
            noise = NoiseCovariance(Diagonal(fill(2.0, 3)))), 4; frequencies = fs3, delays = taus3)
        V3 = zeros(ComplexF64, 3, 3, 2)
        JC.evaluatecovariance!(V3, uncorrelated, ws3[end - 1:end])
        @test all(V3[:, :, i] ≈ 2I for i in 1:2)
    end
    @test_throws ArgumentError RationalScattering(block, 4; frequencies = fs, delays = (tau,))
    @test_throws ArgumentError RationalScattering(block, 4; frequencies = fs, delays = (tau, -1.0))
    # a covariance the delays rotate is rotated; one they do not is not.
    # The cable's gain reaches one at zero frequency, so a covariance
    # which is to be one the block can emit over the whole band needs
    # its diagonal above what the commutator asks there
    C = [2.0 0.5; 0.5 2.0]
    stated = ScatteringParameters(cable; nports = 2, zref = 50.0,
        noise = NoiseCovariance((2pi .* fs, repeat(complex(C), 1, 1, length(fs)))))
    even = RationalScattering(stated, 4; frequencies = fs, delays = taus, passivity = nothing)
    V = zeros(ComplexF64, 2, 2, 2)
    JC.evaluateprovider!(V, even.noise.provider, 2pi .* fs[2:3])
    @test maximum(abs.(V .- repeat(complex(C), 1, 1, 2))) < 1e-14
    uneven = RationalScattering(stated, 4; frequencies = fs, delays = (tau, 0.0), passivity = nothing)
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
        fit = RationalScattering(held, 4; frequencies = fs, delays = (tau, 0.0), passivity = nothing)
        JC.evaluateprovider!(V, fit.noise.provider, between)
        @test all(abs(V[1, 2, i] - C[1, 2]*cis(between[i]*tau)) < 1e-14 for i in 1:2)
    end
    # and beyond a table it keeps turning with the frequency, being a
    # phase on the covariance the extrapolation states there
    narrow = ScatteringParameters(cable; nports = 2, zref = 50.0,
        noise = NoiseCovariance((2pi .* fs[1:20], repeat(complex(C), 1, 1, 20)); extrapolation = :constant))
    beyond = RationalScattering(narrow, 4; frequencies = fs, delays = (tau, 0.0), passivity = nothing)
    out = 2pi .* [fs[60], fs[80]]
    JC.evaluateprovider!(V, beyond.noise.provider, out)
    @test all(abs(V[1, 2, i] - C[1, 2]*cis(out[i]*tau)) < 1e-14 for i in 1:2)
    # a turned covariance is still data the commutation relations can be
    # checked against: this block has a power gain of four and emits the
    # fitted block's noise, which is far less than it must
    @test_throws ArgumentError ScatteringParameters(ComplexF64[0 2; 2 0]; zref = 50.0,
        noise = uneven.noise)
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
    f = RationalScattering(b, 2; frequencies = collect(range(0.1e9, 0.9e9; length = 10)), tol = 10.0, noisetol = 10.0)
    H = zeros(ComplexF64, 1, 1, 1)
    JosephsonCircuits.evaluateprovider!(H, f.providers[1], [0.0])
    @test H[1] ≈ 1
    @test f.dcmodel isa OpenDC
    @test f.atol == b.atol && f.noise.completed && f.noise.atol == b.noise.atol
    # a band which leaves the unconverted response no sample is refused
    lossless = LinearizedScattering([w -> fill((wp - im*w)/(wp + im*w), 1, 1), zero1], wp; harmonics = [0, 1], nports = 1, dcmodel = OpenDC())
    @test_throws ArgumentError RationalScattering(lossless, 2; frequencies = [0.1e9, 0.2e9, 0.3e9], band = (0.5e9, 0.9e9))
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
