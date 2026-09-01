isdefined(Main, :testjpacircuit) || include(joinpath(@__DIR__, "testcircuits.jl"))
using JosephsonCircuits
using LinearAlgebra
using SparseArrays
using Random
using Test

const JC = JosephsonCircuits
hbnlp_relerr(a, b) = norm(a - b)/max(norm(a), norm(b), eps())

# Three configurations covering the cases which behave differently. The
# dc+even case puts a self-conjugate mode in the layout; the two tone case
# is the only one with FCONJ conjugate symmetry slots, so the sign in the
# transposed forward map is exercised nowhere else.
function hbnlp_testproblems()
    circuit, defs = testjpacircuitnumeric()
    return [
        ("1 tone", JC.hbnonlinearproblem((2*pi*4.75001e9,), (8,),
            [(mode=(1,),port=1,current=0.02e-6)], circuit, defs)),
        ("1 tone dc+even", JC.hbnonlinearproblem((2*pi*4.75001e9,), (6,),
            [(mode=(1,),port=1,current=0.02e-6)], circuit, defs;
            dc=true, even=true)),
        ("2 tone", JC.hbnonlinearproblem((2*pi*4.65e9, 2*pi*4.85e9), (3,3),
            [(mode=(1,0),port=1,current=1e-8),(mode=(0,1),port=1,current=1e-8)],
            circuit, defs)),
    ]
end

@testset verbose=true "hbnonlinearproblem" begin

@testset "residual and jacobian-vector product" begin
    Random.seed!(20260828)
    for (nm, prob) in hbnlp_testproblems()
        @testset "$nm" begin
            n = length(prob)
            u = 0.05 .* randn(n)
            Jr = copy(prob.jacobian); JC.hbjacobian!(Jr, prob, u)
            for _ in 1:3
                v = randn(n)
                @test hbnlp_relerr(JC.hbjvp!(zeros(n), prob, u, v), Jr*v) < 1e-12
            end
            # against central differences of the residual
            v = randn(n); v ./= norm(v); h = 1e-6*max(norm(u), 1)
            fd = (JC.hbresidual!(zeros(n), prob, u .+ h.*v) .-
                  JC.hbresidual!(zeros(n), prob, u .- h.*v))./(2h)
            @test hbnlp_relerr(JC.hbjvp!(zeros(n), prob, u, v), fd) < 1e-5
        end
    end
end

@testset "transposed product" begin
    Random.seed!(20260828)
    for (nm, prob) in hbnlp_testproblems()
        @testset "$nm" begin
            n = length(prob)
            u = 0.05 .* randn(n)
            Jr = copy(prob.jacobian); JC.hbjacobian!(Jr, prob, u)
            for _ in 1:3
                v = randn(n); w = randn(n)
                Jtw = JC.hbvjp!(zeros(n), prob, u, w)
                @test hbnlp_relerr(Jtw, transpose(Jr)*w) < 1e-10
                # the adjoint identity, which catches a transpose that is
                # self consistent but wrong
                Jv = JC.hbjvp!(zeros(n), prob, u, v)
                @test isapprox(dot(w, Jv), dot(Jtw, v); rtol = 1e-10)
            end
        end
    end
end

@testset "transposed product: the conjugate multiplicity is load bearing" begin
    # The adjoint of the unnormalized inverse real transform is not the
    # forward transform: a stored bin whose conjugate partner is not also
    # stored represents two harmonics and is counted twice. Getting this
    # wrong gives a vjp which is right on self-conjugate modes and wrong by
    # a factor of two elsewhere, so the test must fail when it is removed.
    Random.seed!(20260828)
    _, prob = hbnlp_testproblems()[3]          # two tone, has FCONJ slots
    n = length(prob)
    u = 0.05 .* randn(n); w = randn(n)
    Jr = copy(prob.jacobian); JC.hbjacobian!(Jr, prob, u)
    ref = transpose(Jr)*w
    @test hbnlp_relerr(JC.hbvjp!(zeros(n), prob, u, w), ref) < 1e-10

    tp = JC.transposeplan!(prob)
    saved = copy(tp.gtscale)
    fill!(tp.gtscale, 1.0)
    @test hbnlp_relerr(JC.hbvjp!(zeros(n), prob, u, w), ref) > 1e-3
    copyto!(tp.gtscale, saved)
    @test hbnlp_relerr(JC.hbvjp!(zeros(n), prob, u, w), ref) < 1e-10
end

@testset "second and third derivatives" begin
    Random.seed!(20260828)
    for (nm, prob) in hbnlp_testproblems()
        @testset "$nm" begin
            n = length(prob)
            u = 0.05 .* randn(n)
            v = randn(n); w = randn(n); z = randn(n)./sqrt(n)
            # symmetric in its directions
            @test hbnlp_relerr(JC.hbd2F!(zeros(n), prob, u, v, w),
                         JC.hbd2F!(zeros(n), prob, u, w, v)) < 1e-12
            a = JC.hbd3F!(zeros(n), prob, u, v, w, z)
            @test hbnlp_relerr(a, JC.hbd3F!(zeros(n), prob, u, z, v, w)) < 1e-12
            @test hbnlp_relerr(a, JC.hbd3F!(zeros(n), prob, u, w, z, v)) < 1e-12
            # d2F against central differences of the product, d3F of d2F
            zh = z./norm(z); h = 1e-6*max(norm(u), 1)
            fd2 = (JC.hbjvp!(zeros(n), prob, u .+ h.*zh, v) .-
                   JC.hbjvp!(zeros(n), prob, u .- h.*zh, v))./(2h)
            @test hbnlp_relerr(JC.hbd2F!(zeros(n), prob, u, v, zh), fd2) < 1e-5
            fd3 = (JC.hbd2F!(zeros(n), prob, u .+ h.*zh, v, w) .-
                   JC.hbd2F!(zeros(n), prob, u .- h.*zh, v, w))./(2h)
            @test hbnlp_relerr(JC.hbd3F!(zeros(n), prob, u, v, w, zh), fd3) < 1e-4
        end
    end
end

@testset "drive scaling and its derivative" begin
    Random.seed!(20260828)
    for (nm, prob) in hbnlp_testproblems()
        @testset "$nm" begin
            n = length(prob)
            u = 0.05 .* randn(n)
            # dF/ds is exactly -b, so it does not depend on u
            dp = JC.hbdFdp!(zeros(n), prob)
            hs = 1e-6
            fd = (JC.drivenresidual!(zeros(n), prob, u, 1.0+hs) .-
                  JC.drivenresidual!(zeros(n), prob, u, 1.0-hs))./(2hs)
            @test hbnlp_relerr(dp, fd) < 1e-6
            @test hbnlp_relerr(dp, JC.hbdFdp!(zeros(n), prob)) < 1e-14
            # zero drive: the residual is purely the undriven system, so the
            # zero state is a solution
            JC.setdrive!(prob, 0)
            @test norm(JC.hbresidual!(zeros(n), prob, zeros(n))) < 1e-14
            # and the scale is relative, so 1 restores what was built
            JC.setdrive!(prob, 1)
            F1 = JC.hbresidual!(zeros(n), prob, u)
            JC.setdrive!(prob, 0.5); JC.setdrive!(prob, 1)
            @test hbnlp_relerr(JC.hbresidual!(zeros(n), prob, u), F1) < 1e-14
            # the residual is affine in the drive scale
            JC.setdrive!(prob, 1)
            Fa = JC.hbresidual!(zeros(n), prob, u)
            JC.setdrive!(prob, 3)
            Fb = JC.hbresidual!(zeros(n), prob, u)
            @test hbnlp_relerr(Fb .- Fa, 2 .* dp) < 1e-12
            JC.setdrive!(prob, 1)
        end
    end
end

@testset "JacobianOperator" begin
    Random.seed!(20260828)
    _, prob = hbnlp_testproblems()[1]
    n = length(prob)
    u = 0.05 .* randn(n); v = randn(n); w = randn(n)
    Jr = copy(prob.jacobian); JC.hbjacobian!(Jr, prob, u)
    J = JC.JacobianOperator(prob, u)
    @test size(J) == (n, n)
    @test size(J, 1) == n
    @test eltype(J) == Float64
    @test hbnlp_relerr(J*v, Jr*v) < 1e-12
    @test hbnlp_relerr(mul!(zeros(n), J, v), Jr*v) < 1e-12
    @test hbnlp_relerr(J'*w, transpose(Jr)*w) < 1e-10
    @test hbnlp_relerr(transpose(J)*w, transpose(Jr)*w) < 1e-10
    # the five argument form
    y = ones(n)
    mul!(y, J, v, 2.0, 3.0)
    @test hbnlp_relerr(y, 2 .* (Jr*v) .+ 3 .* ones(n)) < 1e-12
    @test JC.jacobianprototype(prob) !== prob.jacobian
    @test nnz(JC.jacobianprototype(prob)) == nnz(prob.jacobian)
end

@testset "matrix-free construction skips the Jacobian" begin
    circuit, defs = testjpacircuitnumeric()
    args = ((2*pi*4.75001e9,), (8,), [(mode=(1,),port=1,current=0.02e-6)],
            circuit, defs)
    pf = JC.hbnonlinearproblem(args...; assemblejacobian = false)
    pa = JC.hbnonlinearproblem(args...; assemblejacobian = true)
    @test isnothing(pf.jacobian)
    @test !isnothing(pa.jacobian)
    @test isnothing(JC.jacobianprototype(pf))
    # the products are the same either way
    Random.seed!(1); n = length(pf); u = 0.05 .* randn(n); v = randn(n)
    @test hbnlp_relerr(JC.hbjvp!(zeros(n), pf, u, v),
                 JC.hbjvp!(zeros(n), pa, u, v)) < 1e-14
    # the transposed product needs no assembled Jacobian either
    w = randn(n)
    @test hbnlp_relerr(JC.hbvjp!(zeros(n), pf, u, w),
                 JC.hbvjp!(zeros(n), pa, u, w)) < 1e-14
end

@testset "the products do not allocate" begin
    # The matrix-free products run inside the innermost Krylov loop, so an
    # allocation there is paid once per Jacobian application. `@allocated`
    # inside a testset closure measures boxing of captured variables rather
    # than the callee, so the measurement has to happen inside a function.
    function jvpallocs(prob, u, v)
        y = zeros(length(u)); JC.hbjvp!(y, prob, u, v)
        return @allocated JC.hbjvp!(y, prob, u, v)
    end
    function vjpallocs(prob, u, w)
        y = zeros(length(u)); JC.hbvjp!(y, prob, u, w)
        return @allocated JC.hbvjp!(y, prob, u, w)
    end
    Random.seed!(20260828)
    for (nm, prob) in hbnlp_testproblems()
        @testset "$nm" begin
            n = length(prob)
            u = 0.05 .* randn(n); v = randn(n); w = randn(n)
            @test jvpallocs(prob, u, v) == 0
            @test vjpallocs(prob, u, w) == 0
        end
    end
end

@testset "sparsity is a function of the structure, not the state" begin
    # Everything downstream -- assembly plans, cached symbolic
    # factorizations, preallocated workspaces -- assumes the pattern of the
    # assembled Jacobian is fixed. If it ever moves, those are silently
    # invalid rather than wrong in a way that throws.
    Random.seed!(20260828)
    for (nm, prob) in hbnlp_testproblems()
        @testset "$nm" begin
            n = length(prob)
            J = copy(prob.jacobian)
            JC.hbjacobian!(J, prob, zeros(n))
            colptr = copy(J.colptr); rowval = copy(J.rowval); nz = nnz(J)
            for scale in (1e-9, 1e-3, 0.5)
                JC.hbjacobian!(J, prob, scale .* randn(n))
                @test J.colptr == colptr
                @test J.rowval == rowval
                @test nnz(J) == nz
            end
            # and under a change of drive, which only touches the source
            for s in (0.0, 0.5, 3.0)
                JC.setdrive!(prob, s)
                JC.hbjacobian!(J, prob, 0.05 .* randn(n))
                @test J.colptr == colptr
                @test nnz(J) == nz
            end
            JC.setdrive!(prob, 1)
        end
    end
end

@testset "preconditioner calling conventions" begin
    # The ecosystem disagrees: LinearSolve.jl calls the three argument
    # `ldiv!`, IterativeSolvers.jl the two argument in-place form (on a
    # view), Krylov.jl and KrylovKit `mul!`. All must agree.
    Random.seed!(20260828)
    _, prob = hbnlp_testproblems()[1]
    n = length(prob)
    u = zeros(n)
    P = JC.preconditioner(prob, u)
    @test size(P) == (n, n)
    @test eltype(P) == Float64
    r = randn(n)
    z3 = ldiv!(zeros(n), P, r)
    @test hbnlp_relerr(copy(r) |> x -> (ldiv!(P, x); x), z3) < 1e-14
    @test hbnlp_relerr(mul!(zeros(n), P, r), z3) < 1e-14
    @test hbnlp_relerr(P \ r, z3) < 1e-14
    # in place on a view, which is what IterativeSolvers.jl hands it
    M = zeros(n, 2); M[:,1] .= r
    ldiv!(P, view(M, :, 1))
    @test hbnlp_relerr(M[:,1], z3) < 1e-14
    # a preconditioner which did nothing would pass everything above, so
    # check it actually improves the conditioning
    Jr = copy(prob.jacobian); JC.hbjacobian!(Jr, prob, u)
    PJ = reduce(hcat, [ldiv!(zeros(n), P, Vector(Jr[:,j])) for j in 1:n])
    @test cond(PJ) < cond(Matrix(Jr))
end

@testset "solver protocol" begin
    circuit, defs = testjpacircuitnumeric()
    wp = (2*pi*4.75001e9,); src = [(mode=(1,),port=1,current=0.00565e-6)]
    @test JC.solvermethod(JC.NewtonKrylov()) == :newtonkrylov
    @test JC.solvermethod(JC.Newton()) == :newton
    @test JC.solvermethod(JC.QuasiNewton()) == :quasinewton
    @test JC.solvermethod(:newtonkrylov) == :newtonkrylov

    ref = JC.hbnlsolve(wp, (8,), src, circuit, defs;
                       keyedarrays=false, ftol=1e-14)
    @test ref.solverinfo.converged
    for m in (JC.NewtonKrylov(), JC.Newton(), :newtonkrylov, :newton)
        s = JC.hbnlsolve(wp, (8,), src, circuit, defs;
                         keyedarrays=false, ftol=1e-14, method=m)
        @test s.solverinfo.converged
        @test isapprox(maximum(abs.(s.nodeflux)),
                       maximum(abs.(ref.nodeflux)); rtol=1e-8)
    end

    # a caller supplied solver, using only the public substrate
    ext = JC.ExternalSolver() do prob, u0
        u = copy(u0); F = similar(u)
        JC.hbresidual!(F, prob, u); r0 = norm(F)
        Jm = copy(prob.jacobian)
        for k in 1:60
            norm(F) <= 1e-13*max(r0,1) && return (u, true)
            JC.hbjacobian!(Jm, prob, u)
            u .-= Jm \ F
            JC.hbresidual!(F, prob, u)
        end
        return (u, norm(F) <= 1e-13*max(r0,1))
    end
    s = JC.hbnlsolve(wp, (8,), src, circuit, defs;
                     keyedarrays=false, method=ext)
    @test s.solverinfo.converged
    @test isapprox(maximum(abs.(s.nodeflux)),
                   maximum(abs.(ref.nodeflux)); rtol=1e-8)
end

end

# The augmented problem: with an explicit direct current block the system
# posed is not the harmonic one, and the whole derivative surface has to be
# the augmented system's or a caller ends up differentiating a different
# problem from the one that was solved. Every derivative below is checked
# against something independent -- the assembled Jacobian, its transpose, or
# a finite difference of the derivative one order below.
@testset verbose=true "the direct current augmented problem" begin
    circuit = Circuit(
        [:p1 => Port(1; Z0 = 50.0), :cc => Capacitor(100e-15),
         :jj => JosephsonJunction(1000e-12), :cj => Capacitor(1000e-15)],
        [[(:p1,1),(:cc,1)], [(:cc,2),(:jj,1),(:cj,1)],
         [(:p1,2),(:jj,2),(:cj,2), Ground]])
    src = [(mode = (1,), port = 1, current = 1.2e-6),
           (mode = (0,), port = 1, current = 1.0e-7)]
    wp = (2*pi*4.75e9,)
    kw = (; dc = true, odd = true, even = true)

    prob = JC.hbnonlinearproblem(wp, (4,), src, circuit, Dict{Any,Any}();
        kw...)
    @test JC.isaugmented(prob)
    L = prob.augmentation.work.layout
    @test length(prob) == JC.canonicaldim(L)
    @test L.nvdc > 0

    n = length(prob)
    Random.seed!(20)
    u = randn(n); v = randn(n); w = randn(n); z = randn(n)

    J = JC.hbjacobian!(copy(prob.jacobian), prob, u)
    Jv = zeros(n); JC.hbjvp!(Jv, prob, u, v)
    @test isapprox(J*v, Jv; rtol = 1e-12)

    Jt = zeros(n); JC.hbvjp!(Jt, prob, u, z)
    @test isapprox(transpose(Matrix(J))*z, Jt; rtol = 1e-12)

    # the operator sets its point once and must agree with the assembled one
    Jop = JC.JacobianOperator(prob, u)
    @test size(Jop) == (n, n)
    @test isapprox(Jop*v, J*v; rtol = 1e-12)
    @test isapprox(transpose(Jop)*z, Jt; rtol = 1e-12)

    # the preconditioner is the canonical wrapper, and applies
    pc = JC.preconditioner(prob, u)
    y = similar(v); ldiv!(y, pc, v)
    @test all(isfinite, y)
    @test !all(iszero, y)

    # the residual against the product, and the higher derivatives against
    # finite differences of the one below
    h = 1e-6
    F0 = zeros(n); JC.hbresidual!(F0, prob, u)
    F1 = zeros(n); JC.hbresidual!(F1, prob, u .+ h.*v)
    @test hbnlp_relerr((F1 .- F0)./h, Jv) < 1e-5

    d2 = zeros(n); JC.hbd2F!(d2, prob, u, v, w)
    a0 = zeros(n); JC.hbjvp!(a0, prob, u, v)
    a1 = zeros(n); JC.hbjvp!(a1, prob, u .+ h.*w, v)
    @test hbnlp_relerr((a1 .- a0)./h, d2) < 1e-4

    d3 = zeros(n); JC.hbd3F!(d3, prob, u, v, w, z)
    b0 = zeros(n); JC.hbd2F!(b0, prob, u, v, w)
    b1 = zeros(n); JC.hbd2F!(b1, prob, u .+ h.*z, v, w)
    @test hbnlp_relerr((b1 .- b0)./h, d3) < 1e-4

    # the drive derivative is exact, and the direct current injection moves
    # with the drive rather than staying at the scale it was built with
    dp = zeros(n); JC.setdrive!(prob, 1.0); JC.hbdFdp!(dp, prob)
    g1 = zeros(n); JC.drivenresidual!(g1, prob, u, 1.0 + h)
    g0 = zeros(n); JC.drivenresidual!(g0, prob, u, 1.0)
    @test hbnlp_relerr((g1 .- g0)./h, dp) < 1e-6
    JC.setdrive!(prob, 1.0)

    # and an external solver driving the augmented problem reaches the point
    # the internal one does
    @testset "an external solver sees the augmented system" begin
        ext = JC.ExternalSolver() do p, u0
            u = copy(u0); F = similar(u)
            JC.hbresidual!(F, p, u); r0 = norm(F)
            Jm = copy(p.jacobian)
            for k in 1:60
                norm(F) <= 1e-13*max(r0,1) && return (u, true)
                JC.hbjacobian!(Jm, p, u)
                u .-= Jm \ F
                JC.hbresidual!(F, p, u)
            end
            return (u, norm(F) <= 1e-13*max(r0,1))
        end
        a = JC.hbnlsolve(wp, (4,), src, circuit, Dict{Any,Any}();
            kw..., keyedarrays = false, method = ext)
        b = JC.hbnlsolve(wp, (4,), src, circuit, Dict{Any,Any}();
            kw..., keyedarrays = false, method = :newton, rtol = 1e-13)
        @test a.solverinfo.converged
        @test isapprox(a.dcnodevoltage, b.dcnodevoltage; rtol = 1e-8)
        @test maximum(abs, a.S .- b.S) < 1e-10
        # the node fluxes differ by whole flux quanta, which is the additive
        # static gauge and not a different state
        turns = (a.nodeflux .- b.nodeflux) ./ (2*pi)
        @test all(x -> isapprox(x, round(real(x)); atol = 1e-8), turns)
    end
end
