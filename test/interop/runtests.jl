# The interoperability test suite. Deliberately NOT part of Pkg.test:
# Krylov.jl and SciMLBase live in this directory's own environment so the
# main test suite carries no resolve or precompile cost for them. Run it
# directly:
#
#     julia test/interop/runtests.jl
#
# The script activates its own environment, develops the package into it
# on first run, and instantiates; afterwards it is a plain test run.
import Pkg
Pkg.activate(@__DIR__)
let pkgpath = normpath(joinpath(@__DIR__, "..", ".."))
    deps = Pkg.project().dependencies
    haskey(deps, "JosephsonCircuits") ||
        Pkg.develop(Pkg.PackageSpec(path = pkgpath))
    Pkg.instantiate()
end

using JosephsonCircuits
isdefined(Main, :testjpacircuit) || include(joinpath(@__DIR__, "..", "testcircuits.jl"))
using Krylov
using LinearAlgebra
using SparseArrays
using Random
using Test

const JCX = JosephsonCircuits
ext_relerr(a, b) = norm(a - b)/max(norm(a), norm(b), eps())

function ext_problem(; assemble = true)
    circuit, defs = testjpacircuitnumeric()
    return (JCX.hbnonlinearproblem((2*pi*4.75001e9,), (8,),
                [(mode=(1,),port=1,current=0.00565e-6)], circuit, defs;
                assemblejacobian = assemble),
            circuit, defs)
end

@testset verbose=true "extensions" begin

# An extension which fails to precompile is a warning, not an error, so it
# can rot silently while every test still passes. Assert it loaded before
# testing anything that depends on it.
@testset "the extensions load" begin
    using SciMLBase
    for name in (:JosephsonCircuitsKrylovExt, :JosephsonCircuitsSciMLBaseExt)
        @test !isnothing(Base.get_extension(JosephsonCircuits, name))
    end
end

@testset "SciMLBase.NonlinearProblem construction" begin
    prob, _, _ = ext_problem()
    np = SciMLBase.NonlinearProblem(prob)
    n = length(prob)
    u = 0.05 .* collect(range(-1, 1; length = n))
    F = zeros(n)
    np.f(F, u, np.p)
    @test F == JCX.hbresidual!(zeros(n), prob, u)
    Jv = zeros(n); v = collect(range(1, 2; length = n))
    np.f.jvp(Jv, v, u, np.p)
    @test Jv == JCX.hbjvp!(zeros(n), prob, u, v)
    @test np.f.jac_prototype isa SparseArrays.SparseMatrixCSC
    npf = SciMLBase.NonlinearProblem(ext_problem(assemble = false)[1])
    @test isnothing(npf.f.jac_prototype)
end

@testset "Krylov.jl drives the operator with no adapter" begin
    Random.seed!(20260828)
    prob, circuit, defs = ext_problem()
    n = length(prob)
    u = 0.05 .* randn(n)
    F = JCX.hbresidual!(zeros(n), prob, u)
    Jr = copy(prob.jacobian); JCX.hbjacobian!(Jr, prob, u)
    exact = Jr \ (-F)
    J = JCX.JacobianOperator(prob, u)
    P = JCX.preconditioner(prob, u)

    for (nm, f) in (("gmres", Krylov.gmres), ("fgmres", Krylov.fgmres),
                    ("bicgstab", Krylov.bicgstab))
        @testset "$nm" begin
            x, st = f(J, -F; rtol = 1e-12, atol = 0.0, itmax = 8n)
            @test st.solved
            @test ext_relerr(x, exact) < 1e-6
            # the preconditioner goes straight in: Krylov.jl applies `N`
            # with mul!, which AbstractPreconditioner defines
            xp, stp = f(J, -F; N = P, rtol = 1e-12, atol = 0.0, itmax = 8n)
            @test stp.solved
            @test ext_relerr(xp, exact) < 1e-6
        end
    end

    # preconditioning must actually reduce the work, or the object is
    # plumbed in without doing anything
    _, s0 = Krylov.gmres(J, -F; rtol = 1e-10, atol = 0.0, itmax = 8n)
    _, s1 = Krylov.gmres(J, -F; N = P, rtol = 1e-10, atol = 0.0, itmax = 8n)
    @test s1.niter < s0.niter

    # a solve against the adjoint exercises the matrix-free vjp
    b = randn(n)
    x, st = Krylov.gmres(J', b; rtol = 1e-12, atol = 0.0, itmax = 8n)
    @test st.solved
    @test ext_relerr(x, transpose(Matrix(Jr)) \ b) < 1e-6
end

@testset "a Newton-Krylov loop in user code" begin
    # no assembled Jacobian anywhere in this test
    prob, circuit, defs = ext_problem(; assemble = false)
    @test isnothing(prob.jacobian)
    n = length(prob)
    u = zeros(n); F = JCX.hbresidual!(zeros(n), prob, u); r0 = norm(F)
    its = 0
    for k in 1:40
        norm(F) <= 1e-13*max(r0,1) && break
        J = JCX.JacobianOperator(prob, u)
        P = JCX.preconditioner(prob, u)
        d, st = Krylov.gmres(J, -F; N = P, rtol = 1e-10, atol = 0.0, itmax = 4n)
        st.solved || break
        u .+= d; JCX.hbresidual!(F, prob, u); its += 1
    end
    @test norm(F) <= 1e-13*max(r0,1)
    @test its < 15
    ref = JCX.hbnlsolve((2*pi*4.75001e9,), (8,),
        [(mode=(1,),port=1,current=0.00565e-6)], circuit, defs;
        keyedarrays = false, ftol = 1e-14)
    @test ref.solverinfo.converged
    urefc = vec(collect(ref.nodeflux))
    uref = JCX.complex_to_real!(zeros(prob.modelayout.rdim), urefc,
        prob.modelayout.isreal)
    @test ext_relerr(u[1:length(uref)], uref) < 1e-7
end

@testset "KrylovJL as the linear solver of the package's own loop" begin
    _, circuit, defs = ext_problem()
    wp = (2*pi*4.75001e9,); src = [(mode=(1,),port=1,current=0.00565e-6)]
    ref = JCX.hbnlsolve(wp, (8,), src, circuit, defs;
                        keyedarrays = false, ftol = 1e-14)
    @test ref.solverinfo.converged
    for m in (:gmres, :fgmres, :bicgstab)
        s = JCX.hbnlsolve(wp, (8,), src, circuit, defs; keyedarrays = false,
                          ftol = 1e-14, linearsolver = JCX.KrylovJL(m))
        @test s.solverinfo.converged
        @test isapprox(maximum(abs.(s.nodeflux)),
                       maximum(abs.(ref.nodeflux)); rtol = 1e-8)
    end
    # with mode coupling and with deflation recycling, which is the path a
    # dead reference in the extension once hid because every benchmark left
    # the preconditioner nothing
    for (cm, rec) in ((:all, 0), (:none, 12))
        s = JCX.hbnlsolve(wp, (8,), src, circuit, defs; keyedarrays = false,
            ftol = 1e-14, linearsolver = JCX.KrylovJL(:gmres),
            krylovcouplingmodes = cm, krylovrecycle = rec)
        @test s.solverinfo.converged
        @test isapprox(maximum(abs.(s.nodeflux)),
                       maximum(abs.(ref.nodeflux)); rtol = 1e-8)
    end
end

end
