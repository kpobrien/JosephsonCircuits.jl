using JosephsonCircuits
using LinearAlgebra
using SparseArrays
using Test

@testset verbose=true "hbsystem" begin

    @variables Rleft Cc Lj Cj
    circuitjpa = Tuple{String,String,String,Num}[]
    push!(circuitjpa,("P1","1","0",1))
    push!(circuitjpa,("R1","1","0",Rleft))
    push!(circuitjpa,("C1","1","2",Cc))
    push!(circuitjpa,("Lj1","2","0",Lj))
    push!(circuitjpa,("C2","2","0",Cj))
    circuitdefsjpa = Dict(Lj=>1000.0e-12, Cc=>100.0e-15, Cj=>1000.0e-15,
        Rleft=>50.0)

    # single-tone and two-tone (the latter has self-conjugate modes and
    # negative frequencies from the multi-dimensional RDFT), and a chain
    # with junctions between pairs of ungrounded nodes.
    @variables Ljx Cg Cjx Rl
    Ncells = 4
    circuitchain = Tuple{String,String,String,Num}[]
    push!(circuitchain, ("P1","1","0",1))
    push!(circuitchain, ("R1","1","0",Rl))
    for i in 1:Ncells
        push!(circuitchain, ("Lj$(i)","$(i)","$(i+1)",Ljx))
        push!(circuitchain, ("C$(i)","$(i+1)","0",Cg))
        push!(circuitchain, ("Cj$(i)","$(i)","$(i+1)",Cjx))
    end
    push!(circuitchain, ("R2","$(Ncells+1)","0",Rl))
    circuitdefschain = Dict(Ljx=>200.0e-12, Cg=>50.0e-15, Cjx=>100.0e-15,
        Rl=>50.0)

    # a circuit with mutually coupled inductors, which are promoted to
    # auxiliary branch currents by the modified nodal analysis formulation
    @variables Lla Llb Kab Rl2
    circuitmutual = Tuple{String,String,String,Num}[]
    push!(circuitmutual, ("P1","1","0",1))
    push!(circuitmutual, ("R1","1","0",Rl2))
    push!(circuitmutual, ("C1","1","2",Cc))
    push!(circuitmutual, ("Lj1","2","0",Lj))
    push!(circuitmutual, ("C2","2","0",Cj))
    push!(circuitmutual, ("L1","2","0",Lla))
    push!(circuitmutual, ("L2","3","0",Llb))
    push!(circuitmutual, ("P2","3","0",2))
    push!(circuitmutual, ("R2","3","0",Rl2))
    push!(circuitmutual, ("K1","L1","L2",Kab))
    circuitdefsmutual = Dict(Lj=>500.0e-12, Cc=>100.0e-15, Cj=>1000.0e-15,
        Lla=>300.0e-12, Llb=>300.0e-12, Rl2=>50.0, Kab=>0.99)

    testcases = (
        ("single-tone JPA", (2*pi*4.75001*1e9,),
            [(mode=(1,),port=1,current=0.00565e-6)], (8,), circuitjpa,
            circuitdefsjpa),
        ("two-tone JPA", (2*pi*4.65001*1e9, 2*pi*4.85001*1e9),
            [(mode=(1,0),port=1,current=0.00565e-6*1.7),
             (mode=(0,1),port=1,current=0.00565e-6*1.7)], (4,4), circuitjpa,
            circuitdefsjpa),
        ("two-tone JPA larger grid", (2*pi*4.65001*1e9, 2*pi*4.85001*1e9),
            [(mode=(1,0),port=1,current=0.00565e-6*1.7),
             (mode=(0,1),port=1,current=0.00565e-6*1.7)], (7,5), circuitjpa,
            circuitdefsjpa),
        ("JJ chain", (2*pi*5.0*1e9,),
            [(mode=(1,),port=1,current=1.0e-6)], (8,), circuitchain,
            circuitdefschain),
        ("mutual inductor", (2*pi*4.75001*1e9,),
            [(mode=(1,),port=1,current=1.0e-6)], (6,), circuitmutual,
            circuitdefsmutual),
        )

    for (name, wp, sources, Nharmonics, circuit, circuitdefs) in testcases
        @testset "$name" begin

            d = JosephsonCircuits.hbnlsolve(wp, Nharmonics, sources,
                circuit, circuitdefs; debugJacobian=true)
            sys = d.sys
            nr = length(d.xr)

            # the residual through the system equals the residual through
            # the solver function
            xr = 0.4*randn(nr)
            Fr1 = zeros(nr)
            Fr2 = zeros(nr)
            d.fjreal(Fr1, nothing, copy(xr))
            JosephsonCircuits.setpoint!(sys, xr)
            JosephsonCircuits.residual!(Fr2, sys)
            @test Fr1 == Fr2

            # the matrix-free Jacobian-vector product equals the assembled
            # real Jacobian applied to the same vector, for single-tone and
            # multi-tone problems alike. The mode coupling indices include
            # the couplings which alias back onto the sampled grid, matching
            # the cyclic Fourier transforms of the residual, so the
            # assembled Jacobian is the exact derivative (the matrix-free
            # product, which matches finite differences of the residual
            # below, is the ground truth).
            JosephsonCircuits.setpoint!(sys, xr)
            JosephsonCircuits.jacobian!(d.Jr, sys)
            Jvr = zeros(nr)
            for trial in 1:5
                vr = randn(nr)
                JosephsonCircuits.jacobianvectorproduct!(Jvr, sys, vr)
                @test isapprox(Jvr, d.Jr*vr, atol = 1e-11,
                    norm = v->maximum(abs,v))
            end

            # the matrix-free Jacobian-vector product equals a central
            # finite difference of the residual
            vr = randn(nr)
            JosephsonCircuits.setpoint!(sys, xr)
            JosephsonCircuits.jacobianvectorproduct!(Jvr, sys, vr)
            h = 1e-7
            Fp = zeros(nr); Fm = zeros(nr)
            JosephsonCircuits.setpoint!(sys, xr .+ h.*vr)
            JosephsonCircuits.residual!(Fp, sys)
            JosephsonCircuits.setpoint!(sys, xr .- h.*vr)
            JosephsonCircuits.residual!(Fm, sys)
            @test isapprox(Jvr, (Fp .- Fm)./(2*h), atol = 1e-4,
                norm = v->maximum(abs,v))

            # the matrix-free Hessian-vector product equals a central
            # finite difference of the Jacobian-vector product, and is
            # symmetric in the two directions
            wr = randn(nr)
            Hvw = zeros(nr)
            Hwv = zeros(nr)
            JosephsonCircuits.setpoint!(sys, xr)
            JosephsonCircuits.hessianvectorproduct!(Hvw, sys, vr, wr)
            JosephsonCircuits.hessianvectorproduct!(Hwv, sys, wr, vr)
            @test isapprox(Hvw, Hwv, atol = 1e-12,
                norm = v->maximum(abs,v))
            Jp = zeros(nr); Jm = zeros(nr)
            JosephsonCircuits.setpoint!(sys, xr .+ h.*wr)
            JosephsonCircuits.jacobianvectorproduct!(Jp, sys, vr)
            JosephsonCircuits.setpoint!(sys, xr .- h.*wr)
            JosephsonCircuits.jacobianvectorproduct!(Jm, sys, vr)
            @test isapprox(Hvw, (Jp .- Jm)./(2*h), atol = 1e-4,
                norm = v->maximum(abs,v))

            # the complex representation entry points agree with the real
            # ones through the layout conversions
            JosephsonCircuits.setpoint!(sys, xr)
            x = JosephsonCircuits.real_to_complex(xr, d.modelayout.isreal)
            v = JosephsonCircuits.real_to_complex(vr, d.modelayout.isreal)
            F = zeros(Complex{Float64}, length(d.x))
            JosephsonCircuits.residual!(F, sys)
            @test isapprox(JosephsonCircuits.complex_to_real(F,
                d.modelayout.isreal), Fr2, atol = 1e-13)
            Jv = zeros(Complex{Float64}, length(d.x))
            JosephsonCircuits.jacobianvectorproduct!(Jvr, sys, vr)
            JosephsonCircuits.jacobianvectorproduct!(Jv, sys, v)
            @test isapprox(JosephsonCircuits.complex_to_real(Jv,
                d.modelayout.isreal), Jvr, atol = 1e-12)
            w = JosephsonCircuits.real_to_complex(wr, d.modelayout.isreal)
            Hc = zeros(Complex{Float64}, length(d.x))
            JosephsonCircuits.hessianvectorproduct!(Hvw, sys, vr, wr)
            JosephsonCircuits.hessianvectorproduct!(Hc, sys, v, w)
            @test isapprox(JosephsonCircuits.complex_to_real(Hc,
                d.modelayout.isreal), Hvw, atol = 1e-12)

            # the assembled complex Jacobian through the system matches the
            # solver function path
            Jx2 = copy(d.Jx)
            d.fj(nothing, d.Jx, copy(x))
            JosephsonCircuits.setpoint!(sys, x)
            JosephsonCircuits.jacobian!(Jx2, sys)
            @test d.Jx == Jx2
        end
    end

    @testset "two-tone with dc exact Jacobian" begin

        # a two-tone problem with the self-conjugate dc mode retained: the
        # assembled real Jacobian must equal the exact matrix-free product
        # in the presence of self-conjugate modes as well.

        d = JosephsonCircuits.hbnlsolve(
            (2*pi*4.65001*1e9, 2*pi*4.85001*1e9), (4,4),
            [(mode=(1,0),port=1,current=0.00565e-6*1.7),
             (mode=(0,1),port=1,current=0.00565e-6*1.7)],
            circuitjpa, circuitdefsjpa; dc = true, even = true,
            debugJacobian = true)
        sys = d.sys
        nr = length(d.xr)
        xr = 0.3*randn(nr)
        JosephsonCircuits.setpoint!(sys, xr)
        JosephsonCircuits.jacobian!(d.Jr, sys)
        Jvr = zeros(nr)
        for trial in 1:5
            vr = randn(nr)
            JosephsonCircuits.jacobianvectorproduct!(Jvr, sys, vr)
            @test isapprox(Jvr, d.Jr*vr, atol = 1e-11,
                norm = v->maximum(abs,v))
        end
    end

    @testset "the backend kwarg reaches the solver" begin
        # CPU() is the default, so naming it explicitly must change nothing.
        # This is what checks the kwarg is threaded through both hbnlsolve
        # methods rather than silently ignored.
        args = ((2*pi*4.75001*1e9,), (8,),
                [(mode=(1,),port=1,current=0.00565e-6)], circuitjpa,
                circuitdefsjpa)
        for m in (:quasinewton, :newtonkrylov)
            a = JosephsonCircuits.hbnlsolve(args...; method = m)
            b = JosephsonCircuits.hbnlsolve(args...; method = m,
                backend = JosephsonCircuits.CPU())
            @test a.nodeflux == b.nodeflux
            @test a.solverinfo.converged && b.solverinfo.converged
        end
    end

    @testset "both source vectors are on the backend" begin
        # the complex representation residual subtracts sys.bnm from a vector
        # which lives on the backend, so bnm has to be moved there the way
        # bnmr already was. on CPU() tobackend adopts an Array, so the only
        # host-visible consequence is that a non-Array argument is
        # materialized -- which is exactly what distinguishes a moved field
        # from one stored as handed in.
        d = JosephsonCircuits.hbnlsolve((2*pi*4.75001*1e9,), (2,),
            [(mode=(1,),port=1,current=0.00565e-6)], circuitjpa,
            circuitdefsjpa; debugJacobian=true)
        mk = bnm -> JosephsonCircuits.HBSystem(d.sys.Rbnm,
            d.sys.invLnm, d.sys.Gnm, d.sys.Cnm, d.sys.wmodesm,
            d.sys.wmodes2m, bnm, d.sys.Ljb, d.sys.Ljbm, d.sys.Lmean,
            d.Nbranches, d.sys.freqindexmap, d.sys.conjsourceindices,
            d.sys.conjtargetindices, copy(d.sys.phimatrix),
            copy(d.sys.worktd), d.sys.irfftplan, d.sys.rfftplan,
            d.sys.modelayout, nothing, nothing)

        sys = mk(view(d.sys.bnm, :))
        @test sys.bnm isa Array
        @test sys.bnm == d.sys.bnm
        # the two representations of the source agree on where they live
        @test JosephsonCircuits.KernelAbstractions.get_backend(sys.bnm) ===
            JosephsonCircuits.KernelAbstractions.get_backend(sys.bnmr)

        # and the residual is unchanged by the move
        sysa = mk(copy(d.sys.bnm))
        F1 = zeros(ComplexF64, length(d.x)); F2 = similar(F1)
        JosephsonCircuits.setpoint!(sysa, copy(d.x))
        JosephsonCircuits.residual!(F1, sysa)
        JosephsonCircuits.setpoint!(d.sys, copy(d.x))
        JosephsonCircuits.residual!(F2, d.sys)
        @test F1 == F2
    end

    @testset "jacobian! without a plan throws" begin
        d = JosephsonCircuits.hbnlsolve((2*pi*4.75001*1e9,), (2,),
            [(mode=(1,),port=1,current=0.00565e-6)], circuitjpa,
            circuitdefsjpa; debugJacobian=true)
        sys = JosephsonCircuits.HBSystem(d.sys.Rbnm,
            d.sys.invLnm, d.sys.Gnm, d.sys.Cnm, d.sys.wmodesm,
            d.sys.wmodes2m, d.sys.bnm, d.sys.Ljb, d.sys.Ljbm, d.sys.Lmean,
            d.Nbranches, d.sys.freqindexmap, d.sys.conjsourceindices,
            d.sys.conjtargetindices, copy(d.sys.phimatrix),
            copy(d.sys.worktd), d.sys.irfftplan, d.sys.rfftplan,
            d.sys.modelayout, nothing, nothing)
        JosephsonCircuits.setpoint!(sys, copy(d.xr))
        @test_throws ArgumentError JosephsonCircuits.jacobian!(d.Jr, sys)
        @test_throws ArgumentError JosephsonCircuits.jacobian!(d.Jx, sys)
    end


    @testset "applyffttranspose! bilinear transpose property" begin
        # the defining property, <alpha, td> == <P, applyfft(td)> under the
        # non conjugating bilinear pairing, on randomized data over one,
        # two and three tone grids (odd first dimension, as plan_applynl
        # builds), including a singleton trailing tone dimension. this pins
        # the transpose against the forward transform it claims to
        # transpose, for every grid shape the reverse contraction can see.
        for Nw in ((5,), (5,4), (4,3,3), (5,1))
            NLj = 2
            fd = zeros(Complex{Float64}, Nw..., NLj)
            td, irfftplan, rfftplan = JosephsonCircuits.plan_applynl(fd)
            fftplan = JosephsonCircuits.plan_applyffttranspose(fd, td)
            padded = zeros(Complex{Float64}, size(td))
            alpha = zeros(Complex{Float64}, size(td))
            for trial in 1:3
                td .= randn.()
                JosephsonCircuits.applyfft!(fd, td, rfftplan)
                P = randn(Complex{Float64}, size(fd))
                JosephsonCircuits.applyffttranspose!(alpha, P, padded,
                    fftplan)
                lhs = sum(P .* fd)
                rhs = sum(alpha .* td)
                @test isapprox(lhs, rhs, rtol = 1e-12)
            end
        end
    end
end
