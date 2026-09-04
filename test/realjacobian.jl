using JosephsonCircuits
using LinearAlgebra
using SparseArrays
using Test


@testset verbose=true "realjacobian" begin

    @testset "plan assembled Jr matches the matrix-free Jacobian" begin

        JosephsonCircuits.@params Rleft Cc Lj Cj
        circuit = Tuple{String,String,String,Any}[]
        push!(circuit,("P1","1","0",1))
        push!(circuit,("R1","1","0",Rleft))
        push!(circuit,("C1","1","2",Cc))
        push!(circuit,("Lj1","2","0",Lj))
        push!(circuit,("C2","2","0",Cj))
        circuitdefs = Dict(
            Lj =>1000.0e-12,
            Cc => 100.0e-15,
            Cj => 1000.0e-15,
            Rleft => 50.0,
        )

        # single-tone and two-tone (the latter has self-conjugate modes and
        # negative frequencies from the multi-dimensional RDFT).
        for (wp, sources, Nharmonics) in (
            ((2*pi*4.75001*1e9,),
                [(mode=(1,),port=1,current=0.00565e-6)], (8,)),
            ((2*pi*4.65001*1e9, 2*pi*4.85001*1e9),
                [(mode=(1,0),port=1,current=0.00565e-6*1.7),
                 (mode=(0,1),port=1,current=0.00565e-6*1.7)], (4,4)),
            )

            d = JosephsonCircuits.hbnlsolve(wp, Nharmonics, sources, circuit,
                circuitdefs; debugJacobian=true)

            # the assembled Jacobian is validated against the exact
            # matrix-free Jacobian-vector product of HBSystem (which itself
            # matches finite differences of the residual, see
            # test/hbsystem.jl), for combined and Jacobian-only
            # evaluations, at random points.
            nr = length(d.xr)
            Jvr = zeros(nr)
            for trial in 1:3
                xr = 0.5*randn(nr)
                Fr = similar(d.Fr)
                d.fjreal(Fr, d.Jr, copy(xr))
                Jr1 = copy(d.Jr)
                d.fjreal(nothing, d.Jr, copy(xr))
                @test Jr1 == d.Jr
                JosephsonCircuits.setpoint!(d.sys, xr)
                for vtrial in 1:5
                    vr = randn(nr)
                    JosephsonCircuits.jacobianvectorproduct!(Jvr, d.sys, vr)
                    @test isapprox(Jvr, d.Jr*vr, atol = 1e-11,
                        norm = v->maximum(abs,v))
                end
            end
        end
    end

    @testset "real Jacobian matches finite differences" begin

        JosephsonCircuits.@params Rleft Cc Lj Cj
        circuit = Tuple{String,String,String,Any}[]
        push!(circuit,("P1","1","0",1))
        push!(circuit,("R1","1","0",Rleft))
        push!(circuit,("C1","1","2",Cc))
        push!(circuit,("Lj1","2","0",Lj))
        push!(circuit,("C2","2","0",Cj))
        circuitdefs = Dict(Lj=>1000.0e-12, Cc=>100.0e-15, Cj=>1000.0e-15,
            Rleft=>50.0)
        sources = [(mode=(1,),port=1,current=0.00565e-6)]

        d = JosephsonCircuits.hbnlsolve((2*pi*4.75001*1e9,), (8,), sources,
            circuit, circuitdefs; debugJacobian=true)

        xr0 = 0.3*randn(length(d.xr))
        d.fjreal(nothing, d.Jr, copy(xr0))
        Jfd = zeros(length(d.Fr), length(xr0))
        h = 1e-7
        Fp = similar(d.Fr)
        Fm = similar(d.Fr)
        for j in eachindex(xr0)
            xp = copy(xr0); xp[j] += h
            xm = copy(xr0); xm[j] -= h
            d.fjreal(Fp, nothing, xp)
            d.fjreal(Fm, nothing, xm)
            Jfd[:,j] = (Fp .- Fm) ./ (2*h)
        end
        @test isapprox(Matrix(d.Jr), Jfd, atol = 1e-4,
            norm = x->maximum(abs,x))
    end

    @testset "two-tone Jacobian matches finite differences" begin

        # For multi-tone problems the exact Jacobian couples modes whose
        # differences alias back onto the sampled grid, because the
        # residual is computed with cyclic Fourier transforms. The mode
        # coupling indices of fourierindices include those couplings
        # (hbmatind with alias = true), so the assembled real Jacobian is
        # the exact derivative of the residual for multi-tone problems
        # with negative frequency modes as well. (The historical
        # construction without aliasing was not exact here, with an error
        # proportional to the drive amplitude.)

        JosephsonCircuits.@params Rleft Cc Lj Cj
        circuit = Tuple{String,String,String,Any}[]
        push!(circuit,("P1","1","0",1))
        push!(circuit,("R1","1","0",Rleft))
        push!(circuit,("C1","1","2",Cc))
        push!(circuit,("Lj1","2","0",Lj))
        push!(circuit,("C2","2","0",Cj))
        circuitdefs = Dict(Lj=>1000.0e-12, Cc=>100.0e-15, Cj=>1000.0e-15,
            Rleft=>50.0)
        wp = (2*pi*4.65001*1e9, 2*pi*4.85001*1e9)
        sources = [(mode=(1,0),port=1,current=0.00565e-6*1.7),
            (mode=(0,1),port=1,current=0.00565e-6*1.7)]

        d = JosephsonCircuits.hbnlsolve(wp, (4,4), sources, circuit,
            circuitdefs; debugJacobian=true)

        xr0 = 0.4*randn(length(d.xr))
        d.fjreal(nothing, d.Jr, copy(xr0))
        h = 1e-7
        Fp = similar(d.Fr)
        Fm = similar(d.Fr)
        vr = randn(length(d.xr))
        xp = xr0 .+ h.*vr
        xm = xr0 .- h.*vr
        d.fjreal(Fp, nothing, xp)
        d.fjreal(Fm, nothing, xm)
        @test isapprox(d.Jr*vr, (Fp .- Fm)./(2*h), atol = 1e-4,
            norm = x->maximum(abs,x))
    end

    @testset "branch between internal nodes" begin

        # Josephson junctions between pairs of ungrounded nodes exercise all
        # four incidence matrix scatter targets of the direct assembly.
        JosephsonCircuits.@params Ljx Cg Cjx Rl
        Ncells = 5
        circuit = Tuple{String,String,String,Any}[]
        push!(circuit, ("P1","1","0",1))
        push!(circuit, ("R1","1","0",Rl))
        for i in 1:Ncells
            push!(circuit, ("Lj$(i)","$(i)","$(i+1)",Ljx))
            push!(circuit, ("C$(i)","$(i+1)","0",Cg))
            push!(circuit, ("Cj$(i)","$(i)","$(i+1)",Cjx))
        end
        push!(circuit, ("R2","$(Ncells+1)","0",Rl))
        circuitdefs = Dict(Ljx=>200.0e-12, Cg=>50.0e-15, Cjx=>100.0e-15,
            Rl=>50.0)
        sources = [(mode=(1,),port=1,current=1.0e-6)]

        d = JosephsonCircuits.hbnlsolve((2*pi*5.0*1e9,), (8,), sources,
            circuit, circuitdefs; debugJacobian=true)

        nr = length(d.xr)
        Jvr = zeros(nr)
        for trial in 1:3
            xr = 0.3*randn(nr)
            d.fjreal(nothing, d.Jr, copy(xr))
            JosephsonCircuits.setpoint!(d.sys, xr)
            for vtrial in 1:5
                vr = randn(nr)
                JosephsonCircuits.jacobianvectorproduct!(Jvr, d.sys, vr)
                @test isapprox(Jvr, d.Jr*vr, atol = 1e-11,
                    norm = v->maximum(abs,v))
            end
        end

        # solutions from the two methods should agree closely
        sol_qn = JosephsonCircuits.hbnlsolve((2*pi*5.0*1e9,), (8,), sources,
            circuit, circuitdefs; ftol=1e-12, method = QuasiNewton())
        sol_n = JosephsonCircuits.hbnlsolve((2*pi*5.0*1e9,), (8,), sources,
            circuit, circuitdefs; ftol=1e-12, method = Newton())
        @test isapprox(sol_qn.nodeflux[:], sol_n.nodeflux[:], atol=1e-10)
    end



end
