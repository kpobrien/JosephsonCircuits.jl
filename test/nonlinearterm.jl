using Symbolics
using JosephsonCircuits
using LinearAlgebra
using SparseArrays
using Test

@testset verbose=true "nonlinearterm" begin

    @variables Rleft Cc Lj Cj
    circuitjpa = Tuple{String,String,String,Num}[]
    push!(circuitjpa,("P1","1","0",1))
    push!(circuitjpa,("R1","1","0",Rleft))
    push!(circuitjpa,("C1","1","2",Cc))
    push!(circuitjpa,("Lj1","2","0",Lj))
    push!(circuitjpa,("C2","2","0",Cj))
    circuitdefsjpa = Dict(Lj=>1000.0e-12, Cc=>100.0e-15, Cj=>1000.0e-15,
        Rleft=>50.0)

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
    # auxiliary branch currents by the modified nodal analysis formulation,
    # so the incidence matrix gains structurally empty columns
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
            nc = length(d.x)
            mask = sys.modelayout.isreal

            vr = randn(nr)
            wr = randn(nr)
            vc = JosephsonCircuits.real_to_complex(vr, mask)
            wc = JosephsonCircuits.real_to_complex(wr, mask)

            # the real representation entry points go through the index
            # mapped kernels of the plan, the complex ones do not, so the
            # complex path is the independent reference for all three.
            JosephsonCircuits.setpoint!(sys, d.xr)
            Fc = zeros(ComplexF64, nc)
            JosephsonCircuits.residual!(Fc, sys)
            Fr = zeros(nr)
            JosephsonCircuits.setpoint!(sys, d.xr)
            JosephsonCircuits.residual!(Fr, sys)
            @test isapprox(Fr,
                JosephsonCircuits.complex_to_real(Fc, mask), rtol=1e-12)

            JosephsonCircuits.setpoint!(sys, d.xr)
            Jc = zeros(ComplexF64, nc)
            JosephsonCircuits.jacobianvectorproduct!(Jc, sys, vc)
            Jr = zeros(nr)
            JosephsonCircuits.setpoint!(sys, d.xr)
            JosephsonCircuits.jacobianvectorproduct!(Jr, sys, vr)
            @test isapprox(Jr,
                JosephsonCircuits.complex_to_real(Jc, mask), rtol=1e-12)

            JosephsonCircuits.setpoint!(sys, d.xr)
            Hc = zeros(ComplexF64, nc)
            JosephsonCircuits.hessianvectorproduct!(Hc, sys, vc, wc)
            Hr = zeros(nr)
            JosephsonCircuits.setpoint!(sys, d.xr)
            JosephsonCircuits.hessianvectorproduct!(Hr, sys, vr, wr)
            @test isapprox(Hr,
                JosephsonCircuits.complex_to_real(Hc, mask), rtol=1e-12)

            # both representations must agree with the assembled exact real
            # Jacobian, which is built by an independent code path
            # (planrealjacobian). Since the two representations now share the
            # plan, this and the central finite difference checks in
            # test/hbsystem.jl are what makes them independently verified
            # rather than merely consistent with each other.
            Jasm = copy(d.Jr)
            JosephsonCircuits.setpoint!(sys, d.xr)
            JosephsonCircuits.jacobian!(Jasm, sys)
            @test isapprox(Jr, Jasm*vr, rtol=1e-10)
            @test isapprox(JosephsonCircuits.complex_to_real(Jc, mask),
                Jasm*vr, rtol=1e-10)

            # the residual of the complex representation likewise, through
            # the assembled Jacobian at a nearby point: a central finite
            # difference of the residual along vr reproduces the product
            h = 1e-6
            Fp = zeros(nr); Fm = zeros(nr)
            JosephsonCircuits.setpoint!(sys, d.xr .+ h.*vr)
            JosephsonCircuits.residual!(Fp, sys)
            JosephsonCircuits.setpoint!(sys, d.xr .- h.*vr)
            JosephsonCircuits.residual!(Fm, sys)
            @test isapprox((Fp .- Fm)./(2h), Jr, rtol=1e-5)

            # the point is stored in both representations and they agree
            JosephsonCircuits.setpoint!(sys, d.xr)
            @test isapprox(sys.xr, d.xr, rtol=1e-12)
            @test isapprox(sys.x,
                JosephsonCircuits.real_to_complex(d.xr, mask), rtol=1e-12)

            # setting the point from either representation gives the same
            # residual
            Fr2 = zeros(nr)
            JosephsonCircuits.setpoint!(sys, sys.x)
            JosephsonCircuits.residual!(Fr2, sys)
            @test isapprox(Fr, Fr2, rtol=1e-12)

            # the collapsed linear term matrix reproduces the three separate
            # frequency dependent terms
            z = randn(ComplexF64, nc)
            @test isapprox(sys.Knm*z,
                sys.invLnm*z + im*(sys.Gnm*(sys.wmodesm*z)) -
                sys.Cnm*(sys.wmodes2m*z), rtol=1e-12)

            # # the real representation entry points do not allocate
            # JosephsonCircuits.setpoint!(sys, d.xr)
            # JosephsonCircuits.jacobianvectorproduct!(Jr, sys, vr)
            # @test (@allocated JosephsonCircuits.jacobianvectorproduct!(Jr, sys, vr)) == 0
            # JosephsonCircuits.residual!(Fr, sys)
            # @test (@allocated JosephsonCircuits.residual!(Fr, sys)) == 0
        end
    end

    @testset "plan rejects an incidence row with three entries" begin
        d = JosephsonCircuits.hbnlsolve((2*pi*4.75001*1e9,), (8,),
            [(mode=(1,),port=1,current=0.00565e-6)], circuitjpa,
            circuitdefsjpa; debugJacobian=true)
        sys = d.sys
        # a branch touching three nodes cannot arise from an incidence
        # matrix, and the flat two entry forward map must say so rather
        # than silently dropping the third contribution. build an Rbnm with
        # three entries in every branch row, keeping the mode-diagonal
        # structure branchnodesandsigns requires.
        Nmodes = length(sys.freqindexmap)
        Nbranches = size(sys.Rbnm, 1) ÷ Nmodes
        I = Int[]; J = Int[]; V = Float64[]
        for b in 1:Nbranches, m in 1:Nmodes, node in 1:3
            push!(I, (b-1)*Nmodes + m)
            push!(J, (node-1)*Nmodes + m)
            push!(V, 1.0)
        end
        # the column count must cover the three fake node blocks even when
        # the system itself has fewer columns (nothing is promoted, so
        # there are no auxiliary columns to borrow)
        Rbnm3 = sparse(I, J, V, size(sys.Rbnm,1),
            max(size(sys.Rbnm,2), 3*Nmodes))
        @test_throws ArgumentError JosephsonCircuits.plannonlinearterm(
            Rbnm3, sys.Ljb, sys.Lmean, Nbranches, sys.freqindexmap,
            sys.conjsourceindices, sys.conjtargetindices, sys.phimatrix,
            sys.Knm, sys.modelayout)
    end
    @testset "tobackend adopts a host vector instead of copying it" begin
        # on CPU() the vector is already what the plan wants, so it is taken
        # as is rather than duplicated. anything which is not a Vector is
        # still materialized into one.
        v = [1.0, 2.0, 3.0]
        @test JosephsonCircuits.tobackend(JosephsonCircuits.CPU(), v) === v
        r = JosephsonCircuits.tobackend(JosephsonCircuits.CPU(), 1.0:3.0)
        @test r isa Vector{Float64}
        @test r == [1.0, 2.0, 3.0]
    end

    @testset "the real form of the linear term is optional" begin
        d = JosephsonCircuits.hbnlsolve((2*pi*4.75001*1e9,), (8,),
            [(mode=(1,),port=1,current=0.00565e-6)], circuitjpa,
            circuitdefsjpa; debugJacobian=true)
        sys = d.sys
        Nmodes = length(sys.freqindexmap)
        args = (sys.Rbnm, sys.Ljb, sys.Lmean, size(sys.Rbnm,1) ÷ Nmodes,
            sys.freqindexmap, sys.conjsourceindices, sys.conjtargetindices,
            sys.phimatrix, sys.Knm, sys.modelayout)
        full = JosephsonCircuits.plannonlinearterm(args...)
        lean = JosephsonCircuits.plannonlinearterm(args...; realbackward=false)

        @test JosephsonCircuits.hasrealbackward(full)
        @test !JosephsonCircuits.hasrealbackward(lean)
        @test isempty(lean.kptr) && isempty(lean.kidx) && isempty(lean.kcoef)

        # the complex representation is untouched by the omission: the two
        # plans produce the same backward map, bit for bit
        nc = full.ncomplex
        xc = randn(ComplexF64, nc)
        JosephsonCircuits.applyforwardterm!(sys.phimatrix, full, xc)
        outfull = zeros(ComplexF64, nc)
        outlean = zeros(ComplexF64, nc)
        JosephsonCircuits.applybackwardterm!(outfull, full, sys.phimatrix, xc)
        JosephsonCircuits.applybackwardterm!(outlean, lean, sys.phimatrix, xc)
        @test outfull == outlean

        # and the real representation says so rather than reading the empty
        # arrays, with or without the linear term
        xr = randn(sys.modelayout.rdim)
        outr = zeros(sys.modelayout.rdim)
        @test_throws ArgumentError JosephsonCircuits.applybackwardterm!(
            outr, lean, sys.phimatrix, xr)
        @test_throws ArgumentError JosephsonCircuits.applybackwardterm!(
            outr, lean, sys.phimatrix, xr; addlinearterm = false)
    end
end
