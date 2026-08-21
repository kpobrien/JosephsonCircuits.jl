using JosephsonCircuits
using LinearAlgebra
using SparseArrays
using Test

# Reference construction of the exact real Jacobian at the complex nodeflux
# x, by assembling the complex Jacobians Jx (holomorphic) and Jxconj
# (non-holomorphic) with the Josephson branch matrices and the incidence
# matrix triple products, then converting them to real form with
# complex_to_real_sum, zeroing the non-holomorphic contribution to the
# self-conjugate (eg. DC) columns. This reproduces the construction
# hbnlsolve used for method=:newton before the real Jacobian was assembled
# directly with a precomputed plan, except that the difference couplings use
# the aliased mode coupling indices (hbmatind with alias = true), which the
# plan also uses, so the assembled Jacobian is the exact derivative of the
# residual for multi-tone problems as well. Kept as an independent check.
function realjacobianreference(d, x)
    d.cosphimatrix(x)
    Am = copy(d.phimatrix)
    AoLjbm = JosephsonCircuits.calcAoLjbm2(Am, d.Amatrixindicesaliased, d.Ljb,
        d.Lmean, d.Nmodes, d.Nbranches)
    AoLjbmconj = JosephsonCircuits.calcAoLjbm2(Am, d.Amatrixconjindices,
        d.Ljb, d.Lmean, d.Nmodes, d.Nbranches)
    Rbnmt = sparse(transpose(d.Rbnm))
    Jx = Matrix(Rbnmt*AoLjbm*d.Rbnm .+ d.invLnm .+ im*d.Gnm*d.wmodesm .-
        d.Cnm*d.wmodes2m)
    Jxconj = Matrix(Rbnmt*AoLjbmconj*d.Rbnm)
    return JosephsonCircuits.complex_to_real_sum(Jx, Jxconj,
        d.modelayout.isreal, d.modelayout.isreal; realcolscale_b=0)
end

# Reference sparsity structure of the real Jacobian (including stored
# numerical zeros), from a random, everywhere nonzero, coefficient matrix,
# reproducing the construction from the complex Jacobian sparsity structures.
function realjacobianreferencepattern(d)
    Am = rand(Complex{Float64}, size(d.phimatrix))
    AoLjbm = JosephsonCircuits.calcAoLjbm2(Am, d.Amatrixindicesaliased, d.Ljb,
        d.Lmean, d.Nmodes, d.Nbranches)
    AoLjbmconj = JosephsonCircuits.calcAoLjbm2(Am, d.Amatrixconjindices,
        d.Ljb, d.Lmean, d.Nmodes, d.Nbranches)
    Rbnmt = sparse(transpose(d.Rbnm))
    Jx = JosephsonCircuits.spaddkeepzeros(
        JosephsonCircuits.spaddkeepzeros(
            JosephsonCircuits.spaddkeepzeros(Rbnmt*AoLjbm*d.Rbnm,
                d.invLnm), im*d.Gnm), -d.Cnm)
    Jxconj = Rbnmt*AoLjbmconj*d.Rbnm
    return JosephsonCircuits.complex_to_real_sum(Jx, Jxconj, d.modelayout,
        d.modelayout)
end

@testset verbose=true "realjacobian" begin

    @testset "plan assembled Jr equals reference" begin

        @variables Rleft Cc Lj Cj
        circuit = Tuple{String,String,String,Num}[]
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

            # the sparsity structure (including stored zeros) should match
            # the reference construction so the KLU symbolic factorization
            # behaves the same.
            Jrpat = realjacobianreferencepattern(d)
            @test SparseArrays.getcolptr(d.Jr) == SparseArrays.getcolptr(Jrpat)
            @test rowvals(d.Jr) == rowvals(Jrpat)

            # the assembled Jacobian should agree with the reference at
            # random points, for combined and Jacobian-only evaluations
            for trial in 1:3
                xr = 0.5*randn(length(d.xr))
                x = JosephsonCircuits.real_to_complex(xr,
                    d.modelayout.isreal)
                Fr = similar(d.Fr)
                d.fjreal(Fr, d.Jr, copy(xr))
                Jrref = realjacobianreference(d, x)
                @test isapprox(Matrix(d.Jr), Jrref, atol = 1e-12,
                    rtol = 1e-12)
            end
            xr = 0.5*randn(length(d.xr))
            x = JosephsonCircuits.real_to_complex(xr, d.modelayout.isreal)
            d.fjreal(nothing, d.Jr, copy(xr))
            @test isapprox(Matrix(d.Jr), realjacobianreference(d, x),
                atol = 1e-12, rtol = 1e-12)
        end
    end

    @testset "real Jacobian matches finite differences" begin

        @variables Rleft Cc Lj Cj
        circuit = Tuple{String,String,String,Num}[]
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

        @variables Rleft Cc Lj Cj
        circuit = Tuple{String,String,String,Num}[]
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
        @variables Ljx Cg Cjx Rl
        Ncells = 5
        circuit = Tuple{String,String,String,Num}[]
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

        Jrpat = realjacobianreferencepattern(d)
        @test SparseArrays.getcolptr(d.Jr) == SparseArrays.getcolptr(Jrpat)
        @test rowvals(d.Jr) == rowvals(Jrpat)

        for trial in 1:3
            xr = 0.3*randn(length(d.xr))
            x = JosephsonCircuits.real_to_complex(xr, d.modelayout.isreal)
            d.fjreal(nothing, d.Jr, copy(xr))
            @test isapprox(Matrix(d.Jr), realjacobianreference(d, x),
                atol = 1e-12, rtol = 1e-12)
        end

        # solutions from the two methods should agree closely
        sol_qn = JosephsonCircuits.hbnlsolve((2*pi*5.0*1e9,), (8,), sources,
            circuit, circuitdefs; ftol=1e-12, method=:quasinewton)
        sol_n = JosephsonCircuits.hbnlsolve((2*pi*5.0*1e9,), (8,), sources,
            circuit, circuitdefs; ftol=1e-12, method=:newton)
        @test isapprox(sol_qn.nodeflux[:], sol_n.nodeflux[:], atol=1e-10)
    end

    @testset "realsparseadd! matches complex_to_real" begin

        # random complex sparse matrix, mixed real/complex modes
        Nmodes = 3
        Nnodes = 4
        dim = Nmodes*Nnodes
        isrealmask = [true, false, false]
        rl = JosephsonCircuits.ModeLayout(isrealmask, dim)

        M = sprand(ComplexF64, dim, dim, 0.4)
        # real form of the holomorphic map x -> M*x
        Mr = JosephsonCircuits.complex_to_real(M, isrealmask)

        # destination with the pattern of Mr
        Jr = copy(Mr)
        fill!(nonzeros(Jr), 0)
        indexmap = JosephsonCircuits.realsparseaddmap(Jr, M, rl, rl)
        JosephsonCircuits.realsparseadd!(Jr, 1, M, indexmap)
        @test isapprox(Jr, Mr, atol=1e-15)

        # with a scalar factor and a diagonal
        D = JosephsonCircuits.LinearAlgebra.Diagonal(randn(dim))
        Mr2 = JosephsonCircuits.complex_to_real(sparse((im*2.0)*M*D), isrealmask)
        fill!(nonzeros(Jr), 0)
        JosephsonCircuits.realsparseadd!(Jr, im*2.0, M, D, indexmap)
        # patterns can differ (M*D may drop entries) so compare densely
        @test isapprox(Matrix(Jr), Matrix(Mr2), atol=1e-13)
    end

end
