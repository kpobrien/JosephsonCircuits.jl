using JosephsonCircuits
using LinearAlgebra
using SparseArrays
using Test

# Test-local oracle: independent construction of the Josephson branch matrix
# from the mode coupling index matrix and the Fourier coefficients of
# cos(phi(t)). Entries follow Amatrixindices (negative entries denote complex
# conjugation, zeros denote dropped couplings and produce no stored entry),
# scaled by Lmean/Lj per branch.
function testAoLjbm(phimatrix, Amatrixindices, Ljb, Lmean, Nmodes, Nbranches)
    Nfreq = prod(size(phimatrix)[1:end-1])
    I = Int[]
    J = Int[]
    V = Complex{Float64}[]
    for (bi, b) in enumerate(Ljb.nzind)
        offset = (bi-1)*Nfreq
        for j in 1:Nmodes, i in 1:Nmodes
            ind = Amatrixindices[i, j]
            if ind != 0
                c = ind > 0 ? phimatrix[ind + offset] :
                    conj(phimatrix[-ind + offset])
                push!(I, (b-1)*Nmodes + i)
                push!(J, (b-1)*Nmodes + j)
                push!(V, (Lmean/Ljb.nzval[bi])*c)
            end
        end
    end
    return sparse(I, J, V, Nbranches*Nmodes, Nbranches*Nmodes)
end

@testset verbose=true "complexjacobian" begin

    @testset "plan assembled Jx matches the holomorphic derivative" begin

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
            n = length(d.x)

            for trial in 1:3
                xr = 0.5*randn(length(d.xr))
                x = JosephsonCircuits.real_to_complex(xr,
                    d.modelayout.isreal)

                # combined and Jacobian-only evaluations assemble the same
                # matrix
                F = similar(d.F)
                d.fj(F, d.Jx, copy(x))
                Jx1 = copy(d.Jx)
                d.fj(nothing, d.Jx, copy(x))
                @test Jx1 == d.Jx

                # for a single tone there are no mode couplings which alias
                # back onto the sampled grid, so the holomorphic Jacobian
                # equals the exact holomorphic (Wirtinger) derivative of
                # the residual, extracted column by column from the
                # matrix-free Jacobian-vector product: with jvp(v) =
                # A*v + B*conj(v), A_j = (jvp(e_j) - im*jvp(im*e_j))/2.
                # for multiple tones the holomorphic Jacobian deliberately
                # keeps the truncated (non-aliased) couplings (it is an
                # approximation either way; see fourierindices), so its
                # exactness is a single-tone property; the exact multi-tone
                # real Jacobian is validated against the matrix-free
                # product and finite differences in test/realjacobian.jl
                # and test/hbsystem.jl.
                if length(wp) == 1 && trial == 1
                    JosephsonCircuits.setpoint!(d.sys, xr)
                    A = zeros(Complex{Float64}, n, n)
                    Jv1 = zeros(Complex{Float64}, n)
                    Jv2 = zeros(Complex{Float64}, n)
                    ev = zeros(Complex{Float64}, n)
                    for j in 1:n
                        fill!(ev, 0)
                        ev[j] = 1
                        JosephsonCircuits.jacobianvectorproduct!(Jv1, d.sys,
                            ev)
                        ev[j] = im
                        JosephsonCircuits.jacobianvectorproduct!(Jv2, d.sys,
                            ev)
                        A[:, j] = (Jv1 .- im.*Jv2)./2
                    end
                    @test isapprox(Matrix(d.Jx), A, atol = 1e-11)
                end
            end
        end
    end

    @testset "hblinsolve signal grid matches reference" begin

        # The linearized harmonic balance solver builds its system matrix on
        # the signal frequency grid from the pump solution with the same
        # plan machinery used for the Jacobians of the nonlinear system,
        # augmented with the modified nodal analysis (MNA) equations of the
        # promoted resistors. Check the sparsity structure, the Josephson
        # (pump modulation) contribution, and the full per-frequency
        # assembly against the reference construction from the branch
        # matrices, the incidence matrix products, and the MNA stamps.
        @variables Rleft Cc Lj Cj
        circuit = Tuple{String,String,String,Num}[]
        push!(circuit,("P1","1","0",1))
        push!(circuit,("R1","1","0",Rleft))
        push!(circuit,("C1","1","2",Cc))
        push!(circuit,("Lj1","2","0",Lj))
        push!(circuit,("C2","2","0",Cj))
        circuitdefs = Dict(Lj=>1000.0e-12, Cc=>100.0e-15, Cj=>1000.0e-15,
            Rleft=>50.0)
        wp = (2*pi*4.75001*1e9,)
        sources = [(mode=(1,),port=1,current=0.00565e-6)]

        nonlinear = JosephsonCircuits.hbnlsolve(wp, (8,), sources, circuit,
            circuitdefs)

        # replicate the hblinsolve setup on the signal grid
        psc = JosephsonCircuits.parsesortcircuit(circuit)
        cg = JosephsonCircuits.calccircuitgraph(psc)
        signalfreq = JosephsonCircuits.truncfreqs(
            JosephsonCircuits.calcfreqsdft((6,)); dc=true, odd=false,
            even=true, maxintermodorder=Inf)
        Nsignalmodes = length(signalfreq.modes)
        signalnm = JosephsonCircuits.numericmatrices(psc, cg, circuitdefs,
            Nmodes = Nsignalmodes)

        pumpfreq = nonlinear.frequencies
        allpumpfreq = JosephsonCircuits.calcfreqsrdft(pumpfreq.Nharmonics)
        pumpindices = JosephsonCircuits.fourierindices(pumpfreq)
        Amatrixmodes, Amatrixindices = JosephsonCircuits.hbmatind(
            allpumpfreq, signalfreq)

        Nwtuple = (pumpfreq.Nw..., length(nonlinear.Ljb.nzval))
        phimatrix = zeros(Complex{Float64}, Nwtuple)
        phimatrixtd, irfftplan, rfftplan =
            JosephsonCircuits.plan_applynl(phimatrix)
        branchflux = nonlinear.Rbnm*nonlinear.nodeflux[:]
        JosephsonCircuits.phivectortomatrix!(
            branchflux[nonlinear.Ljbm.nzind], phimatrix,
            pumpindices.vectomatmap, pumpindices.conjsourceindices,
            pumpindices.conjtargetindices, length(nonlinear.Ljb.nzval))
        JosephsonCircuits.applynl!(phimatrix, phimatrixtd, cos, irfftplan,
            rfftplan)

        # the modified nodal analysis augmentation of the promoted resistor
        mnaindices = JosephsonCircuits.mnaportresistorindices(
            psc.componenttypes, psc.nodeindices,
            psc.mutualinductorbranchnames, signalnm.vvn)
        @test length(mnaindices) == 1
        Nauxmna = length(mnaindices)*Nsignalmodes
        Amna0, AmnaG = JosephsonCircuits.calcAmnasplit(mnaindices,
            psc.nodeindices, signalnm.vvn, Nsignalmodes, psc.Nnodes)
        Gnmp = JosephsonCircuits.calcGn(psc.componenttypes[mnaindices],
            psc.nodeindices[:, mnaindices], signalnm.vvn[mnaindices],
            Nsignalmodes, psc.Nnodes)
        Gnmsub = JosephsonCircuits.mnapad(
            JosephsonCircuits.mnasubtractpromoted(signalnm.Gnm, Gnmp),
            Nauxmna)
        invLnmp = JosephsonCircuits.mnapad(signalnm.invLnm, Nauxmna)
        Cnmp = JosephsonCircuits.mnapad(signalnm.Cnm, Nauxmna)
        Rbnmmna = hcat(signalnm.Rbnm, spzeros(eltype(signalnm.Rbnm),
            size(signalnm.Rbnm,1), Nauxmna))

        # the test-local oracle construction
        AoLjbm = testAoLjbm(phimatrix, Amatrixindices,
            signalnm.Ljb, 1, Nsignalmodes, cg.Nbranches)
        Rbnmt = sparse(transpose(Rbnmmna))
        AoLjnm = Rbnmt*AoLjbm*Rbnmmna
        wpumpmodes = JosephsonCircuits.calcmodefreqs(nonlinear.w,
            signalfreq.modes)
        wmodes1 = 2*pi*5.0e9 .+ wpumpmodes
        Cnmcopy = JosephsonCircuits.freqsubst(Cnmp, wmodes1, nothing)
        Gnmcopy = JosephsonCircuits.freqsubst(Gnmsub, wmodes1, nothing)
        invLnmcopy = JosephsonCircuits.freqsubst(invLnmp, wmodes1, nothing)
        Asparseref = JosephsonCircuits.spaddkeepzeros(
            JosephsonCircuits.spaddkeepzeros(
                JosephsonCircuits.spaddkeepzeros(
                    JosephsonCircuits.spaddkeepzeros(
                        JosephsonCircuits.spaddkeepzeros(AoLjnm,
                            invLnmcopy), Gnmcopy), Cnmcopy), Amna0), AmnaG)

        # the plan construction, as used by hblinsolve
        lsys = JosephsonCircuits.HBLinearizedSystem(Amatrixindices,
            signalnm.Ljb, Rbnmmna, Nsignalmodes, cg.Nbranches, phimatrix,
            invLnmcopy, Gnmcopy, Cnmcopy, invLnmp, Gnmsub, Cnmp,
            Int[], Int[], Int[], Amna0, AmnaG, nothing, wpumpmodes,
            psc.Nnodes)
        Asparse = lsys.Asparse

        # identical sparsity structure, including stored zeros
        @test SparseArrays.getcolptr(Asparse) ==
            SparseArrays.getcolptr(Asparseref)
        @test rowvals(Asparse) == rowvals(Asparseref)

        # the Josephson contribution and its conjugate match the reference
        AoLjnmnzval = zeros(Complex{Float64}, nnz(Asparse))
        JosephsonCircuits.addjosephsonterm!(AoLjnmnzval,
            lsys.complexjacobianplan, phimatrix)
        A1 = copy(Asparse)
        copyto!(A1.nzval, AoLjnmnzval)
        @test isapprox(Matrix(A1), Matrix(AoLjnm), atol = 1e-14)

        fill!(AoLjnmnzval, 0)
        JosephsonCircuits.addjosephsonterm!(AoLjnmnzval,
            lsys.complexjacobianplan, phimatrix, true)
        copyto!(A1.nzval, AoLjnmnzval)
        @test isapprox(Matrix(A1), conj.(Matrix(AoLjnm)), atol = 1e-14)

        # the full per-frequency assembly through HBLinearizedSystem
        # matches the reference construction, including the negative
        # frequency mode conjugations and the MNA augmentation, for both
        # the pump modulation matrix and its complex conjugate (the
        # adjoint system).
        for ws in (2*pi*5.0e9, 2*pi*0.3e9)
            wmodes = ws .+ wpumpmodes
            wouter = size(Asparse,1) ÷ length(wmodes)
            wmodesm = JosephsonCircuits.LinearAlgebra.Diagonal(
                repeat(wmodes, outer = wouter))
            wmodes2m = JosephsonCircuits.LinearAlgebra.Diagonal(
                repeat(wmodes.^2, outer = wouter))
            for conjugatepump in (false, true)
                # the reference per-frequency construction
                Aref = copy(Asparse)
                fill!(Aref.nzval, 0)
                AoLjnmuse = conjugatepump ? conj.(AoLjnm) : AoLjnm
                JosephsonCircuits.sparseadd!(Aref, 1, AoLjnmuse,
                    JosephsonCircuits.sparseaddmap(Aref, AoLjnmuse))
                JosephsonCircuits.sparseaddconjsubst!(Aref, -1, Cnmp,
                    wmodes2m, JosephsonCircuits.sparseaddmap(Aref, Cnmp),
                    real.(wmodesm) .< 0, wmodesm, Int[], nothing)
                JosephsonCircuits.sparseaddconjsubst!(Aref, im, Gnmsub,
                    wmodesm, JosephsonCircuits.sparseaddmap(Aref, Gnmsub),
                    real.(wmodesm) .< 0, wmodesm, Int[], nothing)
                JosephsonCircuits.sparseaddconjsubst!(Aref, 1, invLnmp,
                    JosephsonCircuits.LinearAlgebra.Diagonal(
                        ones(size(invLnmp, 1))),
                    JosephsonCircuits.sparseaddmap(Aref, invLnmp),
                    real.(wmodesm) .< 0, wmodesm, Int[], nothing)
                JosephsonCircuits.sparseadd!(Aref, 1, Amna0,
                    JosephsonCircuits.sparseaddmap(Aref, Amna0))
                JosephsonCircuits.sparseaddconjsubst!(Aref, im, AmnaG,
                    wmodesm, JosephsonCircuits.sparseaddmap(Aref, AmnaG),
                    real.(wmodesm) .< 0, wmodesm, Int[], nothing)

                # the shared assembly, both entry points. the reference
                # above uses the Diagonal based sparseaddconjsubst! method,
                # so this also checks the mode indexed method the assembly
                # uses against it.
                A2 = copy(Asparse)
                JosephsonCircuits.assemblesystemmatrix!(A2, lsys, wmodes;
                    conjugatepump = conjugatepump)
                @test isapprox(A2, Aref, atol = 1e-14)
                A3 = copy(Asparse)
                JosephsonCircuits.assemblesystemmatrix!(A3, lsys, ws;
                    conjugatepump = conjugatepump)
                @test A2 == A3
            end
        end
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

        # single tone, so the holomorphic Jacobian equals the exact
        # holomorphic (Wirtinger) derivative extracted from the matrix-free
        # Jacobian-vector product (see the first testset).
        n = length(d.x)
        for trial in 1:3
            x = 0.3*randn(Complex{Float64}, n)
            d.fj(nothing, d.Jx, copy(x))
            xr = JosephsonCircuits.complex_to_real(x, d.modelayout.isreal)
            JosephsonCircuits.setpoint!(d.sys, xr)
            A = zeros(Complex{Float64}, n, n)
            Jv1 = zeros(Complex{Float64}, n)
            Jv2 = zeros(Complex{Float64}, n)
            ev = zeros(Complex{Float64}, n)
            for j in 1:n
                fill!(ev, 0)
                ev[j] = 1
                JosephsonCircuits.jacobianvectorproduct!(Jv1, d.sys, ev)
                ev[j] = im
                JosephsonCircuits.jacobianvectorproduct!(Jv2, d.sys, ev)
                A[:, j] = (Jv1 .- im.*Jv2)./2
            end
            @test isapprox(Matrix(d.Jx), A, atol = 1e-11)
        end
    end

end
