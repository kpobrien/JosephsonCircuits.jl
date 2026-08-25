using JosephsonCircuits
using LinearAlgebra
using SparseArrays
using Test

# The structure aware assembly reads the circuit's structure directly rather
# than a precomputed segmented gather, and is now the only assembly there is.
# Its correctness against the physics is established in test/realjacobian.jl,
# which checks it against the exact matrix-free Jacobian-vector product of
# HBSystem and against central finite differences of the residual.
#
# What is checked here is what those cannot see: that the two orientations
# describe the same matrix, that the transposed assembly really is the row
# major order of the untransposed one, and that the preconditioner restricts
# the Jacobian to the coupling it was asked for. The references are
# `sparse(transpose(.))` and `cscvaluepermutation`, neither of which shares
# code with the assembly.

# the transposed counterpart of `structurejacobian`
function structurejacobiantransposed(d, Ami, Amc, Ljb, Lmean, Rbnm, Nmodes,
    Nbranches, Nfreq, invLnm, Gnm, Cnm, rl, cl)
    P, ns = JosephsonCircuits.realjacobianstructure(Ami, Amc, Ljb, Rbnm,
        Nmodes, Nbranches, invLnm, Gnm, Cnm, rl, cl; transposed = true)
    return P, JosephsonCircuits.planstructurerealjacobian(P, eltype(P), Ami,
        Amc, Ljb, Lmean, ns, d.sys.invLnm, d.sys.Gnm, d.sys.Cnm,
        d.sys.wmodesm, d.sys.wmodes2m, rl, cl, Nmodes, Nfreq,
        JosephsonCircuits.CPU(); transposed = true)
end

@testset verbose=true "structureassembly" begin

    @variables Rleft Rright Cg Lj Cj
    circuit = Tuple{String,String,String,Num}[]
    push!(circuit,("P1_0","1","0",1))
    push!(circuit,("R1_0","1","0",Rleft))
    push!(circuit,("C1_0","1","0",Cg/2))
    push!(circuit,("Lj1_2","1","2",Lj))
    push!(circuit,("C1_2","1","2",Cj))
    for j in 2:6
        push!(circuit,("C$(j)_0","$(j)","0",Cg))
        push!(circuit,("Lj$(j)_$(j+1)","$(j)","$(j+1)",Lj))
        push!(circuit,("C$(j)_$(j+1)","$(j)","$(j+1)",Cj))
    end
    push!(circuit,("C7_0","7","0",Cg/2))
    push!(circuit,("R7_0","7","0",Rright))
    push!(circuit,("P7_0","7","0",2))
    circuitdefs = Dict(Lj => JosephsonCircuits.IctoLj(1e-6), Cg => 45e-15,
        Cj => 55e-15, Rleft => 50.0, Rright => 50.0)

    # single tone, and two tone which has self-conjugate modes and negative
    # frequencies from the multi-dimensional transform
    cases = (((2*pi*5e9,), (4,), [(mode=(1,), port=1, current=1e-6)]),
             ((2*pi*5e9, 2*pi*3e9), (2,2), [(mode=(1,0), port=1, current=1e-6)]))

    @testset "the two orientations describe the same matrix" begin
        for (wp, Nharmonics, sources) in cases
            d = JosephsonCircuits.hbnlsolve(wp, Nharmonics, sources, circuit,
                circuitdefs; debugJacobian = true)
            sys = d.sys; ml = d.modelayout; Nm = d.Nmodes
            # `:none` is the mode block diagonal, `:all` the full Jacobian,
            # and a partial set is the case whose coupling mask is lower
            # triangular and whose pattern is therefore not symmetric
            for couplingmodes in (:none, :all, [2])
                S = couplingmodes === :none ? Int[] :
                    couplingmodes === :all ? collect(1:Nm) : couplingmodes
                mask = JosephsonCircuits.modecouplingmask(Nm, S)
                Ami = JosephsonCircuits.restrictmodecoupling(
                    d.Amatrixindicesaliased, mask)
                Amc = JosephsonCircuits.restrictmodecoupling(
                    d.Amatrixconjindices, mask)
                args = (Ami, Amc, d.Ljb, d.Lmean, d.Rbnm, Nm, d.Nbranches,
                        d.Nfreq, d.invLnm, d.Gnm, d.Cnm, ml, ml)
                J, plan = JosephsonCircuits.structurejacobian(d, args...)
                Jt, plant = structurejacobiantransposed(d, args...)

                ref = sparse(transpose(J))
                @test size(Jt) == reverse(size(J))
                @test SparseArrays.getcolptr(Jt) == SparseArrays.getcolptr(ref)
                @test rowvals(Jt) == rowvals(ref)

                # the precomputed table is the incidence triple product: four
                # entries per two terminal junction, and nothing per
                # contribution
                @test length(plan.pairrow) == 4 * nnz(d.Ljb)

                for scale in (0.0, 0.01, 1.0)
                    xr = scale .* randn(ml.rdim)
                    JosephsonCircuits.setpoint!(sys, xr)
                    JosephsonCircuits._updatecosphimatrix!(sys)
                    a = zeros(nnz(J)); b = zeros(nnz(Jt))
                    JosephsonCircuits.assemblerealjacobian!(a, plan,
                        sys.phimatrix)
                    JosephsonCircuits.assemblerealjacobian!(b, plant,
                        sys.phimatrix)
                    Ja = SparseMatrixCSC(size(J)...,
                        SparseArrays.getcolptr(J), rowvals(J), a)
                    # exactly the row major order of the same matrix
                    @test b == a[JosephsonCircuits.cscvaluepermutation(Ja)]
                    @test b == nonzeros(sparse(transpose(Ja)))
                end
            end
        end
    end

    @testset "the preconditioner assembles the Jacobian it restricts" begin
        for (wp, Nharmonics, sources) in cases
            d = JosephsonCircuits.hbnlsolve(wp, Nharmonics, sources, circuit,
                circuitdefs; debugJacobian = true)
            sys = d.sys; ml = d.modelayout; Nm = d.Nmodes
            for couplingmodes in (:none, :all)
                pc = JosephsonCircuits.ModeCouplingPreconditioner(sys,
                    d.Amatrixindicesaliased, d.Amatrixconjindices, d.Ljb,
                    d.Lmean, d.Rbnm, Nm, d.Nbranches, d.Nfreq, d.invLnm,
                    d.Gnm, d.Cnm, ml; couplingmodes = couplingmodes)
                S = couplingmodes === :none ? Int[] : collect(1:Nm)
                mask = JosephsonCircuits.modecouplingmask(Nm, S)
                Ami = JosephsonCircuits.restrictmodecoupling(
                    d.Amatrixindicesaliased, mask)
                Amc = JosephsonCircuits.restrictmodecoupling(
                    d.Amatrixconjindices, mask)
                Jref, plan = JosephsonCircuits.structurejacobian(d, Ami, Amc,
                    d.Ljb, d.Lmean, d.Rbnm, Nm, d.Nbranches, d.Nfreq,
                    d.invLnm, d.Gnm, d.Cnm, ml, ml)
                @test nnz(pc.P) == nnz(Jref)
                @test size(pc.P) == size(Jref)
                for scale in (0.0, 0.01)
                    xr = scale .* randn(ml.rdim)
                    JosephsonCircuits.updatepreconditioner!(pc, xr)
                    JosephsonCircuits.setpoint!(sys, xr)
                    JosephsonCircuits._updatecosphimatrix!(sys)
                    JosephsonCircuits.assemblerealjacobian!(nonzeros(Jref),
                        plan, sys.phimatrix)
                    @test nonzeros(pc.P) == nonzeros(Jref)
                    z = similar(xr)
                    JosephsonCircuits.applypreconditioner!(z, pc, randn(ml.rdim))
                    @test all(isfinite, z)
                end
            end
        end
    end
end
