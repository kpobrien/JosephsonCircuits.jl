using JosephsonCircuits
using LinearAlgebra
using SparseArrays
using Test

# The structure aware assembly reads the circuit's structure directly rather
# than a precomputed segmented gather. Its correctness against the physics
# is established in test/harmonics/system.jl, which checks it against the
# exact matrix-free Jacobian-vector product of HBSystem and against central
# finite differences of the residual.
#
# What is checked here is what those cannot see: that the two orientations
# describe the same matrix, that the transposed assembly really is the row
# major order of the untransposed one, and that the preconditioner restricts
# the Jacobian to the coupling it was asked for. The references are
# `sparse(transpose(.))` and `cscvaluepermutation`, neither of which shares
# code with the assembly.

# the transposed counterpart of `structurejacobian`
function structurejacobiantransposed(d, Ami, Amc, Ljb, Lscale, Rbnm, Nmodes,
    Nbranches, Nfreq, invLnm, Gnm, Cnm, layout)
    P, _ = JosephsonCircuits.realjacobianstructure(Ami, Amc, Ljb, Rbnm,
        Nmodes, Nbranches, invLnm, Gnm, Cnm, layout; transposed = true)
    junctions = JosephsonCircuits.junctionstructure(eltype(P), Ami, Amc, Ljb,
        Lscale, Rbnm, Nmodes, Nbranches, Nfreq, JosephsonCircuits.CPU())
    return P, JosephsonCircuits.planstructurerealjacobian(P, eltype(P),
        junctions, d.sys.invLnm, d.sys.Gnm, d.sys.Cnm, d.sys.wmodesm,
        d.sys.wmodes2m, layout, JosephsonCircuits.CPU(); transposed = true)
end

@testset verbose=true "structureassembly" begin

    circuit = Any[]
    push!(circuit,("P1_0", "1", "0", Port(1; Z0 = :Rleft)))
    push!(circuit,("C1_0", "1", "0", Capacitor(:Cghalf)))
    push!(circuit,("Lj1_2", "1", "2", JosephsonJunction(:Lj)))
    push!(circuit,("C1_2", "1", "2", Capacitor(:Cj)))
    for j in 2:6
        push!(circuit,("C$(j)_0", "$(j)", "0", Capacitor(:Cg)))
        push!(circuit,("Lj$(j)_$(j+1)", "$(j)", "$(j+1)", JosephsonJunction(:Lj)))
        push!(circuit,("C$(j)_$(j+1)", "$(j)", "$(j+1)", Capacitor(:Cj)))
    end
    push!(circuit,("C7_0", "7", "0", Capacitor(:Cghalf)))
    push!(circuit,("P7_0", "7", "0", Port(2; Z0 = :Rright)))
    circuit = Circuit(circuit)
    circuitdefs = Dict(:Lj => JosephsonCircuits.IctoLj(1e-6), :Cg => 45e-15,
        :Cghalf => 45e-15/2, :Cj => 55e-15, :Rleft => 50.0, :Rright => 50.0)

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
            for spec in (JosephsonCircuits.BlockDiagonal(), JosephsonCircuits.FullJacobian(), JosephsonCircuits.CoupledModes([2]))
                S = spec isa JosephsonCircuits.BlockDiagonal ? Int[] :
                    spec isa JosephsonCircuits.FullJacobian ? collect(1:Nm) :
                    spec.indices
                mask = JosephsonCircuits.modecouplingmask(Nm, S)
                Ami = JosephsonCircuits.restrictmodecoupling(
                    d.Amatrixindicesaliased, mask)
                Amc = JosephsonCircuits.restrictmodecoupling(
                    d.Amatrixconjindices, mask)
                args = (Ami, Amc, d.Ljb, d.Lscale, d.Rbnm, Nm, d.Nbranches,
                        d.Nfreq, d.invLnm, d.Gnm, d.Cnm, ml)
                J, plan = JosephsonCircuits.structurejacobian(d, args...)
                Jt, plant = structurejacobiantransposed(d, args...)

                ref = sparse(transpose(J))
                @test size(Jt) == reverse(size(J))
                @test SparseArrays.getcolptr(Jt) == SparseArrays.getcolptr(ref)
                @test rowvals(Jt) == rowvals(ref)

                # the precomputed table is the incidence triple product: four
                # entries per two terminal junction, and nothing per
                # contribution
                @test length(plan.junctions.pairrow) == 4 * nnz(d.Ljb)

                for scale in (0.0, 0.01, 1.0)
                    xr = scale .* randn(ml.rdim)
                    JosephsonCircuits.setpoint!(sys, xr)
                    cosfd = JosephsonCircuits.cosphimatrix(sys)
                    a = zeros(nnz(J)); b = zeros(nnz(Jt))
                    JosephsonCircuits.assemblerealjacobian!(a, plan, cosfd)
                    JosephsonCircuits.assemblerealjacobian!(b, plant, cosfd)
                    Ja = SparseMatrixCSC(size(J)...,
                        SparseArrays.getcolptr(J), rowvals(J), a)
                    # exactly the row major order of the same matrix
                    @test b == a[JosephsonCircuits.cscvaluepermutation(Ja)]
                    @test b == nonzeros(sparse(transpose(Ja)))
                end
            end
        end
    end

    @testset "a refreshed plan gathers the linear term at the new values" begin
        # With zero cosine coefficients the Josephson term vanishes and a
        # plan assembles its linear term alone, which is compared with the
        # real form of `K = invLnm + im*Gnm*wmodesm - Cnm*wmodes2m` formed by
        # `linearterm` and `complex_to_real`, read at every stored entry
        for (wp, Nharmonics, sources) in cases
            d = JosephsonCircuits.hbnlsolve(wp, Nharmonics, sources, circuit,
                circuitdefs; debugJacobian = true)
            sys = d.sys; ml = d.modelayout
            moved(A) = SparseMatrixCSC(size(A)..., copy(SparseArrays.getcolptr(A)),
                copy(rowvals(A)), nonzeros(A) .* (1 .+ rand(nnz(A))))
            L2, G2, C2 = moved(d.invLnm), moved(d.Gnm), moved(d.Cnm)
            K = JosephsonCircuits.linearterm(L2, G2, C2, sys.wmodesm,
                sys.wmodes2m)
            Kr = JosephsonCircuits.complex_to_real(K, ml, ml)
            zerofd = zero(JosephsonCircuits.cosphimatrix(sys))
            args = (d.Amatrixindicesaliased, d.Amatrixconjindices, d.Ljb,
                d.Lscale, d.Rbnm, d.Nmodes, d.Nbranches, d.Nfreq, d.invLnm,
                d.Gnm, d.Cnm, ml)
            J, plan = JosephsonCircuits.structurejacobian(d, args...)
            Jt, plant = structurejacobiantransposed(d, args...)
            for (P, p, R) in ((J, plan, Kr), (Jt, plant, sparse(transpose(Kr))))
                JosephsonCircuits.refreshvalues!(p, L2, G2, C2, sys.wmodesm,
                    sys.wmodes2m, d.Ljb, d.Lscale)
                a = zeros(nnz(P))
                JosephsonCircuits.assemblerealjacobian!(a, p, zerofd)
                @test a == [R[i, j] for j in axes(P, 2) for i in rowvals(P)[nzrange(P, j)]]
            end
            # the complex plan against `K` itself
            cp = sys.complexjacobianplan
            JosephsonCircuits.refreshvalues!(cp, L2, G2, C2, sys.wmodesm,
                sys.wmodes2m, d.Ljb, d.Lscale)
            c = zeros(ComplexF64, nnz(d.Jx))
            JosephsonCircuits.assemblecomplexjacobian!(c, cp, zerofd)
            @test c == [K[i, j] for j in axes(d.Jx, 2) for i in rowvals(d.Jx)[nzrange(d.Jx, j)]]
        end
    end

    @testset "the assembly kernels on the CPU backend are the host loops" begin
        # A process with one thread assembles on the host with plain loops;
        # with threads a large plan launches the same work as kernels on the
        # CPU backend, and a device runs the per entry kernel, which decodes
        # both sides of every entry. Each is launched here directly and must
        # write exactly what the loops write.
        JC = JosephsonCircuits
        wp, Nharmonics, sources = cases[2]
        d = JC.hbnlsolve(wp, Nharmonics, sources, circuit, circuitdefs;
            debugJacobian = true)
        sys = d.sys; ml = d.modelayout
        JC.setpoint!(sys, 0.3 .* randn(ml.rdim))
        cosfd = JC.cosphimatrix(sys)
        cpu = JC.CPU(); sync() = JC.KernelAbstractions.synchronize(cpu)
        args = (d.Amatrixindicesaliased, d.Amatrixconjindices, d.Ljb, d.Lscale,
            d.Rbnm, d.Nmodes, d.Nbranches, d.Nfreq, d.invLnm, d.Gnm, d.Cnm, ml)
        for (P, plan) in (JC.structurejacobian(d, args...),
                structurejacobiantransposed(d, args...))
            js = plan.junctions
            host = zeros(plan.n); k = zeros(plan.n)
            JC.assemblerealjacobian!(host, plan, cosfd)
            JC.structureassemblykernel!(cpu, 64)(k, plan.colptr, plan.rowval,
                plan.lin, cosfd, js.pairptr, js.pairrow, js.pairjunc,
                js.paircoef, js.lmolj, js.ami, js.amc, plan.linv, plan.lptr,
                js.nmodes, js.nfreq, plan.transposed; ndrange = plan.n)
            sync(); @test isequal(k, host)
            plan.assemble!(k, plan.colptr, plan.rowval, plan.lin, cosfd,
                js.pairptr, js.pairrow, js.pairjunc, js.paircoef, js.lmolj,
                js.ami, js.amc, plan.slots, js.nfreq, plan.transposed;
                ndrange = length(plan.colptr) - 1)
            sync(); @test isequal(k, host)
            g = plan.linear; x = g.inputs
            pos = zero(g.pos)
            JC.storedpositionkernel!(cpu, 64)(pos, plan.colptr, plan.rowval,
                g.row, g.col, g.part, plan.lptr, plan.transposed;
                ndrange = length(pos))
            lin = zeros(plan.n)
            JC.linearkernel!(cpu, 64)(lin, g.pos, g.row, g.col, g.part,
                x.lcolptr, x.lrowval, x.lnzval, x.gcolptr, x.growval, x.gnzval,
                x.wm, x.ccolptr, x.crowval, x.cnzval, x.wm2;
                ndrange = length(g.pos))
            sync()
            @test pos == g.pos
            @test isequal(lin, plan.lin)
        end
        cp = sys.complexjacobianplan
        jp = cp.josephson; js = jp.junctions
        host = zeros(ComplexF64, jp.n); k = zeros(ComplexF64, jp.n)
        JC.addjosephsonterm!(host, jp, cosfd)
        JC.complexjosephsonkernel!(cpu, 64)(k, jp.colptr, jp.rowval, cosfd,
            js.pairptr, js.pairrow, js.pairjunc, js.paircoef, js.lmolj, js.ami,
            js.nmodes, js.nfreq, jp.transposed; ndrange = jp.n)
        sync(); @test isequal(k, host)
        jp.assemble!(k, jp.colptr, jp.rowval, cosfd, js.pairptr, js.pairrow,
            js.pairjunc, js.paircoef, js.lmolj, js.ami, js.nmodes, js.nfreq,
            jp.transposed; ndrange = length(jp.colptr) - 1)
        sync(); @test isequal(k, host)
        g = cp.linear; x = g.inputs
        lin = zeros(ComplexF64, jp.n)
        JC.linearkernel!(cpu, 64)(lin, g.pos, g.row, g.col, g.part, x.lcolptr,
            x.lrowval, x.lnzval, x.gcolptr, x.growval, x.gnzval, x.wm,
            x.ccolptr, x.crowval, x.cnzval, x.wm2; ndrange = length(g.pos))
        sync(); @test isequal(lin, cp.lin)
    end

    @testset "the preconditioner assembles the Jacobian it restricts" begin
        for (wp, Nharmonics, sources) in cases
            d = JosephsonCircuits.hbnlsolve(wp, Nharmonics, sources, circuit,
                circuitdefs; debugJacobian = true)
            sys = d.sys; ml = d.modelayout; Nm = d.Nmodes
            for spec in (JosephsonCircuits.BlockDiagonal(), JosephsonCircuits.FullJacobian())
                pc = JosephsonCircuits.ModeCouplingPreconditioner(sys,
                    d.Amatrixindicesaliased, d.Amatrixconjindices, d.Ljb,
                    d.Lscale, d.Rbnm, Nm, d.Nbranches, d.Nfreq, d.invLnm,
                    d.Gnm, d.Cnm, ml; spec = spec)
                S = spec isa JosephsonCircuits.BlockDiagonal ? Int[] : collect(1:Nm)
                mask = JosephsonCircuits.modecouplingmask(Nm, S)
                Ami = JosephsonCircuits.restrictmodecoupling(
                    d.Amatrixindicesaliased, mask)
                Amc = JosephsonCircuits.restrictmodecoupling(
                    d.Amatrixconjindices, mask)
                Jref, plan = JosephsonCircuits.structurejacobian(d, Ami, Amc,
                    d.Ljb, d.Lscale, d.Rbnm, Nm, d.Nbranches, d.Nfreq,
                    d.invLnm, d.Gnm, d.Cnm, ml)
                @test nnz(pc.P) == nnz(Jref)
                @test size(pc.P) == size(Jref)
                for scale in (0.0, 0.01)
                    xr = scale .* randn(ml.rdim)
                    JosephsonCircuits.updatepreconditioner!(pc, xr)
                    JosephsonCircuits.setpoint!(sys, xr)
                    JosephsonCircuits.assemblerealjacobian!(nonzeros(Jref),
                        plan, JosephsonCircuits.cosphimatrix(sys))
                    @test nonzeros(pc.P) == nonzeros(Jref)
                    z = similar(xr)
                    JosephsonCircuits.applypreconditioner!(z, pc, randn(ml.rdim))
                    @test all(isfinite, z)
                end
            end
        end
    end
end

@testset "hostsparse brings a device valued matrix home" begin
    # A DeviceValuedSparseMatrix holds its structure as the transpose,
    # because compressed sparse row of the transpose is compressed sparse
    # column of the matrix, which is the layout a device direct solver wants.
    # Bringing it home is therefore a copy back plus a sparse transpose, and
    # getting that wrong would transpose the Jacobian a pump operating point
    # retains. The struct is generic over where its parts live, so this runs
    # the same arithmetic the device path does without a device.
    using SparseArrays
    for (m, n) in ((5, 5), (4, 7), (7, 4), (1, 3))
        A = sprandn(m, n, 0.4)
        # a non-symmetric structure, or a transpose would go unnoticed
        At = SparseMatrixCSC(transpose(A))
        dv = JosephsonCircuits.DeviceValuedSparseMatrix(At, nonzeros(At))
        @test size(dv) == (m, n)
        @test nnz(dv) == nnz(A)
        B = JosephsonCircuits.hostsparse(dv)
        @test B isa SparseMatrixCSC{Float64,Int}
        @test B == A
        @test B.colptr == A.colptr && B.rowval == A.rowval
    end
    # a host matrix is returned as it is
    A = sprandn(6, 6, 0.3)
    @test JosephsonCircuits.hostsparse(A) === A
    # and reading an entry of the device valued form is still refused
    At = SparseMatrixCSC(transpose(A))
    dv = JosephsonCircuits.DeviceValuedSparseMatrix(At, nonzeros(At))
    @test_throws ArgumentError dv[1, 1]
end
