using JosephsonCircuits
using LinearAlgebra
using SparseArrays
using Test

@testset verbose=true "batchedblocks" begin

    @testset "layout of a synthetic uniform batch" begin
        # three 2x2 blocks, interleaved so block index runs fastest, which is
        # how the real representation stores modes
        nb, bs = 3, 2
        blockof = [mod(i-1, nb)+1 for i in 1:nb*bs]
        A = spzeros(nb*bs, nb*bs)
        for b in 1:nb
            idx = findall(==(b), blockof)
            A[idx, idx] = [10b 1.0; 2.0 20b]
        end
        L = JosephsonCircuits.blockdiagonallayout(A, blockof, nb)
        @test L isa JosephsonCircuits.BatchedBlockLayout
        @test L.nblocks == nb
        @test L.blocksize == bs
        # each extracted block matches ordinary submatrix extraction
        for b in 1:nb
            idx = findall(==(b), blockof)
            @test Matrix(JosephsonCircuits.blockmatrix(A, L, b)) ==
                Matrix(A[idx, idx])
        end
        # gather / scatter round trips
        v = collect(1.0:nb*bs)
        rhs = zeros(bs, nb)
        JosephsonCircuits.gatherblocks!(rhs, v, L)
        out = zeros(nb*bs)
        JosephsonCircuits.scatterblocks!(out, rhs, L)
        @test out == v
        vals = zeros(size(L.nzindex))
        JosephsonCircuits.gatherblocks!(vals, A, L)
        @test vals[:, 1] == nonzeros(JosephsonCircuits.blockmatrix(A, L, 1))
        @test_throws DimensionMismatch JosephsonCircuits.gatherblocks!(
            zeros(1, 1), A, L)
    end

    @testset "layout refuses what it cannot batch" begin
        nb, bs = 2, 2
        blockof = [1, 2, 1, 2]
        # not block diagonal: an entry couples block 1 to block 2
        A = spzeros(4, 4)
        A[1,1] = 1.0; A[3,3] = 1.0; A[2,2] = 1.0; A[4,4] = 1.0
        A[1,2] = 1.0
        @test isnothing(JosephsonCircuits.blockdiagonallayout(A, blockof, nb))
        # blocks of unequal size, as a dc (self-conjugate) mode produces
        @test isnothing(JosephsonCircuits.blockdiagonallayout(
            sparse(I, 3, 3), [1, 1, 2], 2))
        # blocks with different sparsity patterns
        B = spzeros(4, 4)
        B[1,1] = 1.0; B[3,3] = 1.0; B[3,1] = 1.0     # block 1 has an extra entry
        B[2,2] = 1.0; B[4,4] = 1.0
        @test isnothing(JosephsonCircuits.blockdiagonallayout(B, blockof, nb))
        @test isnothing(JosephsonCircuits.blockdiagonallayout(
            sparse(I, 4, 4), blockof, 0))
    end

    @testset "the real preconditioner is a uniform batch" begin
        circuit = Tuple{String,String,String,Union{Complex{Float64},Symbol,Int64}}[]
        push!(circuit,("P1","1","0",1))
        push!(circuit,("R1","1","0",:Rleft))
        push!(circuit,("C1","1","2",:Cc))
        push!(circuit,("Lj1","2","0",:Lj))
        push!(circuit,("C2","2","0",:Cj))
        defs = Dict{Symbol,Complex{Float64}}(:Lj=>1000e-12, :Cc=>100.0e-15,
            :Cj=>1000e-15, :Rleft=>50.0)
        wp = 2*pi*5e9
        sources = [(mode=(1,),port=1,current=1.0e-6)]
        d = JosephsonCircuits.hbnlsolve((wp,), (8,), sources, circuit, defs;
            debugJacobian = true)
        Nmodes = d.Nmodes
        blockof = JosephsonCircuits.modeslotindex(d.modelayout)

        keep = JosephsonCircuits.modecouplingmask(Nmodes, Int[])
        P, plan = JosephsonCircuits.structurejacobian(d,
            JosephsonCircuits.restrictmodecoupling(d.Amatrixindicesaliased, keep),
            JosephsonCircuits.restrictmodecoupling(d.Amatrixconjindices, keep),
            d.Ljb, d.Lmean, d.Rbnm, Nmodes, d.Nbranches, d.Nfreq,
            d.invLnm, d.Gnm, d.Cnm, d.modelayout, d.modelayout)
        x = 0.3*randn(length(d.xr))
        JosephsonCircuits.setpoint!(d.sys, x)
        JosephsonCircuits.jacobian!(P, plan, d.sys)

        L = JosephsonCircuits.blockdiagonallayout(P, blockof, Nmodes)
        @test L isa JosephsonCircuits.BatchedBlockLayout
        @test L.nblocks == Nmodes
        @test L.blocksize*Nmodes == size(P, 1)

        # the batch reproduces the blocks, and solving blockwise reproduces a
        # solve against the whole preconditioner
        for b in 1:Nmodes
            idx = findall(==(b), blockof)
            @test Matrix(JosephsonCircuits.blockmatrix(P, L, b)) ≈ Matrix(P[idx, idx])
        end
        v = randn(size(P, 1))
        rhs = zeros(L.blocksize, Nmodes)
        JosephsonCircuits.gatherblocks!(rhs, v, L)
        sol = similar(rhs)
        for b in 1:Nmodes
            sol[:, b] = JosephsonCircuits.blockmatrix(P, L, b) \ rhs[:, b]
        end
        z = zeros(size(P, 1))
        JosephsonCircuits.scatterblocks!(z, sol, L)
        @test P*z ≈ v rtol=1e-8

        # the pattern does not change with the operating point, so one layout
        # is valid for the whole Newton solve
        JosephsonCircuits.setpoint!(d.sys, 0.37 .* x)
        JosephsonCircuits.jacobian!(P, plan, d.sys)
        L2 = JosephsonCircuits.blockdiagonallayout(P, blockof, Nmodes)
        @test L2.nzindex == L.nzindex
        @test L2.rowval == L.rowval
        @test L2.colptr == L.colptr

        # a preconditioner with retained coupling is not block diagonal, so
        # the batched path correctly declines it
        keep2 = JosephsonCircuits.modecouplingmask(Nmodes, [1, 2])
        P2, plan2 = JosephsonCircuits.structurejacobian(d,
            JosephsonCircuits.restrictmodecoupling(d.Amatrixindicesaliased, keep2),
            JosephsonCircuits.restrictmodecoupling(d.Amatrixconjindices, keep2),
            d.Ljb, d.Lmean, d.Rbnm, Nmodes, d.Nbranches, d.Nfreq,
            d.invLnm, d.Gnm, d.Cnm, d.modelayout, d.modelayout)
        JosephsonCircuits.jacobian!(P2, plan2, d.sys)
        @test isnothing(JosephsonCircuits.blockdiagonallayout(P2, blockof, Nmodes))
    end

    @testset "a dc mode is correctly refused" begin
        # self-conjugate modes collapse to one real slot, so the blocks are
        # not all the same size and the batched path must decline
        dccircuit = Tuple{String,String,String,Union{Complex{Float64},Symbol,Int64}}[]
        push!(dccircuit,("P1","1","0",1))
        push!(dccircuit,("R1","1","0",:Rleft))
        push!(dccircuit,("L1","1","0",:Lm))
        push!(dccircuit,("K1","L1","L2",:K1))
        push!(dccircuit,("C1","1","2",:Cc))
        push!(dccircuit,("L2","2","3",:Lm))
        push!(dccircuit,("Lj3","3","0",:Lj))
        push!(dccircuit,("Lj4","2","0",:Lj))
        push!(dccircuit,("C2","2","0",:Cj))
        dcdefs = Dict{Symbol,Complex{Float64}}(:Lj=>2000e-12, :Lm=>10e-12,
            :Cc=>200.0e-15, :Cj=>900e-15, :Rleft=>50.0, :Rright=>50.0, :K1=>0.9)
        dcsources = [(mode=(0,),port=1,current=50e-5),
                     (mode=(1,),port=1,current=0.0001e-6)]
        d = JosephsonCircuits.hbnlsolve((2*pi*5e9,), (2,), dcsources,
            dccircuit, dcdefs; dc = true, debugJacobian = true)
        @test any(d.modelayout.isreal)
        Nmodes = d.Nmodes
        blockof = JosephsonCircuits.modeslotindex(d.modelayout)
        keep = JosephsonCircuits.modecouplingmask(Nmodes, Int[])
        P, plan = JosephsonCircuits.structurejacobian(d,
            JosephsonCircuits.restrictmodecoupling(d.Amatrixindicesaliased, keep),
            JosephsonCircuits.restrictmodecoupling(d.Amatrixconjindices, keep),
            d.Ljb, d.Lmean, d.Rbnm, Nmodes, d.Nbranches, d.Nfreq,
            d.invLnm, d.Gnm, d.Cnm, d.modelayout, d.modelayout)
        JosephsonCircuits.setpoint!(d.sys, 0.1*randn(length(d.xr)))
        JosephsonCircuits.jacobian!(P, plan, d.sys)
        @test isnothing(JosephsonCircuits.blockdiagonallayout(P, blockof, Nmodes))
    end

    @testset "csc to csr value permutation" begin
        # a device sparse matrix is CSR and a host one CSC, so their stored
        # values differ by a permutation. CSR of A is CSC of transpose(A),
        # which is what pins the permutation without needing a device.
        for (m, n, d) in ((40, 40, 0.1), (30, 50, 0.2), (7, 7, 0.5), (5, 5, 1.0))
            A = JosephsonCircuits.SparseArrays.sprandn(m, n, d)
            p = JosephsonCircuits.cscvaluepermutation(A)
            @test length(p) == nnz(A)
            @test sort(p) == collect(1:nnz(A))          # it is a permutation
            @test nonzeros(A)[p] ==
                nonzeros(JosephsonCircuits.SparseArrays.sparse(transpose(A)))
        end
        # an empty pattern is still a valid, empty permutation
        @test isempty(JosephsonCircuits.cscvaluepermutation(spzeros(4, 4)))
    end

    @testset "factorization routing" begin
        A = sparse([2.0 0.0; 0.0 3.0])
        gpu = JosephsonCircuits.CUDSSFactorization(batched = true)
        cpu = JosephsonCircuits.KLUfactorization()

        # no GPU factorization supplied: the ordinary CPU path, which is the
        # default everywhere and not a degraded mode
        @test JosephsonCircuits.batchedfactorization(A, [1,2], 2, nothing, cpu) === cpu

        # blocks that are not a uniform batch fall back to the unbatched GPU
        # path rather than failing
        B = sparse(Matrix{Float64}(I, 3, 3))
        @test JosephsonCircuits.batchedfactorization(B, [1,1,2], 2, gpu, cpu) === gpu

        # a genuine uniform batch selects the batched path, which without the
        # extension loaded reports what is missing rather than a MethodError
        @test_throws ArgumentError JosephsonCircuits.batchedfactorization(
            A, [1,2], 2, gpu, cpu)

        # a factorization with no batched form is used unchanged even when the
        # matrix is a uniform batch, which is what makes the batch opt-in: the
        # default cuDSS factorization measured slower batched than whole
        plain = JosephsonCircuits.CUDSSFactorization()
        @test isnothing(plain.batched)
        @test JosephsonCircuits.batchedfactorization(A, [1,2], 2, plain, cpu) === plain
        @test JosephsonCircuits.batchedfactorization(A, [1,2], 2, cpu, cpu) === cpu

        # and the unloaded GPU factorization is constructible and errors clearly
        @test_throws ArgumentError gpu.factorize(A)
        @test_throws ArgumentError JosephsonCircuits.cudssbatchedfactorization(nothing)
    end

    @testset "the gather and scatter run on the layout's backend" begin
        # a block diagonal with three uniform blocks of two
        A = JosephsonCircuits.SparseArrays.sparse(
            [1,2,1,2,3,4,3,4,5,6,5,6], [1,1,2,2,3,3,4,4,5,5,6,6],
            Float64.(1:12))
        blockof = [1,1,2,2,3,3]
        layout = JosephsonCircuits.blockdiagonallayout(A, blockof, 3)
        @test layout !== nothing

        # moving the layout to CPU() leaves it usable and only touches slots
        moved = JosephsonCircuits.tobackend(JosephsonCircuits.CPU(), layout)
        @test moved.slots == layout.slots
        @test moved.nzindex === layout.nzindex
        @test moved.colptr === layout.colptr
        @test moved.rowval === layout.rowval

        # the gather and the scatter are inverse permutations, through the
        # kernel path, and agree with the plain indexing they replaced
        v = randn(6)
        rhs = zeros(layout.blocksize, layout.nblocks)
        JosephsonCircuits.gatherblocks!(rhs, v, moved)
        for b in 1:layout.nblocks, i in 1:layout.blocksize
            @test rhs[i, b] == v[layout.slots[i, b]]
        end
        w = fill(NaN, 6)
        JosephsonCircuits.scatterblocks!(w, rhs, moved)
        @test w == v

        # every slot is written exactly once, which is what makes the scatter
        # safe to run as one work item per slot with no atomics
        @test sort(vec(layout.slots)) == 1:6

        # the shared pattern, and the row major index map a device solver
        # wants: the values it selects are the row major values of each block
        pat = JosephsonCircuits.blockpattern(layout)
        @test size(pat) == (layout.blocksize, layout.blocksize)
        @test nnz(pat) == size(layout.nzindex, 1)
        rowmajor = JosephsonCircuits.blockrowmajorindex(layout)
        @test size(rowmajor) == size(layout.nzindex)
        @test sort(vec(rowmajor)) == sort(vec(layout.nzindex))
        for b in 1:layout.nblocks
            blk = JosephsonCircuits.blockmatrix(A, layout, b)
            @test nonzeros(A)[rowmajor[:, b]] ==
                nonzeros(JosephsonCircuits.SparseArrays.sparse(transpose(blk)))
        end

        # the generic gather, which the batched path uses for both the right
        # hand sides and the values
        got = zeros(size(rowmajor))
        JosephsonCircuits.gathervalues!(got, nonzeros(A), rowmajor)
        @test got == nonzeros(A)[rowmajor]
        @test_throws DimensionMismatch JosephsonCircuits.gathervalues!(
            zeros(1, 1), nonzeros(A), rowmajor)

        JosephsonCircuits.gatherblocks!(rhs, v, moved)
        @test_throws DimensionMismatch JosephsonCircuits.gatherblocks!(
            zeros(layout.blocksize + 1, layout.nblocks), v, moved)
        @test_throws DimensionMismatch JosephsonCircuits.scatterblocks!(
            w, zeros(layout.blocksize + 1, layout.nblocks), moved)
    end

    @testset "tobackend adopts a host matrix instead of copying it" begin
        m = [1 2; 3 4]
        @test JosephsonCircuits.tobackend(JosephsonCircuits.CPU(), m) === m
        t = JosephsonCircuits.tobackend(JosephsonCircuits.CPU(), view(m, :, :))
        @test t isa Matrix{Int} && t == m
    end

end
