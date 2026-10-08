using JosephsonCircuits
using LinearAlgebra
using SparseArrays
using Random
using Test

isdefined(Main, :testjpacircuit) || include(joinpath(@__DIR__, "..", "testcircuits.jl"))

@testset verbose=true "nonlinearterm" begin

    # the nonlinear term plan's maps against the complex representation,
    # the assembled Jacobian and finite differences are checked per circuit
    # class in test/harmonics/system.jl
    circuitjpa, circuitdefsjpa = testjpacircuit()

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
            Rbnm3, sys.Ljb, sys.Lscale, Nbranches, sys.freqindexmap,
            sys.conjsourceindices, sys.conjtargetindices, sys.phimatrix,
            sys.Knm, sys.modelayout)
    end
    @testset "the plan's kernels on the CPU backend are its host loops" begin
        # A process with one thread runs every map of the plan as a plain
        # loop; with threads, a large plan launches the same work items as
        # kernels on the CPU backend. Those are launched here directly and
        # must write exactly what the loops write.
        JC = JosephsonCircuits
        d = JC.hbnlsolve((2*pi*4.75001*1e9,), (8,),
            [(mode=(1,),port=1,current=0.00565e-6)], circuitjpa,
            circuitdefsjpa; debugJacobian=true)
        sys = d.sys; p = sys.nonlineartermplan
        sync() = JC.KernelAbstractions.synchronize(p.backend)
        nr, nc = sys.modelayout.rdim, p.ncomplex
        xr = randn(nr); xc = randn(ComplexF64, nc)
        fd = similar(sys.phimatrix); fk = similar(fd)

        JC.applyforwardterm!(fd, p, xr)
        p.forward!(fk, p.n1, p.s1, p.n2, p.s2, p.flags, xr; ndrange = p.nslots)
        sync(); @test isequal(fk, fd)
        JC.applyforwardterm!(fd, p, xc)
        p.forwardcomplex!(fk, p.cn1, p.s1, p.cn2, p.s2, p.flags, xc;
            ndrange = p.nslots)
        sync(); @test isequal(fk, fd)

        o, ok = zeros(nr), zeros(nr)
        JC.applybackwardterm!(o, p, fd, xr)
        p.backward!(ok, p.bptr, p.bsrc, p.bcoef, fd, p.kptr, p.kidx, p.kcoef,
            xr, p.lptr, p.lwide; ndrange = nc)
        sync(); @test isequal(ok, o)
        oc, okc = zeros(ComplexF64, nc), zeros(ComplexF64, nc)
        JC.applybackwardterm!(oc, p, fd, xc)
        p.backwardcomplex!(okc, p.bptr, p.bsrc, p.bcoef, fd, p.cptr, p.cidx,
            p.ccoef, xc; ndrange = nc)
        sync(); @test isequal(okc, oc)

        JC.applyrealtocomplex!(oc, p, xr)
        p.realtocomplex!(okc, xr, p.lptr, p.lwide; ndrange = nc)
        sync(); @test isequal(okc, oc)
        JC.applycomplextoreal!(o, p, xc)
        p.complextoreal!(ok, xc, p.lptr, p.lwide; ndrange = nc)
        sync(); @test isequal(ok, o)

        # the transposed maps of the problem interface
        tp = JC.plannonlineartermtranspose(p, sys.modelayout, sys.phimatrix,
            sys.phitd)
        JC.applybackwardjosephsontranspose!(fd, tp, p, xr)
        tp.backwardtranspose!(fk, tp.tbptr, tp.tbnode, tp.tbcoef, xr, p.lptr,
            p.lwide, tp.gtscale; ndrange = tp.nslots)
        sync(); @test isequal(fk, fd)
        JC.applyforwardtranspose!(o, tp, fd, xr)
        tp.forwardtranspose!(ok, tp.tfptr, tp.tfslot, tp.tfcoef, tp.tfimag, fd,
            tp.ktptr, tp.ktrow, tp.ktcoef, xr; ndrange = tp.rdim)
        sync(); @test isequal(ok, o)

        # and the refreshes of the values
        maps = JC.valuemaps(sys)
        JC.refreshvalues!(p, maps, sys.Knm, sys.Ljb, sys.Lscale)
        k, c, b = zero(p.kcoef), zero(p.ccoef), zero(p.bcoef)
        JC.refreshrealkernel!(p.backend, 64)(k, maps.knz, maps.kmap, maps.kre,
            maps.kim; ndrange = length(k))
        JC.refreshcomplexkernel!(p.backend, 64)(c, maps.knz, maps.cmap;
            ndrange = length(c))
        JC.refreshjosephsonkernel!(p.backend, 64)(b,
            JC.junctioncoefficients(Float64, sys.Ljb, sys.Lscale), maps.bjunc,
            maps.bsgn; ndrange = length(b))
        sync()
        @test isequal(k, p.kcoef) && isequal(c, p.ccoef) && isequal(b, p.bcoef)
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
        args = (sys.Rbnm, sys.Ljb, sys.Lscale, size(sys.Rbnm,1) ÷ Nmodes,
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

    @testset "a rebound single precision plan is one built at its values" begin
        # a system rebound to new component values holds what a system built
        # at them holds, in single precision too, where `Lscale/Lj` rounded
        # once and the quotient of the rounded `Lscale` and `Lj` differ for
        # these inductances
        JC = JosephsonCircuits
        chain(Ljs) = Circuit(vcat(Any[(:p1, 1, 0, Port(1))],
            [(Symbol(:jj, i), i, i + 1, JosephsonJunction(Ljs[i]))
                for i in eachindex(Ljs)],
            [(Symbol(:c, i), i + 1, 0, Capacitor(40e-15))
                for i in eachindex(Ljs)]))
        system(Ljs) = JC.hbnlsolve((2*pi*5e9,), (2,),
            [(mode = (1,), port = 1, current = 1e-7)], chain(Ljs);
            returnsystem = true, assemblejacobian = false,
            method = NewtonKrylov(precision = Float32)).sys
        s1 = system([100e-12, 110e-12])
        s2 = system([120e-12, 142.9e-12])
        r = JC.rebind!(s1, s2.invLnm, s2.Gnm, s2.Cnm,
            Vector{ComplexF64}(s2.bnm), s2.Ljb, s2.Ljbm, s2.Lscale;
            maps = JC.valuemaps(s1), relations = nothing)
        p, q = r.nonlineartermplan, s2.nonlineartermplan
        @test p.bcoef == q.bcoef && p.kcoef == q.kcoef && p.ccoef == q.ccoef
        x = Float32.(0.1 .* sin.(1:length(r.xr)))
        F1, F2 = similar(x), similar(x)
        JC.residual!(F1, JC.setpoint!(r, x))
        JC.residual!(F2, JC.setpoint!(s2, x))
        @test F1 == F2
    end
end
