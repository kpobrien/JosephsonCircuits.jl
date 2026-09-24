using JosephsonCircuits
using LinearAlgebra
using SparseArrays
using Test

# CPU() through the package, as the other device path tests do, so the test
# environment needs no KernelAbstractions dependency of its own
const CPU = JosephsonCircuits.CPU

# The machinery hblinsolve uses to solve a frequency sweep on a backend. The
# kernels run on CPU() unchanged, so the assembly and the gather are tested
# here against the host assembler; the device solve itself needs cuDSS and is
# exercised where a device is available.

# A host stand-in for the batched cuDSS solver, dispatched on host arrays,
# which no production path hands it: it keeps the options it was made with
# and solves each system of the batch with the host factorization, so the
# sweep can be followed through `solvebatch!` on CPU() and its solutions
# compared with the host's.
struct HostSweep
    n::Int
    colptr::Vector{Int}
    rowval::Vector{Int}
    nzval::Matrix{ComplexF64}
    X::Array{ComplexF64,3}
    B::Array{ComplexF64,3}
    options::NamedTuple
end
function JosephsonCircuits._cudss_sweep(rowptr::Vector{<:Integer},
    colind::Vector{<:Integer}, nzval::Matrix{ComplexF64},
    X::Array{ComplexF64,3}, B::Array{ComplexF64,3}; kwargs...)
    return HostSweep(length(rowptr) - 1, Vector{Int}(rowptr),
        Vector{Int}(colind), nzval, X, B, NamedTuple(kwargs))
end
function JosephsonCircuits._cudss_sweepsolve!(S::HostSweep)
    for k in axes(S.nzval, 2)
        # the structure is compressed sparse row, so read as compressed
        # sparse column it is the transpose
        At = SparseMatrixCSC(S.n, S.n, S.colptr, S.rowval, S.nzval[:, k])
        S.X[:, :, k] .= sparse(transpose(At)) \ S.B[:, :, k]
    end
    return S
end

@testset verbose=true "the device sweep" begin

    # a chain with a Josephson junction per cell, two ports, and a resistor at
    # each end, so the linearized system has a pump modulation contribution,
    # a modified nodal analysis augmentation from the promoted port resistors,
    # and every linear term matrix populated.
    buildcasecache = Dict{Any,Any}()
    buildcase(Nmod, wp, Npump; dft = false) =
        get!(buildcasecache, (Nmod, wp, Npump, dft)) do
            buildcase_(Nmod, wp, Npump, dft)
        end
    function buildcase_(Nmod, wp, Npump, dft)
        circuit = Any[]
        push!(circuit,("P1_0", "1", "0", Port(1; Z0 = :Rleft)))
        push!(circuit,("C1_0", "1", "0", Capacitor(:Cghalf)))
        push!(circuit,("Lj1_2", "1", "2", JosephsonJunction(:Lj)))
        push!(circuit,("C1_2", "1", "2", Capacitor(:Cj)))
        for j in 2:8
            push!(circuit,("C$(j)_0", "$(j)", "0", Capacitor(:Cg)))
            push!(circuit,("Lj$(j)_$(j+1)", "$(j)", "$(j+1)", JosephsonJunction(:Lj)))
            push!(circuit,("C$(j)_$(j+1)", "$(j)", "$(j+1)", Capacitor(:Cj)))
        end
        push!(circuit,("C9_0", "9", "0", Capacitor(:Cghalf)))
        push!(circuit,("P9_0", "9", "0", Port(2; Z0 = :Rright)))
        circuit = Circuit(circuit)
        circuitdefs = Dict(:Lj => JosephsonCircuits.IctoLj(1e-6), :Cg => 45e-15,
            :Cghalf => 45e-15/2, :Cj => 55e-15, :Rleft => 50.0,
            :Rright => 50.0)
        nl = JosephsonCircuits.hbnlsolve(wp, Npump,
            [(mode=ntuple(i->i==1 ? 1 : 0, length(wp)), port=1, current=1e-6)],
            circuit, circuitdefs; keyedarrays=false)
        psc = JosephsonCircuits.compile(circuit)
        # the full grid keeps the modes below the signal, whose frequencies
        # are negative where the sidebands pass zero
        sf = dft ? JosephsonCircuits.truncfreqs(
            JosephsonCircuits.calcfreqsdft(Nmod); dc=true, odd=true,
            even=false, maxintermodorder=Inf) :
            JosephsonCircuits.removeconjfreqs(JosephsonCircuits.truncfreqs(
            JosephsonCircuits.calcfreqsrdft(Nmod); dc=true, odd=true,
            even=false, maxintermodorder=Inf))
        return psc, circuitdefs, sf, nl
    end

    @testset "the sweep assembly matches the host assembler" begin

        # the modes above the signal, and the full grid, whose modes below
        # the signal have negative frequencies, which exercise the
        # conjugation of the stored values
        wp = (2*pi*5e9,)
        for dft in (false, true)
            psc, circuitdefs, sf, nl = buildcase((4,), wp, (6,); dft = dft)
            # off the integer GHz grid, where a signal frequency plus a
            # mode frequency would land on zero
            ws = 2*pi*[0.43e9, 1.37e9, 3.11e9, 6.61e9, 9.29e9, 11.83e9]
            d = JosephsonCircuits.hblinsolve(ws, psc, circuitdefs, sf;
                nonlinear=nl, debuglsys=true)
            lsys = d.lsys
            @test JosephsonCircuits.cansweepondevice(lsys)
            # the full grid must actually reach negative mode frequencies,
            # or the conjugation of the stored values, which is the only
            # branch in the assembly kernel, goes untested
            if dft
                @test any(w -> any(<(0), w .+ d.wpumpmodes), ws)
            end

            A = copy(lsys.Asparse)
            perm = JosephsonCircuits.cscvaluepermutation(A)
            host = Matrix{ComplexF64}(undef, nnz(A), length(ws))
            for (i,w) in enumerate(ws)
                JosephsonCircuits.assemblesystemmatrix!(A, lsys,
                    w .+ d.wpumpmodes)
                host[:,i] .= nonzeros(A)[perm]
            end

            plan, rowptr, colind = JosephsonCircuits.planfrequencysweep(lsys,
                CPU())
            got = Matrix{ComplexF64}(undef, nnz(A), length(ws))
            JosephsonCircuits.assemblesweep!(got, plan, ws)
            # the per-entry quadratic reproduces the host assembly exactly:
            # the terms it folds together are added, and addition of the
            # stored values is what the host assembler does too
            @test got == host

            # the transposed structure the values are ordered for
            At = sparse(transpose(A))
            @test Array(rowptr) == SparseArrays.getcolptr(At)
            @test Array(colind) == rowvals(At)
        end
    end

    @testset "the sweep assembly rejects a symbolic frequency" begin
        JosephsonCircuits.@params Rleft Cc Lj Cj
        w = FrequencyDependent(identity)
        circuit = Any[]
        push!(circuit,("P1", "1", "0", Port(1; Z0 = Rleft)))
        push!(circuit,("C1", "1", "2", Capacitor(Cc*(1+1e-18*w))))
        push!(circuit,("Lj1", "2", "0", JosephsonJunction(Lj)))
        push!(circuit,("C2", "2", "0", Capacitor(Cj)))
        circuit = Circuit(circuit)
        circuitdefs = Dict(Lj=>1000.0e-12, Cc=>100.0e-15, Cj=>1000.0e-15,
            Rleft=>50.0)
        psc = JosephsonCircuits.compile(circuit)
        sf = JosephsonCircuits.removeconjfreqs(JosephsonCircuits.truncfreqs(
            JosephsonCircuits.calcfreqsrdft((2,)); dc=true, odd=true,
            even=false, maxintermodorder=Inf))
        d = JosephsonCircuits.hblinsolve(2*pi*[4.0e9, 5.0e9], psc,
            circuitdefs, sf; debuglsys=true)
        # a component value which depends on the frequency is not a constant
        # coefficient, so there is nothing to precompute and hblinsolve keeps
        # such a sweep on the host
        @test !JosephsonCircuits.cansweepondevice(d.lsys)
        @test_throws ArgumentError JosephsonCircuits.planfrequencysweep(d.lsys,
            CPU())
    end

    @testset "the factorization's options reach the batched solver" begin
        JC = JosephsonCircuits
        # every factorization answers, so a caller which does not know
        # which one it has can forward unconditionally
        @test JC.solverkwargs(nothing) == NamedTuple()
        @test JC.solverkwargs(JC.KLUfactorization()) == NamedTuple()
        @test JC.solverkwargs(JC.CUDSSFactorization()) == NamedTuple()
        @test JC.solverkwargs(JC.CUDSSFactorization(pivot_epsilon = 0.0,
            ir_n_steps = 0)) == (pivot_epsilon = 0.0, ir_n_steps = 0)

        # and the sweep holds them for the call which makes each
        # direction's solver, which is what lets a caller turn the sweep's
        # own defaults off
        wp = (2*pi*5e9,)
        psc, circuitdefs, sf, nl = buildcase((4,), wp, (6,))
        ws = 2*pi*[3.11e9, 6.61e9]
        d = JC.hblinsolve(ws, psc, circuitdefs, sf; nonlinear = nl,
            debuglsys = true)
        spec = (full = true, rows = Int[])
        f = JC.CUDSSFactorization(pivot_epsilon = 0.0, ir_n_steps = 0)
        @test JC.devicesolutions(d.lsys, d.bnm, ws, CPU(), spec;
            factorization = f).solverkwargs ==
            (pivot_epsilon = 0.0, ir_n_steps = 0)
        # both directions are made from the one record
        @test JC.devicesolutions(d.lsys, d.bnm, ws, CPU(), spec, spec;
            factorization = f).solverkwargs ==
            (pivot_epsilon = 0.0, ir_n_steps = 0)
        @test JC.devicesolutions(d.lsys, d.bnm, ws, CPU(),
            spec).solverkwargs == NamedTuple()
    end

    @testset "the sweep followed through solvebatch! on the host" begin
        JC = JosephsonCircuits
        # a circuit with a scattering block, whose auxiliary port current
        # rows the batch carries, and more frequencies than one batch
        # holds, so the last batch is short
        Z0 = 50.0
        circuit = Circuit([
            (:p1, 1, 0, Port(1; Z0 = Z0)),
            (:b, 1, 2, ScatteringParameters(ComplexF64[0 1; 1 0]; zref = Z0,
                noise = Lossless())),
            (:cc, 2, 3, Capacitor(100e-15)),
            (:jj, 3, 0, JosephsonJunction(1000e-12)),
            (:cj, 3, 0, Capacitor(1000e-15))])
        wp = (2*pi*5e9,)
        src = [(mode = (1,), port = 1, current = 0.8e-6)]
        nl = hbnlsolve(wp, (6,), src, circuit; keyedarrays = false)
        psc = JC.compile(circuit)
        sf = JC.removeconjfreqs(JC.truncfreqs(JC.calcfreqsrdft((4,));
            dc = true, odd = true, even = false, maxintermodorder = Inf))
        ws = 2*pi*collect(range(4.0e9, 6.0e9,
            length = JC.uniformbatchlimit(1) + 2))
        d = JC.hblinsolve(ws, psc, Dict{Any,Any}(), sf; nonlinear = nl,
            debuglsys = true)
        lsys = d.lsys
        b = Matrix{ComplexF64}(d.bnm)
        spec = (full = true, rows = Int[])
        f = JC.CUDSSFactorization(pivot_epsilon = 1e-6, ir_n_steps = 3)
        ds = JC.devicesolutions(lsys, b, ws, CPU(), spec, spec;
            factorization = f)
        @test isnothing(ds.blocks)
        @test ds.nb < length(ws)
        # each direction's solution of each frequency is the host's solution
        # of the equations assembled at that frequency, the column scaling
        # of the batch undone
        A = copy(lsys.Asparse)
        X = similar(b)
        for lo in 1:ds.nb:length(ws)
            JC.solvebatch!(ds, lo)
            for i in lo:min(lo + ds.nb - 1, length(ws))
                JC.assemblesystemmatrix!(A, lsys, ws[i] .+ d.wpumpmodes)
                JC.forwardsolution!(X, ds, i)
                @test X ≈ A \ b rtol = 1e-10
                JC.adjointsolution!(X, ds, i)
                @test X ≈ sparse(transpose(A)) \ b rtol = 1e-10
            end
        end
        # the options were handed to the call which made each direction's
        # solver
        for slot in 1:2
            @test ds.sweeps[slot].options ==
                (pivot_epsilon = 1e-6, ir_n_steps = 3)
        end
    end

    @testset "the gathered rows are the ones the ports read" begin
        wp = (2*pi*5e9,)
        psc, circuitdefs, sf, nl = buildcase((4,), wp, (6,))
        ws = 2*pi*collect(range(2.13e9, 9.41e9, length=5))
        d = JosephsonCircuits.hblinsolve(ws, psc, circuitdefs, sf;
            nonlinear=nl, debuglsys=true)
        Nmodes = d.Nmodes
        rows = JosephsonCircuits.portsolutionrows(d.nodeindices,
            d.portindices, Nmodes)
        @test issorted(rows)
        @test allunique(rows)
        @test all(r -> 1 <= r <= size(d.lsys.Asparse,1), rows)

        # the scattering parameters read a solution only through these rows:
        # zeroing every other row leaves them unchanged
        n = size(d.lsys.Asparse,1)
        nrhs = size(d.bnm,2)
        phin = randn(ComplexF64, n, nrhs)
        masked = zeros(ComplexF64, n, nrhs)
        masked[rows,:] .= phin[rows,:]
        wmodes = ws[1] .+ d.wpumpmodes
        drives = Diagonal(JosephsonCircuits.portsourcecurrents(d.bnm,
            d.portindices, d.nodeindices, Nmodes))
        ow = zeros(ComplexF64, nrhs, nrhs); ow2 = similar(ow)
        for (o, p) in ((ow, phin), (ow2, masked))
            JosephsonCircuits.calcoutputwaves!(o, p, drives, d.portindices,
                d.portimpedances, d.nodeindices, d.componenttypes, wmodes)
        end
        @test ow == ow2

        # and the gather picks exactly those rows out of a batch
        nb = 3
        X = randn(ComplexF64, n, nrhs, nb)
        out = zeros(ComplexF64, length(rows), nrhs, nb)
        JosephsonCircuits.gatherportrows!(out, X, rows, CPU())
        @test out == X[rows,:,:]
    end

    @testset "the transposed sweep plan describes the transposed system" begin
        wp = (2*pi*5e9,)
        psc, circuitdefs, sf, nl = buildcase((4,), wp, (6,))
        ws = 2*pi*[0.43e9, 2.17e9, 7.91e9]
        d = JosephsonCircuits.hblinsolve(ws, psc, circuitdefs, sf;
            nonlinear=nl, debuglsys=true)
        lsys = d.lsys
        A = copy(lsys.Asparse)

        n = size(A, 1)
        # the solver reads (rowptr, colind, values) as compressed sparse row,
        # which is the transpose of what those same arrays mean as compressed
        # sparse column
        csr(ptr, ind, val) = sparse(transpose(SparseMatrixCSC(n, n,
            Array(ptr), Array(ind), val)))

        fplan, frowptr, fcolind = JosephsonCircuits.planfrequencysweep(lsys,
            CPU())
        aplan, arowptr, acolind = JosephsonCircuits.planfrequencysweep(lsys,
            CPU(); adjoint = true)
        fgot = Matrix{ComplexF64}(undef, nnz(A), length(ws))
        agot = similar(fgot)
        JosephsonCircuits.assemblesweep!(fgot, fplan, ws)
        JosephsonCircuits.assemblesweep!(agot, aplan, ws)

        for (i,w) in enumerate(ws)
            JosephsonCircuits.assemblesystemmatrix!(A, lsys, w .+ d.wpumpmodes)
            # the forward plan hands the solver the system matrix, the adjoint
            # plan its transpose, whose solutions are the adjoint ones
            @test csr(frowptr, fcolind, fgot[:,i]) == A
            @test csr(arowptr, acolind, agot[:,i]) == sparse(transpose(A))
        end
    end

    @testset "csc to csr value permutation" begin
        # a device sparse matrix is CSR and a host one CSC, so their stored
        # values differ by a permutation. CSR of A is CSC of transpose(A),
        # which is what pins the permutation without needing a device.
        for (m, n, d) in ((40, 40, 0.1), (30, 50, 0.2), (7, 7, 0.5), (5, 5, 1.0))
            A = sprandn(m, n, d)
            p = JosephsonCircuits.cscvaluepermutation(A)
            @test length(p) == nnz(A)
            @test sort(p) == collect(1:nnz(A))          # it is a permutation
            @test nonzeros(A)[p] == nonzeros(sparse(transpose(A)))
        end
        # an empty pattern is still a valid, empty permutation
        @test isempty(JosephsonCircuits.cscvaluepermutation(spzeros(4, 4)))
    end

    @testset "tobackend adopts a host matrix instead of copying it" begin
        m = [1 2; 3 4]
        @test JosephsonCircuits.tobackend(JosephsonCircuits.CPU(), m) === m
        t = JosephsonCircuits.tobackend(JosephsonCircuits.CPU(), view(m, :, :))
        @test t isa Matrix{Int} && t == m
    end

    @testset "the noise outputs formed on the backend are the host's" begin
        # devicenoise forms the noise scattering matrix, its reductions and
        # the added noise covariance where the adjoint solutions are. Run on
        # CPU() against the host sweep, from the same adjoint solutions: the
        # noise ports of a resistor and the channels of two instances of a
        # lossy non-reciprocal block, which share their definition's
        # factors, warm, in a pumped circuit whose modes reach negative
        # frequencies, which sign the columns
        JC = JosephsonCircuits
        wp = (2*pi*4.75001e9,)
        src = [(mode = (1,), port = 1, current = 0.00565e-6)]
        ws = 2*pi*[4.5e9, 4.8e9, 5.0e9]
        T = 0.1
        blk = ScatteringParameters(ComplexF64[0.1 0.2im; 0.8 0.1]; zref = 50.0)
        circuit = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:x, 1, 4, blk),
            (:y, 4, 2, blk),
            (:cc, 2, 3, Capacitor(100e-15)), (:jj, 3, 0, JosephsonJunction(1000e-12)),
            (:cj, 3, 0, Capacitor(1000e-15)), (:r, 3, 0, Resistor(2.0e4))])
        nl = hbnlsolve(wp, (8,), src, circuit; keyedarrays = false)
        psc = JC.compile(circuit)
        sf = JC.truncfreqs(JC.calcfreqsdft((4,)); dc = true, odd = true,
            even = false, maxintermodorder = Inf)
        host = JC.hblinsolve(ws, psc, Dict{Any,Any}(), sf; nonlinear = nl,
            keyedarrays = false, returnSnoise = true, returnCnoise = true,
            temperature = T)
        d = JC.hblinsolve(ws, psc, Dict{Any,Any}(), sf; nonlinear = nl,
            debuglsys = true, temperature = T)
        @test any(w -> any(<(0), w .+ d.wpumpmodes), ws)
        ssys = d.lsys.scattering
        noiseplan = JC.planscatteringnoise(ssys)
        plan = JC.plandevicenoise(d.nodeindices, d.componenttypes,
            d.noiseportimpedanceindices,
            [d.vvn[j] for j in d.noiseportimpedanceindices], d.Nmodes, CPU())
        blockplan = JC.plandeviceblocknoise(ssys, noiseplan, d.Nmodes, CPU())
        providers = JC.plandeviceproviders(ssys, length(ws), CPU(),
            d.wpumpmodes, ssys.scale)
        temperatures = JC.noisechanneltemperatures(psc,
            d.noiseportimpedanceindices, noiseplan, ssys, T)
        A = copy(d.lsys.Asparse)
        adjoint = map(ws) do w
            JC.assemblesystemmatrix!(A, d.lsys, w .+ d.wpumpmodes)
            Matrix(sparse(transpose(A)) \ Matrix(d.bnm))
        end
        nrhs = size(d.bnm, 2)
        noise = JC.devicenoise(plan, blockplan, providers, i -> adjoint[i],
            nrhs, d.wpumpmodes, ws, true, temperatures)
        sources = JC.portsourcecurrents(d.bnm, d.portindices, d.nodeindices,
            d.Nmodes)
        inputwave = zeros(ComplexF64, nrhs)
        for i in eachindex(ws)
            wmodes = ws[i] .+ d.wpumpmodes
            JC.calcinputwaves!(inputwave, sources, d.portindices,
                d.portimpedances, d.componenttypes, wmodes)
            Snoise = similar(host.Snoise[:, :, i])
            Cnoise = similar(host.Cnoise[:, :, i])
            reduction = noise(i, inputwave, Snoise, Cnoise)
            @test isapprox(Snoise, host.Snoise[:, :, i]; rtol = 1e-12)
            @test isapprox(Cnoise, host.Cnoise[:, :, i]; rtol = 1e-12)
            @test isapprox(JC.calcqe(host.S[:, :, i], reduction),
                host.QE[:, :, i]; rtol = 1e-12)
            @test isapprox(JC.calccm(host.S[:, :, i], wmodes, reduction),
                host.CM[:, i]; rtol = 1e-12)
        end
    end

    @testset "the shared impedance" begin
        # the power wave kernels and the host both go through `impedance`,
        # so there is one implementation; these pin its behaviour
        @test JosephsonCircuits.impedancecode(:R) == 1
        @test JosephsonCircuits.impedancecode(:C) == 2
        @test JosephsonCircuits.impedancecode(:L) == 3
        @test_throws ErrorException JosephsonCircuits.impedancecode(:Lj)
        for c in (50.0, 1.0e-12, 3.0e-10, 50.0 - 2.0im)
            for w in (2*pi*5e9, -2*pi*3e9)
                # a negative frequency conjugates the stored value
                @test JosephsonCircuits.calcimpedance(c, :R, w) ==
                    (real(w) >= 0 ? c : conj(c)) + 0.0im
                @test JosephsonCircuits.calcimpedance(c, :C, w) ==
                    1/(im*w*(real(w) >= 0 ? c : conj(c)))
                @test JosephsonCircuits.calcimpedance(c, :L, w) ==
                    im*w*(real(w) >= 0 ? c : conj(c))
            end
        end
        # the wave scale is singular for a purely reactive impedance, which a
        # dissipative noise port never has, and zero at zero frequency
        @test JosephsonCircuits.portwavescale(complex(50.0), 0.0) == 0
        @test JosephsonCircuits.portwavescale(complex(50.0), 2*pi*5e9) ==
            1/sqrt(Complex(50.0))/sqrt(2*pi*5e9)
    end

    @testset "the column equilibration and the solutions it scales" begin
        # the grouping is the stored entries of each column of the structure
        # handed to the solver, which is compressed sparse row, so a column
        # index array with the entries of a column scattered through it
        for (m, dens) in ((40, 0.1), (17, 0.4), (6, 1.0))
            A = sprandn(m, m, dens) + I
            colind = rowvals(sparse(transpose(A)))
            order, segptr = JosephsonCircuits.groupstoredcolumns(colind, m)
            @test sort(order) == collect(1:length(colind))
            @test segptr[1] == 1 && segptr[end] == length(colind) + 1
            for j in 1:m
                @test all(colind[order[q]] == j
                    for q in segptr[j]:(segptr[j+1]-1))
            end
        end

        # a batch of systems whose rows and columns are orders of magnitude
        # apart: each comes out with a largest entry of one in every column
        # and none larger, and the solution of the scaled system, from the
        # right-hand side the equilibration writes, unscales to one of the
        # original
        n, nb, nrhs = 24, 3, 2
        A = sprandn(ComplexF64, n, n, 0.2) + I
        At = sparse(transpose(A))
        rowptr, colind = SparseArrays.getcolptr(At), rowvals(At)
        eq = JosephsonCircuits.planequilibration(rowptr, colind, n, nb, CPU())
        nzval = Matrix{ComplexF64}(undef, nnz(A), nb)
        Rs = [Diagonal(exp10.(range(6, -6, length = n))/k) for k in 1:nb]
        Ds = [Diagonal(exp10.(range(-6, 6, length = n))*k) for k in 1:nb]
        for k in 1:nb
            nzval[:,k] .= nonzeros(sparse(transpose(Rs[k]*A*Ds[k])))
        end
        B0 = randn(ComplexF64, n, nrhs)
        B = Array{ComplexF64}(undef, n, nrhs, nb)
        JosephsonCircuits.equilibrate!(nzval, B, B0, eq, CPU())
        # the structure is compressed sparse row, so read as compressed
        # sparse column it is the transpose
        scaled(k) = sparse(transpose(SparseMatrixCSC(n, n, rowptr, colind,
            nzval[:,k])))
        X = Array{ComplexF64}(undef, n, nrhs, nb)
        for k in 1:nb
            M = scaled(k)
            @test all(j -> maximum(abs, M[:, j]) ≈ 1, 1:n)
            @test maximum(abs, M) <= 1 + 1e-12
            X[:,:,k] .= M\B[:,:,k]
        end
        JosephsonCircuits.unscalesolution!(X, eq, CPU())
        for k in 1:nb
            @test X[:,:,k] ≈ (Rs[k]*A*Ds[k])\B0
        end

        # a row or a column of zeros is left alone rather than divided by
        # nothing
        Z = spzeros(ComplexF64, 3, 3)
        Z[1,1] = 2.0; Z[3,3] = 4.0
        Zt = sparse(transpose(Z))
        eqz = JosephsonCircuits.planequilibration(SparseArrays.getcolptr(Zt),
            rowvals(Zt), 3, 1, CPU())
        nz = reshape(collect(nonzeros(Zt)), :, 1)
        JosephsonCircuits.equilibrate!(nz, zeros(ComplexF64, 3, 1, 1),
            ones(ComplexF64, 3, 1), eqz, CPU())
        @test eqz.rowscale[:,1] == [2.0, 1.0, 4.0]
        @test eqz.scale[:,1] == [1.0, 1.0, 1.0]
    end

    @testset "the batch size cap avoids two cuDSS faults at sixteen" begin
        # cuDSS returns silently wrong solutions from a uniform batch of
        # sixteen or more systems with six or more right hand sides each,
        # and takes about eight times as long for a batch of sixteen as for
        # one of fifteen at every right hand side count, so the cap is
        # fifteen whatever the caller asks for
        for nrhs in (1, 2, 5, 6, 12, 64)
            @test JosephsonCircuits.uniformbatchlimit(nrhs) == 15
        end
    end
end
