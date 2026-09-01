using Symbolics
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

@testset verbose=true "devicelinsolve" begin

    # a chain with a Josephson junction per cell, two ports, and a resistor at
    # each end, so the linearized system has a pump modulation contribution,
    # a modified nodal analysis augmentation from the promoted port resistors,
    # and every linear term matrix populated.
    function buildcase(Nmod, wp, Npump)
        @variables Rleft Rright Cg Lj Cj
        circuit = Tuple{String,String,String,Num}[]
        push!(circuit,("P1_0","1","0",1))
        push!(circuit,("R1_0","1","0",Rleft))
        push!(circuit,("C1_0","1","0",Cg/2))
        push!(circuit,("Lj1_2","1","2",Lj))
        push!(circuit,("C1_2","1","2",Cj))
        for j in 2:8
            push!(circuit,("C$(j)_0","$(j)","0",Cg))
            push!(circuit,("Lj$(j)_$(j+1)","$(j)","$(j+1)",Lj))
            push!(circuit,("C$(j)_$(j+1)","$(j)","$(j+1)",Cj))
        end
        push!(circuit,("C9_0","9","0",Cg/2))
        push!(circuit,("R9_0","9","0",Rright))
        push!(circuit,("P9_0","9","0",2))
        circuitdefs = Dict(Lj => JosephsonCircuits.IctoLj(1e-6), Cg => 45e-15,
            Cj => 55e-15, Rleft => 50.0, Rright => 50.0)
        nl = JosephsonCircuits.hbnlsolve(wp, Npump,
            [(mode=ntuple(i->i==1 ? 1 : 0, length(wp)), port=1, current=1e-6)],
            circuit, circuitdefs; keyedarrays=false)
        psc = JosephsonCircuits.compile(circuit)
        cg = JosephsonCircuits.calccircuitgraph(psc)
        sf = JosephsonCircuits.removeconjfreqs(JosephsonCircuits.truncfreqs(
            JosephsonCircuits.calcfreqsrdft(Nmod); dc=true, odd=true,
            even=false, maxintermodorder=Inf))
        return psc, cg, circuitdefs, sf, nl
    end

    @testset "the sweep assembly matches the host assembler" begin

        # single tone, and two tone, the latter having negative mode
        # frequencies which exercise the conjugation of the stored values
        for (wp, Npump, Nmod) in (((2*pi*5e9,), (6,), (4,)),
                                  ((2*pi*5e9, 2*pi*3e9), (2,2), (2,2)))
            psc, cg, circuitdefs, sf, nl = buildcase(Nmod, wp, Npump)
            # off the integer GHz grid, where a signal frequency plus a
            # mode frequency would land on zero. The lowest is below the
            # magnitude of the two tone case's negative mode offset, so the
            # sweep reaches negative total mode frequencies there.
            ws = 2*pi*[0.43e9, 1.37e9, 3.11e9, 6.61e9, 9.29e9, 11.83e9]
            d = JosephsonCircuits.hblinsolve(ws, psc, cg, circuitdefs, sf;
                nonlinear=nl, debuglsys=true)
            lsys = d.lsys
            @test JosephsonCircuits.cansweepondevice(lsys)
            # the two tone case must actually reach negative mode
            # frequencies, or the conjugation of the stored values, which is
            # the only branch in the assembly kernel, goes untested
            if length(wp) > 1
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
        @variables w Rleft Cc Lj Cj
        circuit = Tuple{String,String,String,Num}[]
        push!(circuit,("P1","1","0",1))
        push!(circuit,("R1","1","0",Rleft))
        push!(circuit,("C1","1","2",Cc*(1+1e-18*w)))
        push!(circuit,("Lj1","2","0",Lj))
        push!(circuit,("C2","2","0",Cj))
        circuitdefs = Dict(Lj=>1000.0e-12, Cc=>100.0e-15, Cj=>1000.0e-15,
            Rleft=>50.0)
        psc = JosephsonCircuits.compile(circuit)
        cg = JosephsonCircuits.calccircuitgraph(psc)
        sf = JosephsonCircuits.removeconjfreqs(JosephsonCircuits.truncfreqs(
            JosephsonCircuits.calcfreqsrdft((2,)); dc=true, odd=true,
            even=false, maxintermodorder=Inf))
        d = JosephsonCircuits.hblinsolve(2*pi*[4.0e9, 5.0e9], psc, cg,
            circuitdefs, sf; symfreqvar=w, debuglsys=true)
        # a component value which depends on the frequency is not a constant
        # coefficient, so there is nothing to precompute and hblinsolve keeps
        # such a sweep on the host
        @test !JosephsonCircuits.cansweepondevice(d.lsys)
        @test_throws ArgumentError JosephsonCircuits.planfrequencysweep(d.lsys,
            CPU())
    end

    @testset "the gathered rows are the ones the ports read" begin
        wp = (2*pi*5e9,)
        psc, cg, circuitdefs, sf, nl = buildcase((4,), wp, (6,))
        ws = 2*pi*collect(range(2.13e9, 9.41e9, length=5))
        d = JosephsonCircuits.hblinsolve(ws, psc, cg, circuitdefs, sf;
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
        iw = zeros(ComplexF64, nrhs, nrhs); ow = similar(iw)
        iw2 = similar(iw); ow2 = similar(iw)
        JosephsonCircuits.calcinputoutput!(iw, ow, phin, d.bnm,
            d.portindices, d.portindices,
            d.portimpedances, d.portimpedances,
            d.nodeindices, d.componenttypes, wmodes, nothing)
        JosephsonCircuits.calcinputoutput!(iw2, ow2, masked, d.bnm,
            d.portindices, d.portindices,
            d.portimpedances, d.portimpedances,
            d.nodeindices, d.componenttypes, wmodes, nothing)
        @test iw == iw2
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
        psc, cg, circuitdefs, sf, nl = buildcase((4,), wp, (6,))
        ws = 2*pi*[0.43e9, 2.17e9, 7.91e9]
        d = JosephsonCircuits.hblinsolve(ws, psc, cg, circuitdefs, sf;
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

    @testset "the noise reduction is what QE and CM read of Snoise" begin
        # calcqe! and calccm! consume the noise scattering matrix only through
        # a sum over the noise index, so the reduced form must give the same
        # answers as the matrix
        for (nport, nnoise, m) in ((2, 5, 3), (3, 40, 4))
            np, nn = nport*m, nnoise*m
            S = randn(ComplexF64, np, np)
            Snoise = randn(ComplexF64, nn, np)
            w = randn(m); w[1] = abs(w[1]); w[end] = -abs(w[end])
            SnoiseT = transpose(Snoise)
            red = JosephsonCircuits.NoiseReduction(
                [sum(abs2, view(SnoiseT, i, :)) for i in 1:np],
                [sum(j -> abs2(SnoiseT[i,j])*sign(w[(j-1) % m + 1]), 1:nn)
                    for i in 1:np])

            qe1 = zeros(Float64, np, np); qe2 = similar(qe1)
            JosephsonCircuits.calcqe!(qe1, S, SnoiseT)
            JosephsonCircuits.calcqe!(qe2, S, red)
            @test isapprox(qe1, qe2; rtol = 1e-12)

            cm1 = zeros(Float64, np); cm2 = similar(cm1)
            JosephsonCircuits.calccm!(cm1, S, SnoiseT, w)
            JosephsonCircuits.calccm!(cm2, S, red, w)
            @test isapprox(cm1, cm2; rtol = 1e-12,
                norm = x -> maximum(abs, x))
        end
    end

    @testset "the shared impedance is what calcimpedance was" begin
        # the power wave kernels and the host both go through `impedance` now,
        # so there is one implementation; these pin its behaviour, which used
        # to be spelled out twice
        @test JosephsonCircuits.impedancecode(:R) == 1
        @test JosephsonCircuits.impedancecode(:C) == 2
        @test JosephsonCircuits.impedancecode(:L) == 3
        @test_throws ErrorException JosephsonCircuits.impedancecode(:Lj)
        for c in (50.0, 1.0e-12, 3.0e-10, 50.0 - 2.0im)
            for w in (2*pi*5e9, -2*pi*3e9)
                # a negative frequency conjugates the stored value
                @test JosephsonCircuits.calcimpedance(c, :R, w, nothing) ==
                    (real(w) >= 0 ? c : conj(c)) + 0.0im
                @test JosephsonCircuits.calcimpedance(c, :C, w, nothing) ==
                    1/(im*w*(real(w) >= 0 ? c : conj(c)))
                @test JosephsonCircuits.calcimpedance(c, :L, w, nothing) ==
                    im*w*(real(w) >= 0 ? c : conj(c))
            end
        end
        # the wave scale is singular for a purely reactive impedance, which a
        # dissipative noise port never has, and zero at zero frequency
        @test JosephsonCircuits.portwavescale(complex(50.0), 0.0) == 0
        @test JosephsonCircuits.portwavescale(complex(50.0), 2*pi*5e9) ==
            1/sqrt(Complex(50.0))/sqrt(2*pi*5e9)
    end

    @testset "the batch size cap avoids the cuDSS wrong answer" begin
        # cuDSS returns silently wrong solutions from a uniform batch of
        # sixteen or more systems with six or more right hand sides each
        @test JosephsonCircuits.uniformbatchlimit(6) == 15
        @test JosephsonCircuits.uniformbatchlimit(12) == 15
        @test JosephsonCircuits.uniformbatchlimit(5) > 15
        @test JosephsonCircuits.uniformbatchlimit(1) > 15
    end
end
