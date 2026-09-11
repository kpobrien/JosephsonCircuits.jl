using JosephsonCircuits, LinearAlgebra, SparseArrays, Random, Test
include("layoutreference.jl")


@testset verbose=true "the mode layout and its real form" begin


    # simple block approach to serve as an independent reference for sparse matrices
    function complex_to_real_ref(A, rl, cl; conj_input = false, rs = 1.0, cs = 1.0)
        I, J, V = Int[], Int[], Float64[]
        Ai, Av = JosephsonCircuits.SparseArrays.rowvals(A), JosephsonCircuits.SparseArrays.nonzeros(A)
        for j in 1:size(A,2), idx in JosephsonCircuits.SparseArrays.nzrange(A, j)
            i = Ai[idx]
            r0, wr = rl.ptr[i], JosephsonCircuits._width(rl, i)
            c0, wc = cl.ptr[j], JosephsonCircuits._width(cl, j)
            a = Av[idx] * ((wr == 1 ? rs : 1.0) * (wc == 1 ? cs : 1.0))
            d = conj_input ? -a : a
            push!(I, r0); push!(J, c0); push!(V, real(a))
            wr == 2 && (push!(I, r0+1); push!(J, c0); push!(V, imag(a)))
            if wc == 2
                push!(I, r0); push!(J, c0+1); push!(V, -imag(d))
                wr == 2 && (push!(I, r0+1); push!(J, c0+1); push!(V, real(d)))
            end
        end
        JosephsonCircuits.SparseArrays.sparse(I, J, V, rl.rdim, cl.rdim)
    end

    canon(xc, L) = [JosephsonCircuits._width(L, i) == 1 ? Complex(real(xc[i]), 0.0) : xc[i] for i in 1:L.dim]
    # A with q dropped wherever both modes are real
    function dropq(A, rl, cl)
        B = copy(A); Bi, Bv = JosephsonCircuits.SparseArrays.rowvals(B), JosephsonCircuits.SparseArrays.nonzeros(B)
        for j in 1:size(B,2), idx in JosephsonCircuits.SparseArrays.nzrange(B, j)
            JosephsonCircuits._width(rl, Bi[idx]) == 1 && JosephsonCircuits._width(cl, j) == 1 &&
                (Bv[idx] = Complex(real(Bv[idx]), 0.0))
        end
        B
    end
    shared(A, v) = JosephsonCircuits.SparseArrays.SparseMatrixCSC(size(A)..., A.colptr, A.rowval, v)   # same pattern object

    CASES = (([true,false,false,false,false], [true,false,false,false,false], 8, 8, 0.1),
                   ([true,false,false,false,false], [true,false,false,false,false], 12, 7, 0.2),
                   ([true,true,false], [false,false,false], 10, 6, 0.3),
                   ([false,false], [true,true], 9, 11, 0.25),
                   ([false], [false], 20, 15, 0.2),          # no real modes at all
                   ([true], [true], 13, 5, 0.4),             # every mode real
                   (rand(Bool, 6), rand(Bool, 6), 5, 4, 0.15))

    @testset "ModeLayout" begin
        L = JosephsonCircuits.ModeLayout([true,false,false,false,false], 10)
        @test L.rdim == 18 && L.nreal == 1 && L.nmodes == 5
        @test L.ptr == [1,2,4,6,8,10,11,13,15,17,19]
        @test L.inv == [1,2,2,3,3,4,4,5,5,6,7,7,8,8,9,9,10,10]
        @test JosephsonCircuits.ModeLayout([1], 5, 10).ptr == L.ptr
        @test JosephsonCircuits.ModeLayout([1,3], 4, 8).rdim == 12
        @test JosephsonCircuits.ModeLayout(falses(4), 8).rdim == 16
        @test JosephsonCircuits.ModeLayout(trues(4), 8).rdim == 8
        @test_throws DimensionMismatch JosephsonCircuits.ModeLayout([true,false], 7)
        @test_throws ArgumentError JosephsonCircuits.ModeLayout([6], 5, 10)
        @test_throws ArgumentError JosephsonCircuits.ModeLayout([1,1], 5, 10)
        # the compact width table agrees with ptr, at every size
        @test L.w isa BitVector
        Lb = JosephsonCircuits.ModeLayout([true,false,false,false,false], 5_000_000)
        @test all(JosephsonCircuits._rowwidth(Lb.w, i) == JosephsonCircuits._width(Lb, i) for i in 1:1000)
        for mask in ([true,true,false], falses(4), trues(3), rand(Bool, 7))
            M = JosephsonCircuits.ModeLayout(mask, length(mask) * 9)
            @test all(JosephsonCircuits._rowwidth(M.w, i) == JosephsonCircuits._width(M, i) for i in 1:M.dim)
        end
    end

    @testset "complex_to_real" begin
        for (rmask, cmask, mmul, nmul, p) in CASES
            m, n = length(rmask)*mmul, length(cmask)*nmul
            rl, cl = JosephsonCircuits.ModeLayout(rmask, m), JosephsonCircuits.ModeLayout(cmask, n)
            A = JosephsonCircuits.SparseArrays.sprand(ComplexF64, m, n, p)
            xc = canon(rand(ComplexF64, n), cl)
            for cj in (false, true)
                Ar = LayoutReference.complex_to_real(A, rl, cl; conj_input = cj)
                @test size(Ar) == (rl.rdim, cl.rdim)
                @test Ar == complex_to_real_ref(A, rl, cl; conj_input = cj)
                @test Ar * LayoutReference.complex_to_real(xc, cl.isreal) ≈ LayoutReference.complex_to_real(A * (cj ? conj(xc) : xc), rl.isreal)
                @test LayoutReference.is_complex_to_real_pattern(Ar, A, rl, cl)
                # in-place agrees and is exact
                S = LayoutReference.complex_to_real(A, rl, cl); fill!(JosephsonCircuits.SparseArrays.nonzeros(S), 0)
                @test LayoutReference.complex_to_real!(S, A, rl, cl; conj_input = cj) === S
                @test S == Ar
                # round trip, and idempotence of the second conversion
                Ac = LayoutReference.real_to_complex(Ar, rl, cl; conj_input = cj)
                @test Ac == dropq(A, rl, cl)
                @test LayoutReference.complex_to_real(Ac, rl, cl; conj_input = cj) == Ar
                C = copy(Ac); JosephsonCircuits.SparseArrays.nonzeros(C) .= 0
                @test LayoutReference.real_to_complex!(C, Ar, rl, cl; conj_input = cj) == Ac
            end
            # both conventions share one pattern
            @test JosephsonCircuits.SparseArrays.rowvals(LayoutReference.complex_to_real(A, rl, cl)) == JosephsonCircuits.SparseArrays.rowvals(LayoutReference.complex_to_real(A, rl, cl; conj_input = true))
            # index type
            @test LayoutReference.complex_to_real(A, rl, cl, Int32) isa JosephsonCircuits.SparseArrays.SparseMatrixCSC{Float64,Int32}
            @test LayoutReference.complex_to_real(A, rl, cl, Int32) == LayoutReference.complex_to_real(A, rl, cl)
        end
        @test LayoutReference.complex_to_real(JosephsonCircuits.SparseArrays.sprand(ComplexF32, 6, 6, 0.4), JosephsonCircuits.ModeLayout([true,false], 6),
                      JosephsonCircuits.ModeLayout([true,false], 6)) isa JosephsonCircuits.SparseArrays.SparseMatrixCSC{Float32,Int}
    end

    @testset "realscale" begin
        for (rmask, cmask, mmul, nmul, p) in CASES
            m, n = length(rmask)*mmul, length(cmask)*nmul
            rl, cl = JosephsonCircuits.ModeLayout(rmask, m), JosephsonCircuits.ModeLayout(cmask, n)
            A = JosephsonCircuits.SparseArrays.sprand(ComplexF64, m, n, p)
            base = LayoutReference.complex_to_real(A, rl, cl)
            @test LayoutReference.complex_to_real(A, rl, cl; realrowscale = 1, realcolscale = 1) == base
            s, xc = 0.5, canon(rand(ComplexF64, n), cl)
            # colscale scales the real modes of the input, rowscale those of the output
            xs = [JosephsonCircuits._width(cl, i) == 1 ? s*xc[i] : xc[i] for i in 1:n]
            @test LayoutReference.complex_to_real(A, rl, cl; realcolscale = s) * LayoutReference.complex_to_real(xc, cl.isreal) ≈ LayoutReference.complex_to_real(A * xs, rl.isreal)
            b  = A * xc
            bs = [JosephsonCircuits._width(rl, i) == 1 ? s*b[i] : b[i] for i in 1:m]
            @test LayoutReference.complex_to_real(A, rl, cl; realrowscale = s) * LayoutReference.complex_to_real(xc, cl.isreal) ≈ LayoutReference.complex_to_real(bs, rl.isreal)
            # against the reference, independently and together
            for (rr, cc) in ((2.0, 1.0), (1.0, 3.0), (2.0, 3.0), (0.0, 1.0), (1.0, 0.0))
                @test LayoutReference.complex_to_real(A, rl, cl; realrowscale = rr, realcolscale = cc) ==
                      complex_to_real_ref(A, rl, cl; rs = rr, cs = cc)
                @test LayoutReference.complex_to_real!(copy(base), A, rl, cl; realrowscale = rr, realcolscale = cc) ==
                      complex_to_real_ref(A, rl, cl; rs = rr, cs = cc)
            end
            Z = LayoutReference.complex_to_real(A, rl, cl; realrowscale = 0, realcolscale = 0)
            @test JosephsonCircuits.SparseArrays.nnz(Z) == JosephsonCircuits.SparseArrays.nnz(base) && Z.colptr == base.colptr && JosephsonCircuits.SparseArrays.rowvals(Z) == JosephsonCircuits.SparseArrays.rowvals(base)
            @test LayoutReference.complex_to_real(A, rl, cl; conj_input = true, realcolscale = s) ==
                  complex_to_real_ref(A, rl, cl; conj_input = true, cs = s)
        end
    end


    # @testset "allocation-free in-place" begin
    #     m = n = 400
    #     rl = JosephsonCircuits.ModeLayout([true,false,false,false,false], m)
    #     A = JosephsonCircuits.SparseArrays.sprand(ComplexF64, m, n, 0.02); B = JosephsonCircuits.SparseArrays.sprand(ComplexF64, m, n, 0.02)
    #     Bs = shared(A, rand(ComplexF64, JosephsonCircuits.SparseArrays.nnz(A)))
    #     Ar = LayoutReference.complex_to_real(A, rl, rl); C = copy(A)
    #     xc = rand(ComplexF64, n); xr = LayoutReference.complex_to_real(xc, rl.isreal)
    #     # warm up every method that is tested below
    #     LayoutReference.complex_to_real!(Ar, A, rl, rl); LayoutReference.complex_to_real!(Ar, A, rl, rl; conj_input = true, realcolscale = 0.5)
    #     LayoutReference.complex_to_real!(xr, xc, rl.isreal); LayoutReference.real_to_complex!(xc, xr, rl.isreal); LayoutReference.is_complex_to_real_pattern(Ar, A, rl, rl)

    #     @test (@allocated LayoutReference.complex_to_real!(Ar, A, rl, rl)) == 0
    #     @test (@allocated LayoutReference.complex_to_real!(Ar, A, rl, rl; conj_input = true, realcolscale = 0.5)) == 0
    #     @test (@allocated LayoutReference.real_to_complex!(C, Ar, rl, rl)) == 0
    #     @test (@allocated LayoutReference.is_complex_to_real_pattern(Ar, A, rl, rl)) == 0
    # end

    @testset "shape checks" begin
        m, n = 20, 15
        rl, cl = JosephsonCircuits.ModeLayout([true,false,false,false,false], m), JosephsonCircuits.ModeLayout([true,false,false,false,false], n)
        A = JosephsonCircuits.SparseArrays.sprand(ComplexF64, m, n, 0.2)
        Ar = LayoutReference.complex_to_real(A, rl, cl)
        @test_throws DimensionMismatch LayoutReference.complex_to_real(A, JosephsonCircuits.ModeLayout([true,false,false,false,false], 25), cl)
        @test_throws DimensionMismatch LayoutReference.complex_to_real!(JosephsonCircuits.SparseArrays.spzeros(rl.rdim+1, cl.rdim), A, rl, cl)
        @test_throws DimensionMismatch LayoutReference.real_to_complex!(A, JosephsonCircuits.SparseArrays.spzeros(rl.rdim, cl.rdim+2), rl, cl)
        @test_throws DimensionMismatch LayoutReference.real_to_complex(JosephsonCircuits.SparseArrays.spzeros(rl.rdim, cl.rdim+2), rl, cl)
        @test !LayoutReference.is_complex_to_real_pattern(LayoutReference.complex_to_real(A + JosephsonCircuits.SparseArrays.sprand(ComplexF64, m, n, 0.2), rl, cl), A, rl, cl)
        @test !LayoutReference.is_complex_to_real_pattern(LayoutReference.complex_to_real(A, rl, cl), A + JosephsonCircuits.SparseArrays.sprand(ComplexF64, m, n, 0.2), rl, cl)
    end

    @testset "mask entry points" begin
        mask = [true,false,false,false,false]
        for (m, n) in ((20, 20), (30, 15), (5, 25))
            A = JosephsonCircuits.SparseArrays.sprand(ComplexF64, m, n, 0.3)
            B = JosephsonCircuits.SparseArrays.sprand(ComplexF64, m, n, 0.3)
            rl, cl = JosephsonCircuits.ModeLayout(mask, m), JosephsonCircuits.ModeLayout(mask, n)
            @test LayoutReference.complex_to_real(A, mask) == LayoutReference.complex_to_real(A, rl, cl)
            @test LayoutReference.complex_to_real(A, mask, Int32) == LayoutReference.complex_to_real(A, rl, cl, Int32)
            @test LayoutReference.complex_to_real(A, mask; conj_input = true, realcolscale = 0.5) ==
                  LayoutReference.complex_to_real(A, rl, cl; conj_input = true, realcolscale = 0.5)
            @test LayoutReference.complex_to_real!(LayoutReference.complex_to_real(A, mask), A, mask; realrowscale = 2) ==
                  LayoutReference.complex_to_real(A, rl, cl; realrowscale = 2)
            Ar = LayoutReference.complex_to_real(A, mask)
            @test LayoutReference.real_to_complex(Ar, mask) == LayoutReference.real_to_complex(Ar, rl, cl)
            @test LayoutReference.real_to_complex!(copy(A), Ar, mask) == LayoutReference.real_to_complex(Ar, rl, cl)
            @test LayoutReference.is_complex_to_real_pattern(Ar, A, mask)
            # the square case shares one layout object between the two axes
            pair = LayoutReference._layouts(mask, m, n)
            m == n ? (@test pair[1] === pair[2]) : (@test pair[1] !== pair[2])
        end
        # dimensions must still divide
        @test_throws DimensionMismatch LayoutReference.complex_to_real(JosephsonCircuits.SparseArrays.sprand(ComplexF64, 7, 10, 0.3), mask)
        @test_throws DimensionMismatch LayoutReference.real_to_complex(JosephsonCircuits.SparseArrays.spzeros(10, 18), mask)  # 10 % 9 != 0
    end


    #  vectors and dense matrices

    disr(mask, i) = mask[(i - 1) % length(mask) + 1]
    dwid(mask, i) = disr(mask, i) ? 1 : 2

    # simple block approach to serve as an independent reference for dense matrices
    function dref(A, rm, cm; conj_input = false, rs = 1.0, cs = 1.0)
        m, n = size(A)
        Ar = zeros(Float64, JosephsonCircuits.realdim(m, rm), JosephsonCircuits.realdim(n, cm))
        c0 = 1
        for j in 1:n
            wc = dwid(cm, j)
            r0 = 1
            for i in 1:m
                wr = dwid(rm, i)
                a = A[i,j] * ((wr == 1 ? rs : 1.0) * (wc == 1 ? cs : 1.0))
                d = conj_input ? -a : a
                Ar[r0, c0] = real(a)
                wr == 2 && (Ar[r0+1, c0] = imag(a))
                if wc == 2
                    Ar[r0, c0+1] = -imag(d)
                    wr == 2 && (Ar[r0+1, c0+1] = real(d))
                end
                r0 += wr
            end
            c0 += wc
        end
        Ar
    end
    dcanon(x, mask) = [disr(mask, i) ? Complex(real(x[i]), 0.0) : x[i] for i in eachindex(x)]
    ddropq(A, rm, cm) = [disr(rm, i) && disr(cm, j) ? Complex(real(A[i,j]), 0.0) : A[i,j]
                        for i in 1:size(A,1), j in 1:size(A,2)]

    DCASES = (([true,false,false,false,false], [true,false,false,false,false], 8, 8),
                   ([true,false,false,false,false], [true,false,false,false,false], 12, 7),
                   ([true,true,false], [false,false,false], 10, 6),
                   ([false,false], [true,true], 9, 11),
                   ([false], [false], 20, 15),               # no real modes
                   ([true], [true], 13, 5),                  # every mode real
                   (rand(Bool, 6), rand(Bool, 6), 5, 4))

    @testset "dense: realdim / complexdim" begin
        @test JosephsonCircuits.realdim(10, [true,false,false,false,false]) == 18
        @test JosephsonCircuits.realdim(8, falses(4)) == 16
        @test JosephsonCircuits.realdim(8, trues(4)) == 8
        @test JosephsonCircuits.complexdim(JosephsonCircuits.realdim(30, [true,false,false]), [true,false,false]) == 30
        @test_throws DimensionMismatch JosephsonCircuits.realdim(7, [true,false])
        @test_throws DimensionMismatch JosephsonCircuits.complexdim(7, [true,false,false])
        @test_throws ArgumentError JosephsonCircuits.realdim(0, Bool[])
    end

    @testset "vectors" begin
        L = [true,false,false,false,false]
        @test LayoutReference.complex_to_real(ComplexF64[1, 2+3im, 4+5im, 6+7im, 8+9im], L) == Float64[1,2,3,4,5,6,7,8,9]
        for (mask, d) in ((( [true,false,false,false,false]), 20), ([false], 7), ([true], 6),
                          ([true,true,false], 12), (rand(Bool, 8), 32))
            xr = rand(JosephsonCircuits.realdim(d, mask))
            @test LayoutReference.complex_to_real(LayoutReference.real_to_complex(xr, mask), mask) == xr
            xc = rand(ComplexF64, d)
            @test LayoutReference.real_to_complex(LayoutReference.complex_to_real(xc, mask), mask) == dcanon(xc, mask)
            @test LayoutReference.real_to_complex!(fill(ComplexF64(NaN, NaN), d), xr, mask) == LayoutReference.real_to_complex(xr, mask)
            # conj_input flips the stored imaginary parts, scale hits the real modes
            @test LayoutReference.complex_to_real(xc, mask; conj_input = true) == LayoutReference.complex_to_real(conj(xc), mask)
            s = 0.5
            xs = [disr(mask, i) ? s*xc[i] : xc[i] for i in 1:d]
            @test LayoutReference.complex_to_real(xc, mask; realscale = s) == LayoutReference.complex_to_real(xs, mask)
            @test LayoutReference.real_to_complex(LayoutReference.complex_to_real(xc, mask), mask; realscale = s) == dcanon(xs, mask)
            @test_throws DimensionMismatch LayoutReference.complex_to_real!(Vector{Float64}(undef, length(xr)+1), xc, mask)
        end
    end

    @testset "dense: complex_to_real" begin
        for (rm, cm, mmul, nmul) in DCASES
            m, n = length(rm)*mmul, length(cm)*nmul
            A = rand(ComplexF64, m, n)
            xc = dcanon(rand(ComplexF64, n), cm)
            for cj in (false, true)
                Ar = LayoutReference.complex_to_real(A, rm, cm; conj_input = cj)
                @test size(Ar) == (JosephsonCircuits.realdim(m, rm), JosephsonCircuits.realdim(n, cm))
                @test Ar == dref(A, rm, cm; conj_input = cj)
                @test Ar * LayoutReference.complex_to_real(xc, cm) ≈ LayoutReference.complex_to_real(A * (cj ? conj(xc) : xc), rm)
                S = fill(NaN, size(Ar))
                @test LayoutReference.complex_to_real!(S, A, rm, cm; conj_input = cj) === S
                @test S == Ar
                Ac = LayoutReference.real_to_complex(Ar, rm, cm; conj_input = cj)
                @test Ac == ddropq(A, rm, cm)
                @test LayoutReference.complex_to_real(Ac, rm, cm; conj_input = cj) == Ar
                C = fill(ComplexF64(NaN, NaN), m, n)
                @test LayoutReference.real_to_complex!(C, Ar, rm, cm; conj_input = cj) == Ac
            end
            @test LayoutReference.complex_to_real(Float32.(real(A)) .+ 0im .|> ComplexF32, rm, cm) isa Matrix{Float32}
        end
    end

    @testset "dense: realscale" begin
        for (rm, cm, mmul, nmul) in DCASES
            m, n = length(rm)*mmul, length(cm)*nmul
            A = rand(ComplexF64, m, n)
            base = LayoutReference.complex_to_real(A, rm, cm)
            @test LayoutReference.complex_to_real(A, rm, cm; realrowscale = 1, realcolscale = 1) == base
            s, xc = 0.5, dcanon(rand(ComplexF64, n), cm)
            xs = [disr(cm, i) ? s*xc[i] : xc[i] for i in 1:n]
            @test LayoutReference.complex_to_real(A, rm, cm; realcolscale = s) * LayoutReference.complex_to_real(xc, cm) ≈ LayoutReference.complex_to_real(A * xs, rm)
            b  = A * xc
            bs = [disr(rm, i) ? s*b[i] : b[i] for i in 1:m]
            @test LayoutReference.complex_to_real(A, rm, cm; realrowscale = s) * LayoutReference.complex_to_real(xc, cm) ≈ LayoutReference.complex_to_real(bs, rm)
            for (rr, cc) in ((2.0, 1.0), (1.0, 3.0), (2.0, 3.0), (0.0, 1.0), (1.0, 0.0))
                @test LayoutReference.complex_to_real(A, rm, cm; realrowscale = rr, realcolscale = cc) ==
                      dref(A, rm, cm; rs = rr, cs = cc)
            end
            @test LayoutReference.complex_to_real(A, rm, cm; conj_input = true, realcolscale = s) ==
                  dref(A, rm, cm; conj_input = true, cs = s)
        end
    end


    # @testset "dense: allocation-free in-place" begin
    #     rm = [true,false,false,false,false]
    #     m, n = 200, 200
    #     A = rand(ComplexF64, m, n); B = rand(ComplexF64, m, n)
    #     xc = rand(ComplexF64, n); xr = LayoutReference.complex_to_real(xc, rm)
    #     LayoutReference.complex_to_real!(Ar, A, rm, rm); LayoutReference.complex_to_real!(Ar, A, rm, rm; conj_input = true, realcolscale = 0.5)
    #     LayoutReference.complex_to_real!(xr, xc, rm); LayoutReference.real_to_complex!(xc, xr, rm)
    #     @test (@allocated LayoutReference.complex_to_real!(Ar, A, rm, rm)) == 0
    #     @test (@allocated LayoutReference.complex_to_real!(Ar, A, rm, rm; conj_input = true, realcolscale = 0.5)) == 0
    #     @test (@allocated LayoutReference.real_to_complex!(C, Ar, rm, rm)) == 0
    #     @test (@allocated LayoutReference.complex_to_real!(xr, xc, rm)) == 0
    #     @test (@allocated LayoutReference.real_to_complex!(xc, xr, rm)) == 0
    # end

    @testset "dense: shape checks" begin
        rm = [true,false,false,false,false]
        m, n = 20, 15
        A = rand(ComplexF64, m, n)
        @test_throws DimensionMismatch LayoutReference.complex_to_real(A, rm, [true,false])          # 15 % 2 != 0
        @test_throws DimensionMismatch LayoutReference.complex_to_real!(zeros(JosephsonCircuits.realdim(m,rm)+1, JosephsonCircuits.realdim(n,rm)), A, rm, rm)
        @test_throws DimensionMismatch LayoutReference.real_to_complex!(A, zeros(JosephsonCircuits.realdim(m,rm), JosephsonCircuits.realdim(n,rm)+3), rm, rm)
        # views work
        Abig = rand(ComplexF64, m+4, n+4)
        @test LayoutReference.complex_to_real(view(Abig, 1:m, 1:n), rm, rm) == LayoutReference.complex_to_real(Abig[1:m, 1:n], rm, rm)
    end

    @testset "real_to_complex rejects an incomplete row mode" begin
        # `complex_to_real` always stores both real slots of a complex row
        # mode. If a caller passes a pattern where the second slot is stored
        # but the first is not, the imaginary part has no entry to be folded
        # into: the write used to land on an unrelated earlier entry, or on
        # index 0 of a zero-length nzval, which `@inbounds` did not catch.
        mask = [false, false]
        Ar = JosephsonCircuits.SparseArrays.sparse([2], [1], [1.0], 4, 4)
        @test_throws ArgumentError LayoutReference.real_to_complex(Ar, mask)

        # the second slot of the second row mode, with the first absent
        Ar2 = JosephsonCircuits.SparseArrays.sparse([1, 4], [1, 1], [1.0, 2.0], 4, 4)
        @test_throws ArgumentError LayoutReference.real_to_complex(Ar2, mask)

        # a well formed pattern still round trips
        A = JosephsonCircuits.SparseArrays.sparse([1, 2], [1, 2], ComplexF64[1.0+2.0im, 3.0-1.0im], 2, 2)
        Arok = LayoutReference.complex_to_real(A, mask)
        @test LayoutReference.real_to_complex(Arok, mask) == A
    end

end

@testset "the package's conversions agree with the reference family" begin
    JC = JosephsonCircuits
    rng = Random.default_rng()
    for isreal in ([true, false, false], [false, false], [true])
        nm = length(isreal)
        nnodes = 5
        xc = randn(rng, ComplexF64, nm*nnodes)
        for t in 1:nm
            isreal[t] && (xc[t:nm:end] .= real.(xc[t:nm:end]))
        end
        xr = JC.complex_to_real(xc, isreal)
        @test xr == LayoutReference.complex_to_real(xc, isreal)
        @test JC.real_to_complex(xr, isreal) == LayoutReference.real_to_complex(xr, isreal)
        @test JC.real_to_complex(xr, isreal) == xc
        A = sprandn(rng, ComplexF64, nm*nnodes, nm*nnodes, 0.3) + I
        L = JC.ModeLayout(isreal, nm*nnodes)
        Ar = JC.complex_to_real(A, L, L)
        @test Ar == LayoutReference.complex_to_real(A, L, L)
        @test Ar*xr ≈ JC.complex_to_real(A*xc, isreal)
        @test JC.complex_to_real(A, L, L, Int32) isa SparseMatrixCSC{Float64,Int32}
    end

    @testset "the gather and scatter kernels are inverse permutations" begin
        # the device side of the canonical layout's permuted copies, run on
        # the CPU backend here
        v = randn(Random.default_rng(), 9)
        index = [4, 1, 9, 2, 7]
        got = zeros(2, 3)
        JosephsonCircuits.gathervalues!(got, v, reshape(vcat(index, 5), 2, 3))
        @test vec(got) == v[vcat(index, 5)]
        @test_throws DimensionMismatch JosephsonCircuits.gathervalues!(
            zeros(1, 1), v, index)
        w = fill(NaN, 9)
        JosephsonCircuits.scattervalues!(w, v[index], index)
        @test w[index] == v[index] && all(isnan, w[setdiff(1:9, index)])
        @test_throws DimensionMismatch JosephsonCircuits.scattervalues!(
            w, v[index], index[1:end-1])
    end
end

# The canonical state layout: which entries of the state a mode owns,
# and the gather and scatter between the full state and the canonical
# one the solvers iterate on.
@testset verbose=true "canonical state layout" begin
    JC = JosephsonCircuits

    @testset "layout of a state with one zero frequency mode" begin
        # 80 nodes, 77 modes, mode 1 is the zero frequency one and the only
        # self conjugate mode: the shape of the two tone hard case
        isdc = [i == 1 for i in 1:77]
        ml = JC.ModeLayout(isdc, 80*77)
        L = JC.compositelayout(ml, isdc)

        @test ml.rdim == 12240
        @test L.ndc == 80                  # one per node
        @test JC.canonicaldim(L) == ml.rdim
        @test L.nvdc == 0
        @test JC.isinternal(L)
        @test JC.nwindow(L) == 80

        # each node contributes 1 + 2*76 = 153 internal entries with the
        # zero frequency one first, so the window names exactly those
        @test L.dcpos == [k*153 + 1 for k in 0:79]
        @test JC.windowindices(L) == L.dcpos
        @test [JC.windowindex(L, k) for k in 1:80] == L.dcpos
    end

    @testset "a self conjugate mode which is not zero frequency stays alternating" begin
        # mode 40 is self conjugate (a Nyquist mode) but not direct current
        isreal = [i == 1 || i == 40 for i in 1:77]
        isdc   = [i == 1 for i in 1:77]
        ml = JC.ModeLayout(isreal, 80*77)
        L = JC.compositelayout(ml, isdc)
        @test L.ndc == 80                  # only the zero frequency mode
        @test JC.canonicaldim(L) == ml.rdim
        @test JC.isinternal(L)
    end

    @testset "the mode tuples pick out the zero frequency mode" begin
        modes = [(0,0), (1,0), (0,1), (1,1)]
        ml = JC.ModeLayout([all(iszero, m) for m in modes], 5*4)
        @test JC.compositelayout(ml, modes).ndc == 5
    end

    @testset "the internal state is the first block, and the voltages follow" begin
        isdc = [i == 1 for i in 1:9]
        ml = JC.ModeLayout(isdc, 7*9)
        L = JC.compositelayout(ml, isdc; nvdc = 2)
        @test JC.canonicaldim(L) == ml.rdim + 2
        @test JC.voltagerange(L) == (ml.rdim + 1):(ml.rdim + 2)
        @test JC.windowindices(L) == vcat(L.dcpos, ml.rdim + 1, ml.rdim + 2)
        @test !JC.isinternal(L)
        r = randn(ml.rdim)
        u = zeros(ml.rdim + 2); u[end-1:end] .= (3.0, 4.0)
        JC.gathercanonical!(u, r, L)
        @test u[1:ml.rdim] == r            # a copy, so bit exact
        @test u[end-1:end] == [3.0, 4.0]   # and the voltages are left alone
        back = similar(r)
        JC.scattercanonical!(back, u, L)
        @test back == r
        @test JC.internalpart(u, L) == r
    end

    @testset "rejected layouts" begin
        isreal = [i == 1 for i in 1:9]
        ml = JC.ModeLayout(isreal, 7*9)
        # a zero frequency mode has no imaginary part, so it must be self
        # conjugate; marking a conjugate pair as direct current is a bug
        @test_throws ArgumentError JC.compositelayout(ml, [i == 2 for i in 1:9])
        @test_throws DimensionMismatch JC.compositelayout(ml, trues(3))
        @test_throws ArgumentError JC.compositelayout(ml, isreal; nvdc = -1)

        L = JC.compositelayout(ml, isreal)
        @test_throws DimensionMismatch JC.gathercanonical!(zeros(3), zeros(ml.rdim), L)
        @test_throws DimensionMismatch JC.gathercanonical!(zeros(ml.rdim), zeros(3), L)
        @test_throws DimensionMismatch JC.scattercanonical!(zeros(ml.rdim), zeros(3), L)
    end

    @testset "residual and product are the same operator in the new basis" begin
        # a small two tone circuit, taken to a point that is not the origin
        circuit = Circuit(
            [:p1 => Port(1), :cc => Capacitor(100e-15),
             :jj => JosephsonJunction(1000e-12), :cj => Capacitor(1000e-15),
             :gnd => Ground()],
            [[(:p1, 1), (:cc, 1)],
             [(:cc, 2), (:jj, 1), (:cj, 1)],
             [(:p1, 2), (:jj, 2), (:cj, 2), (:gnd, 1)]])
        srcs = [(mode = (1,), port = 1, current = 0.5e-6)]
        # three iterations is a truncated solve on purpose: what this needs
        # is a point away from the origin, not the solution, so the solver
        # reports that it did not converge and that is asserted here rather
        # than printed
        s = @test_logs((:warn,), match_mode=:any,
            JC.hbnlsolve((2*pi*4.75e9,), (4,), srcs, circuit;
                dc = true, odd = true, even = true,
                returnoperatingpoint = true, iterations = 3))
        sys = s.operatingpoint.sys
        ml = s.operatingpoint.modelayout
        L = JC.compositelayout(ml, s.frequencies.modes)
        @test JC.isinternal(L)
        @test L.ndc > 0

        x = randn(L.rdim); v = randn(L.rdim)
        u = similar(x); vc = similar(v)
        JC.gathercanonical!(u, x, L); JC.gathercanonical!(vc, v, L)
        work = JC.CanonicalWork(L, x)

        # the residual in the new basis is the old one, gathered
        Fi = zeros(L.rdim)
        JC.setpoint!(sys, x); JC.residual!(Fi, sys)
        want = similar(Fi); JC.gathercanonical!(want, Fi, L)
        got = zeros(L.rdim); uu = copy(u)
        JC.canonicalresidual(
            (F, J, xx) -> (JC.setpoint!(sys, xx);
                isnothing(F) || JC.residual!(F, sys); nothing), work)(got, nothing, uu)
        @test got == want
        @test uu == u                      # the point is handed back unchanged

        # and so is the Jacobian vector product
        Ji = zeros(L.rdim)
        JC.setpoint!(sys, x); JC.jacobianvectorproduct!(Ji, sys, v)
        wantj = similar(Ji); JC.gathercanonical!(wantj, Ji, L)
        gotj = zeros(L.rdim)
        JC.setpoint!(sys, x)
        JC.canonicaljvp((o, vv) -> JC.jacobianvectorproduct!(o, sys, vv),
            work)(gotj, vc)
        @test gotj == wantj
    end

    @testset "the methods which can carry the block, and the one which cannot" begin
        circuit = Circuit(
            [:p1 => Port(1; Z0 = 50.0), :cc => Capacitor(100e-15),
             :jj => JosephsonJunction(1000e-12), :gnd => Ground()],
            [[(:p1, 1), (:cc, 1)], [(:cc, 2), (:jj, 1)],
             [(:p1, 2), (:jj, 2), (:gnd, 1)]])
        srcs = [(mode = (0,), port = 1, current = 1.0e-7)]

        # the two which solve the real system take the explicit block
        for m in (NewtonKrylov(), Newton())
            s = JC.hbnlsolve((2*pi*4.75e9,), (4,), srcs, circuit;
                dc = true, odd = true, method = m, rtol = 1e-12)
            @test s.solverinfo.converged
            @test !isnothing(s.dcnodevoltage)
        end

        # `:quasinewton` solves the complex holomorphic system, which has no
        # place for the real direct current unknowns, and says so
        @test_throws ArgumentError JC.hbnlsolve((2*pi*4.75e9,), (4,), srcs,
            circuit; dc = true, odd = true, method = QuasiNewton())
    end

    @testset "the assembled Jacobian survives a growing internal pattern" begin
        # The canonical Jacobian's pattern is the internal pattern under a
        # permutation plus the direct current block, and none of it moves.
        # It used to be rebuilt from the *values* at each point, and sparse
        # addition prunes exact zeros, so the pattern it produced was the
        # one the starting point happened to have. A junction driven hard
        # enough to develop harmonics fills in mode coupling entries which
        # were zero at the origin, the pattern grew, and the solve stopped
        # with the internal error which guarded against exactly that.
        #
        # The plan is built from the pattern alone, so it holds everywhere.
        # The check is that `:newton` reaches the same point `:newtonkrylov`
        # does, on a drive strong enough for the fill in to happen. The
        # drive is a few times the critical current, so the junction is
        # strongly nonlinear, but not so far past it that the operating
        # point is one of several and the two paths pick different ones.
        circuit = Circuit(
            [:p1 => Port(1; Z0 = 50.0), :cc => Capacitor(100e-15),
             :jj => JosephsonJunction(1000e-12), :cj => Capacitor(1000e-15)],
            [[(:p1,1),(:cc,1)], [(:cc,2),(:jj,1),(:cj,1)],
             [(:p1,2),(:jj,2),(:cj,2), Ground]])
        srcs = [(mode = (1,), port = 1, current = 1.0e-6),
                (mode = (0,), port = 1, current = 1.0e-7)]
        kw = (; dc = true, odd = true, even = true, keyedarrays = false,
              rtol = 1e-12)
        a = JC.hbnlsolve((2*pi*4.75e9,), (8,), srcs, circuit;
            kw..., method = Newton())
        b = JC.hbnlsolve((2*pi*4.75e9,), (8,), srcs, circuit;
            kw..., method = NewtonKrylov())
        @test a.solverinfo.converged
        @test b.solverinfo.converged
        # the scattering parameters and the direct current operating point
        # are the physical content and agree to roundoff
        @test maximum(abs, a.S .- b.S) < 1e-10
        @test isapprox(a.dcnodevoltage, b.dcnodevoltage; rtol = 1e-8)
        # the node fluxes may differ by a whole number of flux quanta, which
        # is the additive static gauge: a junction phase is defined modulo
        # 2*pi and the two methods are free to land on different branches
        turns = (a.nodeflux .- b.nodeflux) ./ (2*pi)
        @test all(x -> isapprox(x, round(real(x)); atol = 1e-8), turns)
    end

    @testset "the assembled Jacobian is the matrix free one" begin
        # A circuit with an explicit direct current block: a scattering
        # block which is a resistor, driven at zero frequency. The assembled
        # canonical Jacobian has to reproduce the matrix free product for
        # every unit vector, which is a sharper check than differencing the
        # residual and is what lets the direct solve methods use the same
        # formulation as the Krylov one.
        blk = ScatteringParameters(
            w -> JC.ABCDtoS(JC.ABCD_seriesZ(100.0 + 0im));
            nports = 2, grounded = true, noise = Lossless())
        c = Circuit(
            [:p1 => Port(1; Z0 = 1.0e9), :b => blk, :c1 => Capacitor(1e-12)],
            [[(:p1,1),(:b,1),(:c1,1)], [(:b,2), Ground],
             [(:p1,2),(:c1,2), Ground]])
        srcs = [(mode = (0,), port = 1, current = 1.0e-6)]
        d = JC.hbnlsolve((2*pi*5e9,), (1,), srcs, c, Dict{Any,Any}();
            dc = true, odd = true, keyedarrays = false, returnsystem = true)

        sys, ml, Nmodes = d.sys, d.modelayout, d.Nmodes
        Nnodes = length(d.dcplan.componentof) + 1
        tr = JC.transportrows(d.dcplan, d.bnmsource, Nmodes)
        L = JC.compositelayout(ml, d.frequencies.modes;
            nvdc = JC.nvoltages(tr))
        psc = JC.compile(c)
        Naux = ml.dim - d.Nnodal
        nauxsc = JC.countscatteringports(psc)*Nmodes
        ssys = JC.scatteringstampsystem(psc.scatteringblocks, Nmodes;
            auxoffset = d.Nnodal + Naux - nauxsc,
            Ntotal = d.Nnodal + Naux, scale = d.Lscale)
        br = JC.dcblockrows(ssys.blocks, d.dcplan.componentof, Nmodes,
            d.dcplan.modeindex, Nnodes - 1, d.Lscale)
        work = JC.CanonicalWork(L, zeros(ml.rdim); transport = tr,
            blockrows = br, nnodaldc = Nnodes - 1)
        @test L.nvdc > 0                    # the block is really explicit

        n = JC.canonicaldim(L)
        JC.setpoint!(sys, zeros(ml.rdim))
        jvp = JC.canonicaljvp(
            (o, v) -> JC.jacobianvectorproduct!(o, sys, v), work)
        free = zeros(n, n); col = zeros(n); out = zeros(n)
        for k in 1:n
            fill!(col, 0.0); col[k] = 1.0
            jvp(out, col); free[:,k] .= out
        end
        # the coupling the preconditioner subtracts is the Jacobian's own:
        # the columns of the direct current coordinates on the nodal zero
        # frequency rows, which is the resistor current the voltages drive
        # and the two terminals each block current enters
        H = JC.dccoupling(work)
        idx = JC.dcsubsystemindices(work)
        nodal = L.dcpos[1:size(H, 1)]
        @test free[nodal, idx] ≈ Matrix(H)
        # and no other zero frequency row sees the direct current
        # coordinates from outside the subsystem: the block is triangular
        below = setdiff(JC.windowindices(L), nodal, idx)
        @test all(iszero, free[below, idx])

        Jint = copy(d.Jr)
        JC.jacobian!(Jint, sys)
        assembled = JC.canonicaljacobian(Jint, work)
        @test size(assembled) == (n, n)
        @test Matrix(assembled) == free      # exactly, not to a tolerance
    end
end
