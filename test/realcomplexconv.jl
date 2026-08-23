using JosephsonCircuits, LinearAlgebra, Test


@testset verbose=true "realcomplexconv" begin


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
                Ar = JosephsonCircuits.complex_to_real(A, rl, cl; conj_input = cj)
                @test size(Ar) == (rl.rdim, cl.rdim)
                @test Ar == complex_to_real_ref(A, rl, cl; conj_input = cj)
                @test Ar * JosephsonCircuits.complex_to_real(xc, cl.isreal) ≈ JosephsonCircuits.complex_to_real(A * (cj ? conj(xc) : xc), rl.isreal)
                @test JosephsonCircuits.is_complex_to_real_pattern(Ar, A, rl, cl)
                # in-place agrees and is exact
                S = JosephsonCircuits.complex_to_real(A, rl, cl); fill!(JosephsonCircuits.SparseArrays.nonzeros(S), 0)
                @test JosephsonCircuits.complex_to_real!(S, A, rl, cl; conj_input = cj) === S
                @test S == Ar
                # round trip, and idempotence of the second conversion
                Ac = JosephsonCircuits.real_to_complex(Ar, rl, cl; conj_input = cj)
                @test Ac == dropq(A, rl, cl)
                @test JosephsonCircuits.complex_to_real(Ac, rl, cl; conj_input = cj) == Ar
                C = copy(Ac); JosephsonCircuits.SparseArrays.nonzeros(C) .= 0
                @test JosephsonCircuits.real_to_complex!(C, Ar, rl, cl; conj_input = cj) == Ac
            end
            # both conventions share one pattern
            @test JosephsonCircuits.SparseArrays.rowvals(JosephsonCircuits.complex_to_real(A, rl, cl)) == JosephsonCircuits.SparseArrays.rowvals(JosephsonCircuits.complex_to_real(A, rl, cl; conj_input = true))
            # index type
            @test JosephsonCircuits.complex_to_real(A, rl, cl, Int32) isa JosephsonCircuits.SparseArrays.SparseMatrixCSC{Float64,Int32}
            @test JosephsonCircuits.complex_to_real(A, rl, cl, Int32) == JosephsonCircuits.complex_to_real(A, rl, cl)
        end
        @test JosephsonCircuits.complex_to_real(JosephsonCircuits.SparseArrays.sprand(ComplexF32, 6, 6, 0.4), JosephsonCircuits.ModeLayout([true,false], 6),
                      JosephsonCircuits.ModeLayout([true,false], 6)) isa JosephsonCircuits.SparseArrays.SparseMatrixCSC{Float32,Int}
    end

    @testset "realscale" begin
        for (rmask, cmask, mmul, nmul, p) in CASES
            m, n = length(rmask)*mmul, length(cmask)*nmul
            rl, cl = JosephsonCircuits.ModeLayout(rmask, m), JosephsonCircuits.ModeLayout(cmask, n)
            A = JosephsonCircuits.SparseArrays.sprand(ComplexF64, m, n, p)
            base = JosephsonCircuits.complex_to_real(A, rl, cl)
            @test JosephsonCircuits.complex_to_real(A, rl, cl; realrowscale = 1, realcolscale = 1) == base
            s, xc = 0.5, canon(rand(ComplexF64, n), cl)
            # colscale scales the real modes of the input, rowscale those of the output
            xs = [JosephsonCircuits._width(cl, i) == 1 ? s*xc[i] : xc[i] for i in 1:n]
            @test JosephsonCircuits.complex_to_real(A, rl, cl; realcolscale = s) * JosephsonCircuits.complex_to_real(xc, cl.isreal) ≈ JosephsonCircuits.complex_to_real(A * xs, rl.isreal)
            b  = A * xc
            bs = [JosephsonCircuits._width(rl, i) == 1 ? s*b[i] : b[i] for i in 1:m]
            @test JosephsonCircuits.complex_to_real(A, rl, cl; realrowscale = s) * JosephsonCircuits.complex_to_real(xc, cl.isreal) ≈ JosephsonCircuits.complex_to_real(bs, rl.isreal)
            # against the reference, independently and together
            for (rr, cc) in ((2.0, 1.0), (1.0, 3.0), (2.0, 3.0), (0.0, 1.0), (1.0, 0.0))
                @test JosephsonCircuits.complex_to_real(A, rl, cl; realrowscale = rr, realcolscale = cc) ==
                      complex_to_real_ref(A, rl, cl; rs = rr, cs = cc)
                @test JosephsonCircuits.JosephsonCircuits.complex_to_real!(copy(base), A, rl, cl; realrowscale = rr, realcolscale = cc) ==
                      complex_to_real_ref(A, rl, cl; rs = rr, cs = cc)
            end
            Z = JosephsonCircuits.complex_to_real(A, rl, cl; realrowscale = 0, realcolscale = 0)
            @test JosephsonCircuits.SparseArrays.nnz(Z) == JosephsonCircuits.SparseArrays.nnz(base) && Z.colptr == base.colptr && JosephsonCircuits.SparseArrays.rowvals(Z) == JosephsonCircuits.SparseArrays.rowvals(base)
            @test JosephsonCircuits.JosephsonCircuits.complex_to_real(A, rl, cl; conj_input = true, realcolscale = s) ==
                  complex_to_real_ref(A, rl, cl; conj_input = true, cs = s)
        end
    end


    # @testset "allocation-free in-place" begin
    #     m = n = 400
    #     rl = JosephsonCircuits.ModeLayout([true,false,false,false,false], m)
    #     A = JosephsonCircuits.SparseArrays.sprand(ComplexF64, m, n, 0.02); B = JosephsonCircuits.SparseArrays.sprand(ComplexF64, m, n, 0.02)
    #     Bs = shared(A, rand(ComplexF64, JosephsonCircuits.SparseArrays.nnz(A)))
    #     Ar = JosephsonCircuits.complex_to_real(A, rl, rl); C = copy(A)
    #     xc = rand(ComplexF64, n); xr = JosephsonCircuits.complex_to_real(xc, rl.isreal)
    #     # warm up every method that is tested below
    #     JosephsonCircuits.complex_to_real!(Ar, A, rl, rl); JosephsonCircuits.complex_to_real!(Ar, A, rl, rl; conj_input = true, realcolscale = 0.5)
    #     JosephsonCircuits.complex_to_real!(xr, xc, rl.isreal); JosephsonCircuits.real_to_complex!(xc, xr, rl.isreal); JosephsonCircuits.is_complex_to_real_pattern(Ar, A, rl, rl)

    #     @test (@allocated JosephsonCircuits.complex_to_real!(Ar, A, rl, rl)) == 0
    #     @test (@allocated JosephsonCircuits.complex_to_real!(Ar, A, rl, rl; conj_input = true, realcolscale = 0.5)) == 0
    #     @test (@allocated JosephsonCircuits.real_to_complex!(C, Ar, rl, rl)) == 0
    #     @test (@allocated JosephsonCircuits.is_complex_to_real_pattern(Ar, A, rl, rl)) == 0
    # end

    @testset "shape checks" begin
        m, n = 20, 15
        rl, cl = JosephsonCircuits.ModeLayout([true,false,false,false,false], m), JosephsonCircuits.ModeLayout([true,false,false,false,false], n)
        A = JosephsonCircuits.SparseArrays.sprand(ComplexF64, m, n, 0.2)
        Ar = JosephsonCircuits.complex_to_real(A, rl, cl)
        @test_throws DimensionMismatch JosephsonCircuits.complex_to_real(A, JosephsonCircuits.ModeLayout([true,false,false,false,false], 25), cl)
        @test_throws DimensionMismatch JosephsonCircuits.complex_to_real!(JosephsonCircuits.SparseArrays.spzeros(rl.rdim+1, cl.rdim), A, rl, cl)
        @test_throws DimensionMismatch JosephsonCircuits.real_to_complex!(A, JosephsonCircuits.SparseArrays.spzeros(rl.rdim, cl.rdim+2), rl, cl)
        @test_throws DimensionMismatch JosephsonCircuits.real_to_complex(JosephsonCircuits.SparseArrays.spzeros(rl.rdim, cl.rdim+2), rl, cl)
        @test !JosephsonCircuits.is_complex_to_real_pattern(JosephsonCircuits.complex_to_real(A + JosephsonCircuits.SparseArrays.sprand(ComplexF64, m, n, 0.2), rl, cl), A, rl, cl)
        @test !JosephsonCircuits.is_complex_to_real_pattern(JosephsonCircuits.complex_to_real(A, rl, cl), A + JosephsonCircuits.SparseArrays.sprand(ComplexF64, m, n, 0.2), rl, cl)
    end

    @testset "mask entry points" begin
        mask = [true,false,false,false,false]
        for (m, n) in ((20, 20), (30, 15), (5, 25))
            A = JosephsonCircuits.SparseArrays.sprand(ComplexF64, m, n, 0.3)
            B = JosephsonCircuits.SparseArrays.sprand(ComplexF64, m, n, 0.3)
            rl, cl = JosephsonCircuits.ModeLayout(mask, m), JosephsonCircuits.ModeLayout(mask, n)
            @test JosephsonCircuits.complex_to_real(A, mask) == JosephsonCircuits.complex_to_real(A, rl, cl)
            @test JosephsonCircuits.complex_to_real(A, mask, Int32) == JosephsonCircuits.complex_to_real(A, rl, cl, Int32)
            @test JosephsonCircuits.complex_to_real(A, mask; conj_input = true, realcolscale = 0.5) ==
                  JosephsonCircuits.complex_to_real(A, rl, cl; conj_input = true, realcolscale = 0.5)
            @test JosephsonCircuits.complex_to_real!(JosephsonCircuits.complex_to_real(A, mask), A, mask; realrowscale = 2) ==
                  JosephsonCircuits.complex_to_real(A, rl, cl; realrowscale = 2)
            Ar = JosephsonCircuits.complex_to_real(A, mask)
            @test JosephsonCircuits.real_to_complex(Ar, mask) == JosephsonCircuits.real_to_complex(Ar, rl, cl)
            @test JosephsonCircuits.real_to_complex!(copy(A), Ar, mask) == JosephsonCircuits.real_to_complex(Ar, rl, cl)
            @test JosephsonCircuits.is_complex_to_real_pattern(Ar, A, mask)
            # the square case shares one layout object between the two axes
            pair = JosephsonCircuits._layouts(mask, m, n)
            m == n ? (@test pair[1] === pair[2]) : (@test pair[1] !== pair[2])
        end
        # dimensions must still divide
        @test_throws DimensionMismatch JosephsonCircuits.complex_to_real(JosephsonCircuits.SparseArrays.sprand(ComplexF64, 7, 10, 0.3), mask)
        @test_throws DimensionMismatch JosephsonCircuits.real_to_complex(JosephsonCircuits.SparseArrays.spzeros(10, 18), mask)  # 10 % 9 != 0
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
        @test JosephsonCircuits.complex_to_real(ComplexF64[1, 2+3im, 4+5im, 6+7im, 8+9im], L) == Float64[1,2,3,4,5,6,7,8,9]
        for (mask, d) in ((( [true,false,false,false,false]), 20), ([false], 7), ([true], 6),
                          ([true,true,false], 12), (rand(Bool, 8), 32))
            xr = rand(JosephsonCircuits.realdim(d, mask))
            @test JosephsonCircuits.complex_to_real(JosephsonCircuits.real_to_complex(xr, mask), mask) == xr
            xc = rand(ComplexF64, d)
            @test JosephsonCircuits.real_to_complex(JosephsonCircuits.complex_to_real(xc, mask), mask) == dcanon(xc, mask)
            @test JosephsonCircuits.real_to_complex!(fill(ComplexF64(NaN, NaN), d), xr, mask) == JosephsonCircuits.real_to_complex(xr, mask)
            # conj_input flips the stored imaginary parts, scale hits the real modes
            @test JosephsonCircuits.complex_to_real(xc, mask; conj_input = true) == JosephsonCircuits.complex_to_real(conj(xc), mask)
            s = 0.5
            xs = [disr(mask, i) ? s*xc[i] : xc[i] for i in 1:d]
            @test JosephsonCircuits.complex_to_real(xc, mask; realscale = s) == JosephsonCircuits.complex_to_real(xs, mask)
            @test JosephsonCircuits.real_to_complex(JosephsonCircuits.complex_to_real(xc, mask), mask; realscale = s) == dcanon(xs, mask)
            @test_throws DimensionMismatch JosephsonCircuits.complex_to_real!(Vector{Float64}(undef, length(xr)+1), xc, mask)
        end
    end

    @testset "dense: complex_to_real" begin
        for (rm, cm, mmul, nmul) in DCASES
            m, n = length(rm)*mmul, length(cm)*nmul
            A = rand(ComplexF64, m, n)
            xc = dcanon(rand(ComplexF64, n), cm)
            for cj in (false, true)
                Ar = JosephsonCircuits.complex_to_real(A, rm, cm; conj_input = cj)
                @test size(Ar) == (JosephsonCircuits.realdim(m, rm), JosephsonCircuits.realdim(n, cm))
                @test Ar == dref(A, rm, cm; conj_input = cj)
                @test Ar * JosephsonCircuits.complex_to_real(xc, cm) ≈ JosephsonCircuits.complex_to_real(A * (cj ? conj(xc) : xc), rm)
                S = fill(NaN, size(Ar))
                @test JosephsonCircuits.complex_to_real!(S, A, rm, cm; conj_input = cj) === S
                @test S == Ar
                Ac = JosephsonCircuits.real_to_complex(Ar, rm, cm; conj_input = cj)
                @test Ac == ddropq(A, rm, cm)
                @test JosephsonCircuits.complex_to_real(Ac, rm, cm; conj_input = cj) == Ar
                C = fill(ComplexF64(NaN, NaN), m, n)
                @test JosephsonCircuits.real_to_complex!(C, Ar, rm, cm; conj_input = cj) == Ac
            end
            @test JosephsonCircuits.complex_to_real(Float32.(real(A)) .+ 0im .|> ComplexF32, rm, cm) isa Matrix{Float32}
        end
    end

    @testset "dense: realscale" begin
        for (rm, cm, mmul, nmul) in DCASES
            m, n = length(rm)*mmul, length(cm)*nmul
            A = rand(ComplexF64, m, n)
            base = JosephsonCircuits.complex_to_real(A, rm, cm)
            @test JosephsonCircuits.complex_to_real(A, rm, cm; realrowscale = 1, realcolscale = 1) == base
            s, xc = 0.5, dcanon(rand(ComplexF64, n), cm)
            xs = [disr(cm, i) ? s*xc[i] : xc[i] for i in 1:n]
            @test JosephsonCircuits.complex_to_real(A, rm, cm; realcolscale = s) * JosephsonCircuits.complex_to_real(xc, cm) ≈ JosephsonCircuits.complex_to_real(A * xs, rm)
            b  = A * xc
            bs = [disr(rm, i) ? s*b[i] : b[i] for i in 1:m]
            @test JosephsonCircuits.complex_to_real(A, rm, cm; realrowscale = s) * JosephsonCircuits.complex_to_real(xc, cm) ≈ JosephsonCircuits.complex_to_real(bs, rm)
            for (rr, cc) in ((2.0, 1.0), (1.0, 3.0), (2.0, 3.0), (0.0, 1.0), (1.0, 0.0))
                @test JosephsonCircuits.complex_to_real(A, rm, cm; realrowscale = rr, realcolscale = cc) ==
                      dref(A, rm, cm; rs = rr, cs = cc)
            end
            @test JosephsonCircuits.complex_to_real(A, rm, cm; conj_input = true, realcolscale = s) ==
                  dref(A, rm, cm; conj_input = true, cs = s)
        end
    end


    # @testset "dense: allocation-free in-place" begin
    #     rm = [true,false,false,false,false]
    #     m, n = 200, 200
    #     A = rand(ComplexF64, m, n); B = rand(ComplexF64, m, n)
    #     xc = rand(ComplexF64, n); xr = JosephsonCircuits.complex_to_real(xc, rm)
    #     JosephsonCircuits.complex_to_real!(Ar, A, rm, rm); JosephsonCircuits.complex_to_real!(Ar, A, rm, rm; conj_input = true, realcolscale = 0.5)
    #     JosephsonCircuits.complex_to_real!(xr, xc, rm); JosephsonCircuits.real_to_complex!(xc, xr, rm)
    #     @test (@allocated JosephsonCircuits.complex_to_real!(Ar, A, rm, rm)) == 0
    #     @test (@allocated JosephsonCircuits.complex_to_real!(Ar, A, rm, rm; conj_input = true, realcolscale = 0.5)) == 0
    #     @test (@allocated JosephsonCircuits.real_to_complex!(C, Ar, rm, rm)) == 0
    #     @test (@allocated JosephsonCircuits.complex_to_real!(xr, xc, rm)) == 0
    #     @test (@allocated JosephsonCircuits.real_to_complex!(xc, xr, rm)) == 0
    # end

    @testset "dense: shape checks" begin
        rm = [true,false,false,false,false]
        m, n = 20, 15
        A = rand(ComplexF64, m, n)
        @test_throws DimensionMismatch JosephsonCircuits.complex_to_real(A, rm, [true,false])          # 15 % 2 != 0
        @test_throws DimensionMismatch JosephsonCircuits.complex_to_real!(zeros(JosephsonCircuits.realdim(m,rm)+1, JosephsonCircuits.realdim(n,rm)), A, rm, rm)
        @test_throws DimensionMismatch JosephsonCircuits.real_to_complex!(A, zeros(JosephsonCircuits.realdim(m,rm), JosephsonCircuits.realdim(n,rm)+3), rm, rm)
        # views work
        Abig = rand(ComplexF64, m+4, n+4)
        @test JosephsonCircuits.complex_to_real(view(Abig, 1:m, 1:n), rm, rm) == JosephsonCircuits.complex_to_real(Abig[1:m, 1:n], rm, rm)
    end

    @testset "real_to_complex rejects an incomplete row mode" begin
        # `complex_to_real` always stores both real slots of a complex row
        # mode. If a caller passes a pattern where the second slot is stored
        # but the first is not, the imaginary part has no entry to be folded
        # into: the write used to land on an unrelated earlier entry, or on
        # index 0 of a zero-length nzval, which `@inbounds` did not catch.
        mask = [false, false]
        Ar = JosephsonCircuits.SparseArrays.sparse([2], [1], [1.0], 4, 4)
        @test_throws ArgumentError JosephsonCircuits.real_to_complex(Ar, mask)

        # the second slot of the second row mode, with the first absent
        Ar2 = JosephsonCircuits.SparseArrays.sparse([1, 4], [1, 1], [1.0, 2.0], 4, 4)
        @test_throws ArgumentError JosephsonCircuits.real_to_complex(Ar2, mask)

        # a well formed pattern still round trips
        A = JosephsonCircuits.SparseArrays.sparse([1, 2], [1, 2], ComplexF64[1.0+2.0im, 3.0-1.0im], 2, 2)
        Arok = JosephsonCircuits.complex_to_real(A, mask)
        @test JosephsonCircuits.real_to_complex(Arok, mask) == A
    end

end
