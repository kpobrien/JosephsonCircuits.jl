using JosephsonCircuits
using LinearAlgebra
using Test

@testset verbose = true "quantumoptics" begin

    @testset "symplectic form" begin

        # Serafini B.2, the symplectic form is a member of the symplectic
        # group
        @test JosephsonCircuits.is_symplectic_block(
            JosephsonCircuits.symplectic_form_block(4),
        )

        @test JosephsonCircuits.is_symplectic_pair(
            JosephsonCircuits.symplectic_form_pair(4),
        )

        # Serafini B.3, the inverse equals the adjoint
        @test isapprox(
            inv(Matrix(JosephsonCircuits.symplectic_form_block(4))),
            adjoint(JosephsonCircuits.symplectic_form_block(4)),
        )

        @test isapprox(
            inv(Matrix(JosephsonCircuits.symplectic_form_pair(4))),
            adjoint(JosephsonCircuits.symplectic_form_pair(4)),
        )

        # Serafini B.3, the adjoint equals the negative of the symplectic
        # form

        @test isapprox(
            adjoint(JosephsonCircuits.symplectic_form_block(4)),
            -JosephsonCircuits.symplectic_form_block(4),
        )

        @test isapprox(
            adjoint(JosephsonCircuits.symplectic_form_pair(4)),
            -JosephsonCircuits.symplectic_form_pair(4),
        )

        # Serafini pg. 31, symplectic form times itself is minus identity
        @test isapprox(
            JosephsonCircuits.symplectic_form_block(4)^2,
            -I(2 * 4),
        )

        @test isapprox(
            JosephsonCircuits.symplectic_form_pair(4)^2,
            -I(2 * 4),
        )

        # Serafini pg. 31, symplectic form times adjoint is the identity
        @test isapprox(
            JosephsonCircuits.symplectic_form_block(4) * JosephsonCircuits.symplectic_form_block(4)',
            I(2 * 4),
        )

        @test isapprox(
            JosephsonCircuits.symplectic_form_pair(4) * JosephsonCircuits.symplectic_form_pair(4)',
            I(2 * 4),
        )

        # convert between the symplectic forms
        @test isapprox(
            JosephsonCircuits.symplectic_form_block(4),
            JosephsonCircuits.pair_to_block(JosephsonCircuits.symplectic_form_pair(4)),
        )

        @test isapprox(
            JosephsonCircuits.symplectic_form_block(4),
            JosephsonCircuits.pair_to_block2(JosephsonCircuits.symplectic_form_pair(4)),
        )

        @test isapprox(
            JosephsonCircuits.symplectic_form_pair(4),
            JosephsonCircuits.block_to_pair(JosephsonCircuits.symplectic_form_block(4)),
        )

        @test isapprox(
            JosephsonCircuits.symplectic_form_pair(4),
            JosephsonCircuits.block_to_pair2(JosephsonCircuits.symplectic_form_block(4)),
        )

        # the permutation matrices reorder a matrix which is not square as
        # the permutations do
        M = rand(6, 4)
        @test JosephsonCircuits.block_to_pair2(M) == JosephsonCircuits.block_to_pair(M)
        @test JosephsonCircuits.pair_to_block2(M) == JosephsonCircuits.pair_to_block(M)

        # test the conversions and their inverses
        @test isapprox(
            JosephsonCircuits.symplectic_form_pair(4),
            JosephsonCircuits.block_to_pair(JosephsonCircuits.pair_to_block(JosephsonCircuits.symplectic_form_pair(4))),
        )

        @test isapprox(
            JosephsonCircuits.symplectic_form_block(4),
            JosephsonCircuits.pair_to_block(JosephsonCircuits.block_to_pair(JosephsonCircuits.symplectic_form_block(4))),
        )

    end

    @testset "indefinite hermitian form" begin

        # imaginary number times indefinite Hermitian form is a member of the
        # bogoliubov group
        @test JosephsonCircuits.is_bogoliubov_block(
            -im*JosephsonCircuits.indefinite_hermitian_form_block(4),
        )

        @test JosephsonCircuits.is_bogoliubov_pair(
            -im*JosephsonCircuits.indefinite_hermitian_form_pair(4),
        )

        # convert symplectic form to indefinite hermitian form
        @test isapprox(
            JosephsonCircuits.quadrature_to_ladder_block(JosephsonCircuits.symplectic_form_block(4)),
            -im*JosephsonCircuits.indefinite_hermitian_form_block(4)
        )

        @test isapprox(
            JosephsonCircuits.quadrature_to_ladder_pair(JosephsonCircuits.symplectic_form_pair(4)),
            -im*JosephsonCircuits.indefinite_hermitian_form_pair(4)
        )

        # convert indefinite hermitian form to symplectic form
        @test isapprox(
            JosephsonCircuits.symplectic_form_block(4),
            JosephsonCircuits.ladder_to_quadrature_block(-im*JosephsonCircuits.indefinite_hermitian_form_block(4))
        )

        @test isapprox(
            JosephsonCircuits.symplectic_form_pair(4),
            JosephsonCircuits.ladder_to_quadrature_pair(-im*JosephsonCircuits.indefinite_hermitian_form_pair(4))
        )

    end

    @testset "is_positive_semi_definite" begin

        @test JosephsonCircuits.is_positive_semi_definite([1 0;0 1])
        @test JosephsonCircuits.is_positive_semi_definite([1 0;0 0])
        @test !JosephsonCircuits.is_positive_semi_definite([1 0;0 -1])
        @test !JosephsonCircuits.is_positive_semi_definite([1 1;0 1])

    end

    @testset "is_positive_definite" begin

        @test JosephsonCircuits.is_positive_definite([1 0;0 1])
        @test !JosephsonCircuits.is_positive_definite([1 0;0 0])
        @test !JosephsonCircuits.is_positive_definite([1 0;0 -1])
        @test !JosephsonCircuits.is_positive_definite([1 1;0 1])

    end

    @testset "is_symplectic_pair and is_symplectic_block" begin

        @test_throws(
            ErrorException,
            JosephsonCircuits.is_symplectic_pair(rand(9,9)),
        )
        @test_throws(
            ErrorException,
            JosephsonCircuits.is_symplectic_block(rand(9,9)),
        )

    end

    @testset "is_cptp" begin

        @test !JosephsonCircuits.is_cptp_quadrature_pair([1 0;0 1],[1 0;1 1])
        @test JosephsonCircuits.is_cptp_quadrature_pair([1 0;0 1],[0 0;0 0])

        # a Gaussian unitary adds no noise, so it is CPTP with Y = 0, where
        # Y + im*(Ω - X*Ω*X') is zero up to rounding of either sign
        n = 2
        Z = zeros(2n, 2n)
        @test all(JosephsonCircuits.is_cptp_quadrature_pair(
            JosephsonCircuits.rand_symplectic_pair(n), Z) for trial in 1:20)
        @test all(JosephsonCircuits.is_cptp_quadrature_block(
            JosephsonCircuits.rand_symplectic_block(n), Z) for trial in 1:20)
        @test all(JosephsonCircuits.is_cptp_ladder_pair(
            JosephsonCircuits.rand_bogoliubov_pair(n), complex(Z)) for trial in 1:20)
        @test all(JosephsonCircuits.is_cptp_ladder_block(
            JosephsonCircuits.rand_bogoliubov_block(n), complex(Z)) for trial in 1:20)

        # a phase-insensitive amplifier of gain G, X = sqrt(G)*I, is CPTP
        # when it adds at least the noise Y = (G-1)*I of the Caves limit,
        # the same in the quadrature and the ladder bases
        G = 3.0
        X = sqrt(G) * Matrix(1.0I, 2n, 2n)
        for is_cptp_form in (JosephsonCircuits.is_cptp_quadrature_pair,
                JosephsonCircuits.is_cptp_quadrature_block,
                JosephsonCircuits.is_cptp_ladder_pair,
                JosephsonCircuits.is_cptp_ladder_block)
            @test is_cptp_form(X, (G - 1) * Matrix(1.0I, 2n, 2n))
            @test !is_cptp_form(X, (G - 1) / 2 * Matrix(1.0I, 2n, 2n))
        end

        # the tolerances can be given: a violation of 1e-6 is refused at an
        # absolute tolerance below it and accepted at one above it
        Y = (G - 1 - 1e-6) * Matrix(1.0I, 2n, 2n)
        @test !JosephsonCircuits.is_cptp_quadrature_pair(X, Y; atol = 1e-8)
        @test JosephsonCircuits.is_cptp_quadrature_pair(X, Y; atol = 1e-4)

    end

    @testset "random positive definite" begin

        @test JosephsonCircuits.is_positive_definite(JosephsonCircuits.rand_positive_definite(4))

    end

    @testset "random unitary" begin

        @test JosephsonCircuits.is_unitary(JosephsonCircuits.rand_unitary(4))

    end

    @testset "random orthogonal" begin

        @test JosephsonCircuits.is_orthogonal(JosephsonCircuits.rand_orthogonal(4))

    end

    @testset "random positive semi-definite" begin

        @test JosephsonCircuits.is_positive_semi_definite(JosephsonCircuits.rand_positive_semi_definite(2,2))

    end

    @testset "random matrices of each group" begin

        # each generator gives a member of its group, in either order, for
        # the default, real and complex element types
        JC = JosephsonCircuits
        for (rand_pair, rand_block, is_pair, is_block) in (
                (JC.rand_symplectic_pair, JC.rand_symplectic_block,
                    JC.is_symplectic_pair, JC.is_symplectic_block),
                (JC.rand_orthogonal_symplectic_pair, JC.rand_orthogonal_symplectic_block,
                    JC.is_orthogonal_symplectic_pair, JC.is_orthogonal_symplectic_block),
                (JC.rand_positive_definite_symplectic_pair, JC.rand_positive_definite_symplectic_block,
                    JC.is_positive_definite_symplectic_pair, JC.is_positive_definite_symplectic_block),
                (JC.rand_conjugate_symplectic_pair, JC.rand_conjugate_symplectic_block,
                    JC.is_conjugate_symplectic_pair, JC.is_conjugate_symplectic_block),
                (JC.rand_bogoliubov_pair, JC.rand_bogoliubov_block,
                    JC.is_bogoliubov_pair, JC.is_bogoliubov_block),
                (JC.rand_orthogonal_bogoliubov_pair, JC.rand_orthogonal_bogoliubov_block,
                    JC.is_orthogonal_bogoliubov_pair, JC.is_orthogonal_bogoliubov_block),
                (JC.rand_pseudo_unitary_pair, JC.rand_pseudo_unitary_block,
                    JC.is_pseudo_unitary_pair, JC.is_pseudo_unitary_block))
            for args in ((4,), (Float64, 4), (Complex{Float64}, 4))
                @test is_pair(rand_pair(args...))
                @test is_block(rand_block(args...))
            end
        end

    end

    @testset "random cptp quadrature" begin

        # symplectic
        @test JosephsonCircuits.is_cptp_quadrature_pair(JosephsonCircuits.rand_cptp_quadrature_pair(4)...)
        @test JosephsonCircuits.is_cptp_quadrature_pair(JosephsonCircuits.rand_cptp_quadrature_pair(Float64,4)...)
        
        @test JosephsonCircuits.is_cptp_quadrature_block(JosephsonCircuits.rand_cptp_quadrature_block(4)...)
        @test JosephsonCircuits.is_cptp_quadrature_block(JosephsonCircuits.rand_cptp_quadrature_block(Float64,4)...)

        # an environment of other than as many modes as the system: the
        # noise it adds has at most its rank
        for nenv in (1, 3)
            X, Y = JosephsonCircuits.rand_cptp_quadrature_pair(Float64, 2; nenv = nenv)
            @test JosephsonCircuits.is_cptp_quadrature_pair(X, Y) && rank(Y) <= 2*nenv
            X, Y = JosephsonCircuits.rand_cptp_quadrature_pair(2; nenv = nenv)
            @test JosephsonCircuits.is_cptp_quadrature_pair(X, Y) && rank(Y) <= 2*nenv
            X, Y = JosephsonCircuits.rand_cptp_quadrature_block(2; nenv = nenv)
            @test JosephsonCircuits.is_cptp_quadrature_block(X, Y) && rank(Y) <= 2*nenv
        end

    end

    @testset "random cptp ladder" begin

        # bogoliubov
        @test JosephsonCircuits.is_cptp_ladder_pair(JosephsonCircuits.rand_cptp_ladder_pair(4)...)
        @test JosephsonCircuits.is_cptp_ladder_pair(JosephsonCircuits.rand_cptp_ladder_pair(Complex{Float64},4)...)
        
        @test JosephsonCircuits.is_cptp_ladder_block(JosephsonCircuits.rand_cptp_ladder_block(4)...)
        @test JosephsonCircuits.is_cptp_ladder_block(JosephsonCircuits.rand_cptp_ladder_block(Complex{Float64},4)...)

        for nenv in (1, 3)
            X, Y = JosephsonCircuits.rand_cptp_ladder_pair(Complex{Float64}, 2; nenv = nenv)
            @test JosephsonCircuits.is_cptp_ladder_pair(X, Y) && rank(Y) <= 2*nenv
            X, Y = JosephsonCircuits.rand_cptp_ladder_pair(2; nenv = nenv)
            @test JosephsonCircuits.is_cptp_ladder_pair(X, Y) && rank(Y) <= 2*nenv
            X, Y = JosephsonCircuits.rand_cptp_ladder_block(2; nenv = nenv)
            @test JosephsonCircuits.is_cptp_ladder_block(X, Y) && rank(Y) <= 2*nenv
        end

    end

    @testset "scattering to ladder and quadrature" begin

        # matrix functions
        S = randn(Complex{Float64},8,8)
        for w in [sign.(randn(4)), sign.(randn(8))]
            @test isapprox(
                JosephsonCircuits.ladder_to_quadrature_pair(JosephsonCircuits.scattering_to_ladder_pair(S,w)),
                JosephsonCircuits.scattering_to_quadrature_pair(S,w),
            )

            @test isapprox(
                JosephsonCircuits.ladder_to_quadrature_block(JosephsonCircuits.scattering_to_ladder_block(S,w)),
                JosephsonCircuits.scattering_to_quadrature_block(S,w),
            )

            @test isapprox(
                JosephsonCircuits.ladder_to_quadrature_pair(JosephsonCircuits.scattering_to_ladder_pair(S,w)),
                JosephsonCircuits.block_to_pair(JosephsonCircuits.scattering_to_quadrature_block(S,w)),
            )
        end

        # vector functions
        @test_throws(
            DimensionMismatch,
            JosephsonCircuits.scattering_to_ladder_pair!(zeros(10),zeros(10),ones(5)),
        )

        @test_throws(
            DimensionMismatch,
            JosephsonCircuits.scattering_to_ladder_block!(zeros(10),zeros(10),ones(5)),
        )

        @test_throws(
            DimensionMismatch,
            JosephsonCircuits.scattering_to_quadrature_pair!(zeros(10),zeros(10),ones(5)),
        )

        @test_throws(
            DimensionMismatch,
            JosephsonCircuits.scattering_to_quadrature_block!(zeros(10),zeros(10),ones(5)),
        )

        @test_throws(
            DimensionMismatch,
            JosephsonCircuits.scattering_to_ladder_pair!(zeros(10),zeros(5),ones(6)),
        )

        @test_throws(
            DimensionMismatch,
            JosephsonCircuits.scattering_to_ladder_block!(zeros(10),zeros(5),ones(6)),
        )

        @test_throws(
            DimensionMismatch,
            JosephsonCircuits.scattering_to_quadrature_pair!(zeros(10),zeros(5),ones(6)),
        )

        @test_throws(
            DimensionMismatch,
            JosephsonCircuits.scattering_to_quadrature_block!(zeros(10),zeros(5),ones(6)),
        )

        # the block form of a non-square matrix holds the entries of its pair
        # form, reordered, and converts back
        for (n, m) in [(2, 1), (1, 2), (2, 3), (3, 2)]
            S = randn(Complex{Float64}, n, m)
            for w in ([1.0, -1.0, 1.0], [-1.0, 1.0, -1.0])
                @test JosephsonCircuits.scattering_to_ladder_block(S, w) ==
                    JosephsonCircuits.pair_to_block(JosephsonCircuits.scattering_to_ladder_pair(S, w))
                @test JosephsonCircuits.scattering_to_quadrature_block(S, w) ==
                    JosephsonCircuits.pair_to_block(JosephsonCircuits.scattering_to_quadrature_pair(S, w))
                @test isapprox(JosephsonCircuits.ladder_to_scattering_block(
                    JosephsonCircuits.scattering_to_ladder_block(S, w), w), S)
                @test isapprox(JosephsonCircuits.quadrature_to_scattering_block(
                    JosephsonCircuits.scattering_to_quadrature_block(S, w), w), S)
                # the quadrature form is the ladder form in the quadrature
                # basis of its rows and of its columns, in either order
                @test isapprox(JosephsonCircuits.ladder_to_quadrature_pair(
                    JosephsonCircuits.scattering_to_ladder_pair(S, w)),
                    JosephsonCircuits.scattering_to_quadrature_pair(S, w))
                @test isapprox(JosephsonCircuits.ladder_to_quadrature_block(
                    JosephsonCircuits.scattering_to_ladder_block(S, w)),
                    JosephsonCircuits.scattering_to_quadrature_block(S, w))
                @test isapprox(JosephsonCircuits.quadrature_to_ladder_pair(
                    JosephsonCircuits.scattering_to_quadrature_pair(S, w)),
                    JosephsonCircuits.scattering_to_ladder_pair(S, w))
                @test isapprox(JosephsonCircuits.quadrature_to_ladder_block(
                    JosephsonCircuits.scattering_to_quadrature_block(S, w)),
                    JosephsonCircuits.scattering_to_ladder_block(S, w))
            end
        end

        # the sign convention, literally: a mode of positive frequency is
        # an annihilation operator, whose amplitude and its conjugate take
        # the first and the second operator of its pair, and a mode of
        # negative frequency a creation operator, which swaps them; the
        # imaginary part of the amplitude of a creation operator is that of
        # a conjugate
        a, b, c, d = 0.3 + 0.4im, 0.1 - 0.7im, -0.2 + 0.5im, 0.6 + 0.1im
        @test JosephsonCircuits.scattering_to_ladder_pair([a b; c d], [1.0, -1.0]) ==
            [a 0 0 b; 0 conj(a) conj(b) 0; 0 conj(c) conj(d) 0; c 0 0 d]
        @test JosephsonCircuits.scattering_to_quadrature_pair(fill(a, 1, 1), [1.0]) ==
            [real(a) -imag(a); imag(a) real(a)]
        @test JosephsonCircuits.scattering_to_quadrature_pair(fill(a, 1, 1), [-1.0]) ==
            [real(a) imag(a); -imag(a) real(a)]

        # each vector conversion is its matrix conversion acting on the
        # vector: f(S*v, w) == f(S, w)*f(v, w)
        S = randn(Complex{Float64}, 6, 6)
        v = randn(Complex{Float64}, 6)
        w = [1.0, -1.0, 1.0]
        for f in (JosephsonCircuits.scattering_to_ladder_pair,
                JosephsonCircuits.scattering_to_ladder_block,
                JosephsonCircuits.scattering_to_quadrature_pair,
                JosephsonCircuits.scattering_to_quadrature_block)
            @test isapprox(f(S * v, w), f(S, w) * f(v, w))
        end
    end

    @testset "negative examples" begin
        # a two-mode squeezer between a mode of positive and one of negative
        # frequency, a signal and its idler, is a Bogoliubov transformation
        # and not a unitary one; between two modes of positive frequency the
        # same matrix would amplify without an idler, and is neither
        # Bogoliubov, pseudo-unitary nor symplectic
        r = 0.5
        S2 = [cosh(r) sinh(r); sinh(r) cosh(r)]
        L = JosephsonCircuits.scattering_to_ladder_pair(S2, [1.0, -1.0])
        @test JosephsonCircuits.is_bogoliubov_pair(L)
        @test !JosephsonCircuits.is_unitary(L)
        @test JosephsonCircuits.is_bogoliubov_block(
            JosephsonCircuits.scattering_to_ladder_block(S2, [1.0, -1.0]))
        w = [1.0, 1.0]
        L = JosephsonCircuits.scattering_to_ladder_pair(S2, w)
        @test !JosephsonCircuits.is_bogoliubov_pair(L)
        @test !JosephsonCircuits.is_pseudo_unitary_pair(L)
        @test !JosephsonCircuits.is_bogoliubov_block(
            JosephsonCircuits.scattering_to_ladder_block(S2, w))
        @test !JosephsonCircuits.is_symplectic_pair(
            JosephsonCircuits.scattering_to_quadrature_pair(S2, w))
        @test !JosephsonCircuits.is_symplectic_block(
            JosephsonCircuits.scattering_to_quadrature_block(S2, w))
        # a shear is neither unitary nor orthogonal
        @test !JosephsonCircuits.is_unitary([1.0 1.0; 0.0 1.0])
        @test !JosephsonCircuits.is_orthogonal([1.0 1.0; 0.0 1.0])
    end

    @testset "sizes the conversions refuse" begin
        # a form of odd size, and no mode frequencies
        @test_throws DimensionMismatch JosephsonCircuits.ladder_to_scattering_pair(rand(5, 5), [1.0])
        @test_throws DimensionMismatch JosephsonCircuits.quadrature_to_scattering_block(rand(5, 5), [1.0])
        @test_throws ArgumentError JosephsonCircuits.scattering_to_ladder_pair(rand(2, 2), Float64[])
    end

    @testset "conversions back to scattering parameters" begin

        JC = JosephsonCircuits
        for (to, back, rand_form) in (
                (JC.scattering_to_ladder_pair, JC.ladder_to_scattering_pair, JC.rand_symplectic_pair),
                (JC.scattering_to_ladder_block, JC.ladder_to_scattering_block, JC.rand_symplectic_block),
                (JC.scattering_to_quadrature_pair, JC.quadrature_to_scattering_pair, JC.rand_symplectic_pair),
                (JC.scattering_to_quadrature_block, JC.quadrature_to_scattering_block, JC.rand_symplectic_block))

            # complex and real floating point input, and the form stored
            # complex
            for T in (Complex{Float64}, Float64)
                X = JC.rand_unitary(T, 10)
                w = sign.(randn(size(X, 1)))
                S = to(X, w)
                @test isapprox(back(S, w), X)
                @test isapprox(back(complex(S), w), X)
            end

            # a matrix which is not a form is refused
            @test_throws ErrorException back(rand_form(5), sign.(randn(10)))

            # complex floating point input carrying rounding is converted
            # within the default tolerance, and the tolerances can be given
            X = JC.rand_unitary(6)
            w = [1.0, -1.0, 1.0]
            U = JC.rand_unitary(12)
            S = to(X, w)*(U*U')
            @test isapprox(back(S, w), X)
            @test isapprox(back(S, w; rtol = 1e-10), X)
        end

        # ladder forms which break one of the relations between their
        # entries: a conjugate of another sign, an entry of the wrong
        # operator, in either order
        for (back, forms) in (
                (JC.ladder_to_scattering_pair, ([1 0 0 0;0 -1 0 0;0 0 1 0;0 0 0 1],
                    [1 1 0 0;0 1 0 0;0 0 1 0;0 0 0 1], [1 0 0 0;1 1 0 0;0 0 1 0;0 0 0 1])),
                (JC.ladder_to_scattering_block, ([1 0 0 0;0 -1 0 0;0 0 1 0;0 0 0 1],
                    [1 0 0 0;0 1 0 1;0 0 1 0;0 0 0 1], [1 0 0 0;0 1 0 0;0 0 1 0;1 0 0 1])))
            for form in forms
                @test_throws ErrorException back(form, [1, 1, 1, 1])
            end
        end

    end

    @testset "port mode conversion" begin

        for (S, Nmodes) in [(rand(Complex{Float64},8,8), 2),
                (rand(Complex{Float64},8,12), 2), (rand(Complex{Float64},12,12), 3),
                (rand(Complex{Float64},12,18), 3), (rand(Complex{Float64},12,12), 2)]
            @test isapprox(JosephsonCircuits.ports_modes_to_modes_ports_scattering(JosephsonCircuits.modes_ports_to_ports_modes_scattering(S,Nmodes),Nmodes),S)
            @test isapprox(JosephsonCircuits.ports_modes_to_modes_ports_pair(JosephsonCircuits.modes_ports_to_ports_modes_pair(S,Nmodes),Nmodes),S)
            @test isapprox(JosephsonCircuits.ports_modes_to_modes_ports_block(JosephsonCircuits.modes_ports_to_ports_modes_block(S,Nmodes),Nmodes),S)
        end

        # permuting the pair or block form of a scattering matrix is the pair
        # or block form of the permuted scattering matrix, for every number
        # of ports and modes; with one port or one mode the two orderings
        # are the same
        for (Nports, Nmodes) in [(2, 3), (3, 2), (1, 3), (3, 1), (1, 1)]
            N = Nports*Nmodes
            S = rand(Complex{Float64}, N, N)
            s = [isodd(i) ? 1.0 : -1.0 for i in 1:N]
            for (perm, fpair, fblock, fscattering) in [
                    (JosephsonCircuits.ports_modes_to_modes_ports_perm,
                        JosephsonCircuits.ports_modes_to_modes_ports_pair,
                        JosephsonCircuits.ports_modes_to_modes_ports_block,
                        JosephsonCircuits.ports_modes_to_modes_ports_scattering),
                    (JosephsonCircuits.modes_ports_to_ports_modes_perm,
                        JosephsonCircuits.modes_ports_to_ports_modes_pair,
                        JosephsonCircuits.modes_ports_to_ports_modes_block,
                        JosephsonCircuits.modes_ports_to_ports_modes_scattering)]
                p = perm(Nports, Nmodes)
                @test isperm(p)
                @test fscattering(S, Nmodes) == S[p, p]
                @test fpair(JosephsonCircuits.scattering_to_ladder_pair(S, s), Nmodes) ==
                    JosephsonCircuits.scattering_to_ladder_pair(S[p, p], s[p])
                @test fblock(JosephsonCircuits.scattering_to_ladder_block(S, s), Nmodes) ==
                    JosephsonCircuits.scattering_to_ladder_block(S[p, p], s[p])
                if Nports == 1 || Nmodes == 1
                    @test fscattering(S, Nmodes) == S
                end
            end
        end

        # an axis which does not hold whole ports
        @test_throws(
            DimensionMismatch,
            JosephsonCircuits.ports_modes_to_modes_ports_scattering(rand(6, 6), 4),
        )
        @test_throws(
            DimensionMismatch,
            JosephsonCircuits.modes_ports_to_ports_modes_block(rand(6, 6), 2),
        )
        @test_throws(
            DimensionMismatch,
            JosephsonCircuits.ports_modes_to_modes_ports_pair(rand(5, 5), 1),
        )

    end

    @testset "polar" begin

        A = randn(Complex{Float64}, 4, 4)
        P, Y = JosephsonCircuits.polar(A)
        @test isapprox(P * Y, A)
        @test isapprox(Y * Y', I(4))

    end

    @testset "williamson pair" begin

        # the values are the symplectic eigenvalues, the moduli of the
        # eigenvalues of im*Ω*M, each twice
        Omega = JosephsonCircuits.symplectic_form_pair(2)
        for M in [JosephsonCircuits.rand_positive_definite(4), JosephsonCircuits.rand_positive_semi_definite(2,2)]
            d, S = JosephsonCircuits.williamson_pair(M)
            @test JosephsonCircuits.is_symplectic_pair(S)
            # literature convention
            @test isapprox(S*Diagonal(d)*transpose(S),M)
            @test isapprox(sort(d), sort(abs.(eigvals(im * Omega * M))))
        end

        # a covariance S*Diagonal(d)*transpose(S), which rounding leaves
        # symmetric only approximately, here by one ulp
        S0 = JosephsonCircuits.rand_symplectic_pair(2)
        d0 = [1.5, 1.5, 3.0, 3.0]
        M = S0 * Diagonal(d0) * transpose(S0)
        M[1, 2] = nextfloat(M[2, 1])
        d, S = JosephsonCircuits.williamson_pair(M)
        @test isapprox(S*Diagonal(d)*transpose(S), M)
        @test isapprox(sort(d), d0)
        d, S = JosephsonCircuits.williamson_block(JosephsonCircuits.pair_to_block(M))
        @test isapprox(S*Diagonal(d)*transpose(S), JosephsonCircuits.pair_to_block(M))

        # a matrix of rank two up to rounding, decomposed in either order
        M = [
            1.166893242623039673e-01  -1.048482006204835837e-02  -3.493446128036903353e-02   1.565614325188508238e-01;
           -1.048482006204835837e-02   4.805855935380488053e-01  -1.069468907857842987e+00  -3.770319099460901491e-01;
           -3.493446128036903353e-02  -1.069468907857842987e+00   2.409089316529508640e+00   7.648117805967328264e-01;
            1.565614325188508238e-01  -3.770319099460901491e-01   7.648117805967328264e-01   4.847266530147685271e-01;
        ]
        for (williamson, is_symplectic, Omega) in (
                (JosephsonCircuits.williamson_pair, JosephsonCircuits.is_symplectic_pair,
                    JosephsonCircuits.symplectic_form_pair(2)),
                (JosephsonCircuits.williamson_block, JosephsonCircuits.is_symplectic_block,
                    JosephsonCircuits.symplectic_form_block(2)))
            d, S = williamson(M)
            @test is_symplectic(S)
            @test isapprox(S*Diagonal(d)*transpose(S), M)
            @test isapprox(sort(d), sort(abs.(eigvals(im * Omega * M))); atol = 1e-12)
        end

        # a positive semi-definite matrix of odd rank, and one whose range
        # holds the positions alone, have no symplectic normal form
        @test_throws ArgumentError JosephsonCircuits.williamson_pair([1.0 0; 0 0])
        @test_throws ArgumentError JosephsonCircuits.williamson_pair(
            Matrix(Diagonal([1.0, 0, 1, 1])))
        @test_throws ArgumentError JosephsonCircuits.williamson_pair(
            Matrix(Diagonal([1.0, 0, 1, 0])))
    end


    @testset "williamson block " begin

        Omega = JosephsonCircuits.symplectic_form_block(2)
        for M in [JosephsonCircuits.rand_positive_definite(4), JosephsonCircuits.rand_positive_semi_definite(2,2)]
            d, S = JosephsonCircuits.williamson_block(M)
            @test JosephsonCircuits.is_symplectic_block(S)
            # literature convention
            @test isapprox(S*Diagonal(d)*transpose(S),M)
            @test isapprox(sort(d), sort(abs.(eigvals(im * Omega * M))))
        end
    end

    @testset "cholesky_williamson" begin

        L1, rankL1 = JosephsonCircuits.cholesky_williamson([1 0;0 1])
        L2, rankL2 = JosephsonCircuits.cholesky_williamson(JosephsonCircuits.SparseArrays.sparse([1 0;0 1]))

        @test isapprox(L1,L2)
        @test isapprox(rankL1,rankL2)

        @test_throws(
            ErrorException,
            JosephsonCircuits.cholesky_williamson([1 0 0 0;0 1 0 0;0 0 0 0;0 0 0 -1]),
        )

        @test_throws(
            ErrorException,
            JosephsonCircuits.cholesky_williamson(JosephsonCircuits.SparseArrays.sparse([1 0 0 0;0 1 0 0;0 0 0 0;0 0 0 -1])),
        )

        # positive semi-definite, of odd rank
        @test_throws ArgumentError JosephsonCircuits.cholesky_williamson([1 0;0 0])
    end

    @testset "autonne_takagi complex" begin

        A = Symmetric(rand(Complex{Float64}, 4, 4))
        values, vectors = JosephsonCircuits.autonne_takagi(A)
        @test isapprox(vectors * Diagonal(values) * transpose(vectors), A)
        @test isapprox(vectors * vectors', I(size(A, 1)))

        # https://github.com/XanaduAI/thewalrus/pull/403
        # https://gist.github.com/tomdodd4598/f1b42a1c491c43c7661b90685160496b/revisions
        A = exp(im * 0.0) * [-1.3197074035840624+3.2134893495780524e-16im -0.05059154327551117+1.6739618537092438e-16im 0.21057507448953267-3.1495068151629784e-16im -0.2805588371720386+2.852187447161915e-14im; -0.05059154327551117+1.6739618537092438e-16im -1.0196377243957104+4.0211032717849745e-16im -0.24013110723237877-1.1134028611135582e-16im -0.11192459833509985+1.1097897832435794e-14im; 0.21057507448953267-3.1495068151629784e-16im -0.24013110723237877-1.1134028611135582e-16im -0.09413344885723335+3.27941577506768e-17im 0.020092625005744727-2.191989489326832e-15im; -0.2805588371720386+2.852187447161915e-14im -0.11192459833509985+1.1097897832435794e-14im 0.020092625005744727-2.191989489326832e-15im -0.0697017023578413+1.4091189735026237e-14im]
        values, vectors = JosephsonCircuits.autonne_takagi(A)
        @test isapprox(vectors * Diagonal(values) * transpose(vectors), A)
        @test isapprox(vectors * vectors', I(size(A, 1)))

        A = rand(Float64, 4, 4)
        Q, R = qr(A)
        A = Symmetric(Q * Diagonal([1.2, 1.2001, 0, 5e-16]) * Q')
        values, vectors = JosephsonCircuits.autonne_takagi(A)
        @test isapprox(vectors * Diagonal(values) * transpose(vectors), A)
        @test isapprox(vectors * vectors', I(size(A, 1)))

        # https://github.com/JLTastet/TakagiFactorization.jl/issues/4
        A = exp(im * 0.0) * ComplexF64[0.925+0.0im 0.0+0.0im 0.0+0.0im 0.0+0.0im; 0.0+0.0im -0.02399982992272+0.0im -0.00489937871047+0.0im -0.00500042517513+0.0im; 0.0+0.0im -0.00489937871047+0.0im -0.00100017007728+0.0im 0.02449481063548+0.0im; 0.0+0.0im -0.00500042517513+0.0im 0.02449481063548+0.0im 0.0+0.0im]
        values, vectors = JosephsonCircuits.autonne_takagi(A)
        @test isapprox(vectors * Diagonal(values) * transpose(vectors), A)
        @test isapprox(vectors * vectors', I(size(A, 1)))

        # real symmetric matrices stored complex, and the same times a
        # phase, whose unitary Z = U'*conj(V) has its eigenvalues clustered
        # on a short arc: the widest gap is the rest of the circle
        reconstructed = true
        for n in (2, 4, 8), trial in 1:10
            Q = Matrix(qr(randn(n, n)).Q)
            for A in (Complex{Float64}.(Q * Diagonal(1 .+ rand(n)) * Q'),
                    cis(0.3) .* (Q * Diagonal(randn(n)) * Q'))
                A = (A + transpose(A)) / 2
                values, vectors = JosephsonCircuits.autonne_takagi(A)
                reconstructed &= isapprox(vectors * Diagonal(values) * transpose(vectors), A) &&
                    isapprox(vectors * vectors', I(n))
            end
        end
        @test reconstructed

        # a matrix symmetric only to rounding, here by one ulp
        U = JosephsonCircuits.rand_unitary(3)
        A = U * Diagonal([1.0, 2.0, 3.0]) * transpose(U)
        A[1, 2] = complex(nextfloat(real(A[2, 1])), imag(A[2, 1]))
        values, vectors = JosephsonCircuits.autonne_takagi(A)
        @test isapprox(values, [3.0, 2.0, 1.0])
        @test isapprox(vectors * Diagonal(values) * transpose(vectors), A)

        # an empty matrix
        values, vectors = JosephsonCircuits.autonne_takagi(zeros(Complex{Float64}, 0, 0))
        @test isempty(values) && size(vectors) == (0, 0)

        @test_throws(
            ErrorException,
            JosephsonCircuits.autonne_takagi(Complex{Float64}[1 1;-1 1]),
        )
    end

    @testset "autonne_takagi real" begin

        A = Symmetric(rand(Float64, 4, 4))
        values, vectors = JosephsonCircuits.autonne_takagi(A)
        @test isapprox(vectors * Diagonal(values) * transpose(vectors), A)
        @test isapprox(vectors * vectors', I(size(A, 1)))

        A = rand(Float64, 4, 4)
        Q, R = qr(A)
        A = Symmetric(Q * Diagonal([1.2, 1.20000000001, 0.1e-16, -0.5e-16]) * Q')
        values, vectors = JosephsonCircuits.autonne_takagi(A)
        @test isapprox(vectors * Diagonal(values) * transpose(vectors), A)
        @test isapprox(vectors * vectors', I(size(A, 1)))

        # an exactly zero eigenvalue
        A = [1.0 0.0; 0.0 0.0]
        values, vectors = JosephsonCircuits.autonne_takagi(A)
        @test isapprox(vectors * Diagonal(values) * transpose(vectors), A)
        @test isapprox(vectors * vectors', I(size(A, 1)))

        # a matrix symmetric only to rounding, here by one ulp
        Q = Matrix(qr(randn(3, 3)).Q)
        A = Q * Diagonal([1.0, -2.0, 3.0]) * Q'
        A[1, 2] = nextfloat(A[2, 1])
        values, vectors = JosephsonCircuits.autonne_takagi(A)
        @test isapprox(values, [1.0, 2.0, 3.0])
        @test isapprox(vectors * Diagonal(values) * transpose(vectors), A)

        @test_throws(
            ErrorException,
            JosephsonCircuits.autonne_takagi(Float64[1 1;-1 1]),
        )
    end

    @testset "bloch_messiah_block" begin

        # blochmessiah doesn't give correct answers in some cases #26
        # https://github.com/apkille/SymplecticFactorizations.jl/issues/26
        z = 0.1404594873693119
        S = Diagonal([exp(-z), 1, exp(z), 1])
        @test JosephsonCircuits.is_symplectic_block(S)
        O, D, Q = JosephsonCircuits.bloch_messiah_block(S)
        @test isapprox(O * Diagonal(D) * Q, S)
        @test JosephsonCircuits.is_symplectic_block(O)
        @test JosephsonCircuits.is_symplectic_block(Diagonal(D))
        @test JosephsonCircuits.is_symplectic_block(Q)

        S = JosephsonCircuits.rand_symplectic_block(Float64, 4)
        @test JosephsonCircuits.is_symplectic_block(S)
        O, D, Q = JosephsonCircuits.bloch_messiah_block(S)
        @test isapprox(O * Diagonal(D) * Q, S)
        @test JosephsonCircuits.is_symplectic_block(O)
        @test JosephsonCircuits.is_symplectic_block(Diagonal(D))
        @test JosephsonCircuits.is_symplectic_block(Q)
        # the outer factors are orthogonal
        @test JosephsonCircuits.is_orthogonal(O) && JosephsonCircuits.is_orthogonal(Q)

        # single mode squeezers in a real orthogonal basis, which leave x
        # and p uncoupled: the Takagi factorization is then of a real
        # symmetric matrix stored complex
        decomposed = true
        for trial in 1:20
            Or = Matrix(qr(randn(3, 3)).Q)
            r = rand(3)
            S = [Or*Diagonal(exp.(r))*Or' zeros(3, 3); zeros(3, 3) Or*Diagonal(exp.(-r))*Or']
            O, D, Q = JosephsonCircuits.bloch_messiah_block(S)
            decomposed &= isapprox(O * Diagonal(D) * Q, S) &&
                JosephsonCircuits.is_symplectic_block(O) &&
                JosephsonCircuits.is_symplectic_block(Q)
        end
        @test decomposed

        # add a test for this
        # Bloch-Messiah returns incorrect results #728
        # https://github.com/XanaduAI/strawberryfields/issues/728


        # Bloch-messiah decomposition sometimes returns decomposed matrices
        # with permuted rows and columns #14
        # https://github.com/XanaduAI/strawberryfields/issues/14
        S = Float64[1 0 0 0;1 1 0 0;0 0 1 -1;0 0 0 1]
        O, D, Q = JosephsonCircuits.bloch_messiah_block(S)
        @test isapprox(O * Diagonal(D) * Q, S)
        @test JosephsonCircuits.is_symplectic_block(O)
        @test JosephsonCircuits.is_symplectic_block(Diagonal(D))
        @test JosephsonCircuits.is_symplectic_block(Q)
        @test JosephsonCircuits.is_orthogonal(O) && JosephsonCircuits.is_orthogonal(Q)

        @test_throws(
            ErrorException,
            JosephsonCircuits.bloch_messiah_block(Float64[1 1;-1 1]),
        )
    end

    @testset "bloch_messiah_pair" begin

        S = JosephsonCircuits.rand_symplectic_pair(Float64, 4)
        @test JosephsonCircuits.is_symplectic_pair(S)
        O, D, Q = JosephsonCircuits.bloch_messiah_pair(S)
        @test isapprox(O * Diagonal(D) * Q, S)
        @test JosephsonCircuits.is_orthogonal_symplectic_pair(O)
        @test JosephsonCircuits.is_symplectic_pair(Diagonal(D))
        @test JosephsonCircuits.is_orthogonal_symplectic_pair(Q)

    end

    @testset "pre_iwasawa_block" begin
        # real, where F is orthogonal as well as symplectic
        S = JosephsonCircuits.rand_symplectic_block(Float64, 4)
        E, D, F = JosephsonCircuits.pre_iwasawa_block(S)
        @test isapprox(S, E * D * F)
        @test JosephsonCircuits.is_orthogonal_symplectic_block(F)

        # complex
        S = JosephsonCircuits.rand_symplectic_block(Complex{Float64}, 4)
        E, D, F = JosephsonCircuits.pre_iwasawa_block(S)
        @test isapprox(S, E * D * F)
        @test JosephsonCircuits.is_symplectic_block(F)

        @test_throws(
            ErrorException,
            JosephsonCircuits.pre_iwasawa_block(Float64[1 1;-1 1]),
        )
    end

    @testset "pre_iwasawa_pair" begin
        # real, where F is orthogonal as well as symplectic
        S = JosephsonCircuits.rand_symplectic_pair(Float64, 4)
        E, D, F = JosephsonCircuits.pre_iwasawa_pair(S)
        @test isapprox(S, E * D * F)
        @test JosephsonCircuits.is_orthogonal_symplectic_pair(F)

        # complex
        S = JosephsonCircuits.rand_symplectic_pair(Complex{Float64}, 4)
        E, D, F = JosephsonCircuits.pre_iwasawa_pair(S)
        @test isapprox(S, E * D * F)
        @test JosephsonCircuits.is_symplectic_pair(F)

    end

    @testset "iwasawa_block" begin

        # real
        S = JosephsonCircuits.rand_symplectic_block(Float64, 2)
        F = JosephsonCircuits.iwasawa_block(S)
        @test isapprox(F.K * F.K', I(4))
        @test isapprox(F.K * F.A * F.N, S)
        @test JosephsonCircuits.is_symplectic_block(F.K)
        @test JosephsonCircuits.is_symplectic_block(F.A)
        @test JosephsonCircuits.is_symplectic_block(F.N)
        # A is diagonal, and N block upper triangular with a unit upper
        # triangular first block
        @test isdiag(F.A)
        @test iszero(F.N[3:4, 1:2])
        @test istriu(F.N[1:2, 1:2]) && isapprox(diag(F.N[1:2, 1:2]), ones(2))

        # complex
        S = JosephsonCircuits.rand_symplectic_block(Complex{Float64}, 2)
        F = JosephsonCircuits.iwasawa_block(S)
        @test isapprox(F.K * F.K', I(4))
        @test isapprox(F.K * F.A * F.N, S)
        @test JosephsonCircuits.is_symplectic_block(F.K)
        @test JosephsonCircuits.is_symplectic_block(F.A)
        @test JosephsonCircuits.is_symplectic_block(F.N)
        @test isdiag(F.A)
        @test iszero(F.N[3:4, 1:2])
        @test istriu(F.N[1:2, 1:2]) && isapprox(diag(F.N[1:2, 1:2]), ones(2))

        @test_throws(
            ErrorException,
            JosephsonCircuits.iwasawa_block(Float64[1 1;-1 1]),
        )
    end

    @testset "iwasawa_pair" begin

        # real
        S = JosephsonCircuits.rand_symplectic_pair(Float64, 2)
        F = JosephsonCircuits.iwasawa_pair(S)
        @test isapprox(F.K * F.K', I(4))
        @test isapprox(F.K * F.A * F.N, S)
        @test JosephsonCircuits.is_symplectic_pair(F.K)
        @test JosephsonCircuits.is_symplectic_pair(F.A)
        @test JosephsonCircuits.is_symplectic_pair(F.N)

        # complex
        S = JosephsonCircuits.rand_symplectic_pair(Complex{Float64}, 2)
        F = JosephsonCircuits.iwasawa_pair(S)
        @test isapprox(F.K * F.K', I(4))
        @test isapprox(F.K * F.A * F.N, S)
        @test JosephsonCircuits.is_symplectic_pair(F.K)
        @test JosephsonCircuits.is_symplectic_pair(F.A)
        @test JosephsonCircuits.is_symplectic_pair(F.N)

    end


    @testset "iwasawa_bogoliubov_block" begin

        # complex
        S = JosephsonCircuits.rand_bogoliubov_block(2)
        F = JosephsonCircuits.iwasawa_bogoliubov_block(S)
        @test isapprox(F.K * F.K', I(4))
        @test isapprox(F.K * F.A * F.N, S)
        @test JosephsonCircuits.is_bogoliubov_block(F.K)
        @test JosephsonCircuits.is_bogoliubov_block(F.A)
        @test JosephsonCircuits.is_bogoliubov_block(F.N)

    end

    @testset "iwasawa_bogoliubov_pair" begin

        # complex
        S = JosephsonCircuits.rand_bogoliubov_pair(2)
        F = JosephsonCircuits.iwasawa_bogoliubov_pair(S)
        @test isapprox(F.K * F.K', I(4))
        @test isapprox(F.K * F.A * F.N, S)
        @test JosephsonCircuits.is_bogoliubov_pair(F.K)
        @test JosephsonCircuits.is_bogoliubov_pair(F.A)
        @test JosephsonCircuits.is_bogoliubov_pair(F.N)

    end

    @testset "symplectic_normal_form_pair" begin

        A = randn(Float64, 4, 4)
        Aa = (A - A') / 2
        Q = JosephsonCircuits.symplectic_normal_form_pair(Aa)
        Omega = JosephsonCircuits.symplectic_form_pair(2)
        @test isapprox(Aa, Q * Omega * Q')

        # a singular matrix, of rank two by construction
        O = Matrix(qr(randn(4, 4)).Q)
        Aa1 = O * [0 1.3 0 0; -1.3 0 0 0; 0 0 0 0; 0 0 0 0] * transpose(O)
        Q1 = JosephsonCircuits.symplectic_normal_form_pair(Aa1)
        @test isapprox(Aa1, Q1 * Omega * Q1')

        # errors
        @test_throws(
            ErrorException,
            JosephsonCircuits.symplectic_normal_form_pair([1 1;1 1]),
        )
        @test_throws(
            ErrorException,
            JosephsonCircuits.symplectic_normal_form_pair([0 0 1;0 0 0;-1 0 0]),
        )

    end

    @testset "symplectic_normal_form_block" begin

        A = randn(Float64, 4, 4)
        Aa = (A - A') / 2
        Q = JosephsonCircuits.symplectic_normal_form_block(Aa)
        Omega = JosephsonCircuits.symplectic_form_block(2)
        @test isapprox(Aa, Q * Omega * Q')

    end

    @testset "halmos dilation" begin

        S = [0.1 0;0 0.1]
        # test that it gives the same result as the simpler formula
        @test isapprox(JosephsonCircuits.halmos_dilation(S),[S sqrt(I(size(S,1)) - S*S');sqrt(I(size(S,1)) - S*S') -S'])

        # test that the Halmos dilation produces a symplectic matrix after
        # conversion of the scattering parameter matrix to a symplectic matrix
        @test JosephsonCircuits.is_symplectic_block(JosephsonCircuits.scattering_to_ladder_block(JosephsonCircuits.halmos_dilation(S),[1,1,1,1]))
    
        # test that the Halmos dilation produces a unitary scattering
        # parameter matrix
        @test JosephsonCircuits.is_unitary(JosephsonCircuits.halmos_dilation([0.1 0;0 0.1]))

        # the closed form of the dilation of a passive matrix which is not
        # normal
        S = 0.5*JosephsonCircuits.rand_unitary(3)*Diagonal([0.2, 0.9, 1.0])*JosephsonCircuits.rand_unitary(3)
        @test isapprox(JosephsonCircuits.halmos_dilation(S),
            [S sqrt(Hermitian(I - S*S')); sqrt(Hermitian(I - S'*S)) -S'])

        # a lossless matrix, whose singular values are one up to rounding,
        # dilates to a unitary one; a matrix with gain is refused
        dilated = true
        for trial in 1:20
            U = JosephsonCircuits.rand_unitary(4)
            dilated &= JosephsonCircuits.is_unitary(JosephsonCircuits.halmos_dilation(U))
        end
        @test dilated
        @test_throws(
            ArgumentError,
            JosephsonCircuits.halmos_dilation([2.0 0;0 0.5]),
        )

        # a rectangular contraction dilates to a unitary matrix of the sum
        # of its dimensions, of the same closed form
        for (n, m) in ((2, 3), (3, 2))
            S = randn(Complex{Float64}, n, m)
            S = 0.9 * S / opnorm(S)
            U = JosephsonCircuits.halmos_dilation(S)
            @test JosephsonCircuits.is_unitary(U)
            @test isapprox(U, [S sqrt(Hermitian(I - S*S')); sqrt(Hermitian(I - S'*S)) -S'])
        end

    end

    @testset "Ymin_from_X_quadrature_pair and Ymin_from_X_quadrature_block" begin

        X = rand(Float64,4,4)
        for method in 1:3
            @test JosephsonCircuits.is_cptp_quadrature_pair(X,JosephsonCircuits.Ymin_from_X_quadrature_pair(X;method=method))
            @test JosephsonCircuits.is_cptp_quadrature_block(X,JosephsonCircuits.Ymin_from_X_quadrature_block(X;method=method))
        end
        # the three methods are three routes to |im*(Ω - X*Ω*X')|
        for method in 2:3
            @test isapprox(JosephsonCircuits.Ymin_from_X_quadrature_pair(X; method = method),
                JosephsonCircuits.Ymin_from_X_quadrature_pair(X; method = 1))
            @test isapprox(JosephsonCircuits.Ymin_from_X_quadrature_block(X; method = method),
                JosephsonCircuits.Ymin_from_X_quadrature_block(X; method = 1))
        end

        @test_throws(
            ErrorException,
            JosephsonCircuits.is_cptp_quadrature_pair(X,JosephsonCircuits.Ymin_from_X_quadrature_pair(X;method=4)),
            )
    end

    @testset "A_B_to_symplectic" begin

        S0 = JosephsonCircuits.rand_symplectic_pair(Float64, 3)
        A = S0[1:2, 1:2]
        B = S0[1:2, 3:end]
        S1 = JosephsonCircuits.A_B_to_symplectic_pair(A, B)
        @test JosephsonCircuits.is_symplectic_pair(S1)

        nsys = 2
        nenv = 2 * nsys
        S0 = JosephsonCircuits.rand_symplectic_pair(Float64, nsys + nenv)
        A = S0[1:2*nsys, 1:2*nsys]
        B = S0[1:2*nsys, 2*nsys+1:end]

        X = A
        Y = B * B'
        B1 = JosephsonCircuits.B_from_X_Y_quadrature_pair(X, Y)
        S1 = JosephsonCircuits.A_B_to_symplectic_pair(A, B1)
        @test JosephsonCircuits.is_symplectic_pair(S1)

        # rows of rank below 2n cannot be completed
        @test_throws(ArgumentError,
            JosephsonCircuits.A_B_to_symplectic_pair(zeros(2, 2), zeros(2, 4)))
    end

    @testset "B_from_X_Y_quadrature_block" begin

        @test_throws(
            ErrorException,
            JosephsonCircuits.B_from_X_Y_quadrature_block([1 0;0 1],[1 0;0 -1]),
        )
    end

    # the matrix S of the system and its environment realizes the map: its
    # block on the system is X, and the environment in the vacuum, whose
    # covariance is the identity, adds the noise B*B' = Y through its block
    # B from the environment to the system

    @testset "X_Y_to_symplectic_pair" begin

        X = rand(Float64,4,4)
        Y = JosephsonCircuits.Ymin_from_X_quadrature_pair(X)
        S = JosephsonCircuits.X_Y_to_symplectic_pair(X,Y)
        @test JosephsonCircuits.is_symplectic_pair(S)
        @test S[1:4, 1:4] == X
        @test isapprox(S[1:4, 5:end] * S[1:4, 5:end]', Y)

        X, Y = JosephsonCircuits.rand_cptp_quadrature_pair(2)
        S = JosephsonCircuits.X_Y_to_symplectic_pair(X,Y)
        @test JosephsonCircuits.is_symplectic_pair(S)
        @test S[1:4, 1:4] == X
        @test isapprox(S[1:4, 5:end] * S[1:4, 5:end]', Y)

    end

    @testset "X_Y_to_symplectic_block" begin

        # the system, two of six modes, in block order
        sys = [1:2; 7:8]
        env = setdiff(1:12, sys)

        X = rand(Float64,4,4)
        Y = JosephsonCircuits.Ymin_from_X_quadrature_block(X)
        S = JosephsonCircuits.X_Y_to_symplectic_block(X,Y)
        @test JosephsonCircuits.is_symplectic_block(S)
        @test S[sys, sys] == X
        @test isapprox(S[sys, env] * S[sys, env]', Y)

        X, Y = JosephsonCircuits.rand_cptp_quadrature_block(2)
        S = JosephsonCircuits.X_Y_to_symplectic_block(X,Y)
        @test JosephsonCircuits.is_symplectic_block(S)
        @test S[sys, sys] == X
        @test isapprox(S[sys, env] * S[sys, env]', Y)

    end

    @testset "X_Y_to_bogoliubov_pair" begin

        X, Y = JosephsonCircuits.rand_cptp_ladder_pair(2)
        S = JosephsonCircuits.X_Y_to_bogoliubov_pair(X,Y)
        @test JosephsonCircuits.is_bogoliubov_pair(S)
        @test isapprox(S[1:4, 1:4], X)
        @test isapprox(S[1:4, 5:end] * S[1:4, 5:end]', Y)

    end

    @testset "X_Y_to_bogoliubov_block" begin

        sys = [1:2; 7:8]
        env = setdiff(1:12, sys)
        X, Y = JosephsonCircuits.rand_cptp_ladder_block(2)
        S = JosephsonCircuits.X_Y_to_bogoliubov_block(X,Y)
        @test JosephsonCircuits.is_bogoliubov_block(S)
        @test isapprox(S[sys, sys], X)
        @test isapprox(S[sys, env] * S[sys, env]', Y)

    end

    @testset "wmatrix" begin
        # the in-place form fills and returns a matrix of the size it checks
        w = zeros(2, 3)
        @test JosephsonCircuits.wmatrix!(w, 0.1:0.1:0.3, (1.0,), [(1,), (-1,)]) === w
        @test w == JosephsonCircuits.wmatrix(0.1:0.1:0.3, (1.0,), [(1,), (-1,)])
        @test_throws(DimensionMismatch,
            JosephsonCircuits.wmatrix!(zeros(3, 3), 0.1:0.1:0.3, (1.0,), [(1,), (-1,)]))
        # signal frequencies in a vector, and an integer pump
        @test isapprox(JosephsonCircuits.wmatrix([0.1, 0.2], (1,), [(1,), (-1,)]),
            [1.1 1.2; -0.9 -0.8])
    end

    @testset "interpolate_scattering" begin

        # a notch, whose transmission zero at w = 5 lies between samples,
        # alone and behind a delay whose phase winds three times across the
        # band, against the exact response at the midpoints
        w0 = collect(4.0:0.01:6.0)
        wm = (w0[1:end-1] .+ w0[2:end]) ./ 2
        for tau in (0.0, 20.0)
            notch(w) = cis(-w * tau) * (w - 5) / (w - 5 + 0.05im)
            Sm = JosephsonCircuits.interpolate_scattering(w0,
                reshape(notch.(w0), 1, 1, :), wm)
            @test maximum(abs, Sm[1, 1, :] .- notch.(wm)) < 1e-3
        end

        # the idlers of a 4 to 8 GHz table pumped at 12 GHz, which rounding
        # puts up to an ulp outside the table, take the values at its edges
        w0 = collect(2pi .* (4e9:0.01e9:8e9))
        w = JosephsonCircuits.wmatrix(w0, (2pi * 12e9,), [(0,), (-1,)])
        Sd = reshape(cis.(-w0 ./ 1e10), 1, 1, :)
        @test isapprox(JosephsonCircuits.interpolate_scattering(w0, Sd, w)[1, 1, :, :],
            cis.(-w ./ 1e10))
        @test isapprox(JosephsonCircuits.interpolate_scattering(w0, Sd, w;
            extrap = true)[1, 1, :, :], cis.(-w ./ 1e10))

        # real samples give complex values
        Sr = JosephsonCircuits.interpolate_scattering([1.0, 2.0, 3.0],
            reshape([1.0, -1.0, 1.0], 1, 1, :), [1.0, 2.0])
        @test Sr[1, 1, :] ≈ [1.0, -1.0] && eltype(Sr) == Complex{Float64}

        # fewer than three samples, and samples out of order
        @test_throws ArgumentError JosephsonCircuits.interpolate_scattering(
            [1.0, 2.0], ones(1, 1, 2), [1.5])
        @test_throws ArgumentError JosephsonCircuits.interpolate_scattering(
            [3.0, 2.0, 1.0], ones(1, 1, 3), [1.5])

        # a matched delay line, whose phase winds four times across the
        # band, between its samples: the exact response, conjugated at
        # negative frequencies, and the value given outside the band
        w0 = collect(range(1.0, 10.0, length = 101))
        tau = 3.0
        Sd = zeros(Complex{Float64}, 2, 2, length(w0))
        Sd[2, 1, :] .= cis.(-w0 .* tau)
        Sd[1, 2, :] .= Sd[2, 1, :]
        wm = (w0[1:end-1] .+ w0[2:end]) ./ 2
        Sm = JosephsonCircuits.interpolate_scattering(w0, Sd, wm)
        @test isapprox(Sm[2, 1, :], cis.(-wm .* tau))
        @test isapprox(Sm[1, 2, :], cis.(-wm .* tau))
        Sm = JosephsonCircuits.interpolate_scattering(w0, Sd, -wm)
        @test isapprox(Sm[2, 1, :], cis.(wm .* tau))
        Sm = JosephsonCircuits.interpolate_scattering(w0, Sd, [20.0, -20.0];
            extrap = true, extrap_value = 0.5im)
        @test Sm[2, 1, :] == [0.5im, -0.5im]

        # test with incorrect dimensions
        w = 0.01:0.01:1.0
        @test_throws(
            ErrorException,
            JosephsonCircuits.interpolate_scattering(w,randn(Complex{Float64},2),w),
        )

        @test_throws(
            ErrorException,
            JosephsonCircuits.interpolate_scattering(w,randn(Complex{Float64},2,2,2*length(w)),w),
        )
    end

end
