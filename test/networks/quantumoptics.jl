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
            ErrorException(lazy"The dimensions of the input matrix must be even."),
            JosephsonCircuits.is_symplectic_pair(rand(9,9)),
        )
        @test_throws(
            ErrorException(lazy"The dimensions of the input matrix must be even."),
            JosephsonCircuits.is_symplectic_block(rand(9,9)),
        )

    end

    @testset "is_cptp" begin

        @test !JosephsonCircuits.is_cptp_quadrature_pair([1 0;0 1],[1 0;1 1])
        @test JosephsonCircuits.is_cptp_quadrature_pair([1 0;0 1],[0 0;0 0])

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

    @testset "random symplectic" begin

        @test JosephsonCircuits.is_symplectic_pair(JosephsonCircuits.rand_symplectic_pair(4))
        @test JosephsonCircuits.is_symplectic_pair(JosephsonCircuits.rand_symplectic_pair(Float64,4))
        @test JosephsonCircuits.is_symplectic_pair(JosephsonCircuits.rand_symplectic_pair(Complex{Float64},4))

        @test JosephsonCircuits.is_symplectic_block(JosephsonCircuits.rand_symplectic_block(4))
        @test JosephsonCircuits.is_symplectic_block(JosephsonCircuits.rand_symplectic_block(Float64,4))
        @test JosephsonCircuits.is_symplectic_block(JosephsonCircuits.rand_symplectic_block(Complex{Float64},4))

    end

    @testset "random orthogonal symplectic" begin

        @test JosephsonCircuits.is_orthogonal_symplectic_pair(JosephsonCircuits.rand_orthogonal_symplectic_pair(4))
        @test JosephsonCircuits.is_orthogonal_symplectic_pair(JosephsonCircuits.rand_orthogonal_symplectic_pair(Float64,4))
        @test JosephsonCircuits.is_orthogonal_symplectic_pair(JosephsonCircuits.rand_orthogonal_symplectic_pair(Complex{Float64},4))

        @test JosephsonCircuits.is_orthogonal_symplectic_block(JosephsonCircuits.rand_orthogonal_symplectic_block(4))
        @test JosephsonCircuits.is_orthogonal_symplectic_block(JosephsonCircuits.rand_orthogonal_symplectic_block(Float64,4))
        @test JosephsonCircuits.is_orthogonal_symplectic_block(JosephsonCircuits.rand_orthogonal_symplectic_block(Complex{Float64},4))

    end

    @testset "random positive definite symplectic" begin

        @test JosephsonCircuits.is_positive_definite_symplectic_pair(JosephsonCircuits.rand_positive_definite_symplectic_pair(4))
        @test JosephsonCircuits.is_positive_definite_symplectic_pair(JosephsonCircuits.rand_positive_definite_symplectic_pair(Float64,4))
        @test JosephsonCircuits.is_positive_definite_symplectic_pair(JosephsonCircuits.rand_positive_definite_symplectic_pair(Complex{Float64},4))

        @test JosephsonCircuits.is_positive_definite_symplectic_block(JosephsonCircuits.rand_positive_definite_symplectic_block(4))
        @test JosephsonCircuits.is_positive_definite_symplectic_block(JosephsonCircuits.rand_positive_definite_symplectic_block(Float64,4))
        @test JosephsonCircuits.is_positive_definite_symplectic_block(JosephsonCircuits.rand_positive_definite_symplectic_block(Complex{Float64},4))

    end

    @testset "random conjugate symplectic" begin

        @test JosephsonCircuits.is_conjugate_symplectic_pair(JosephsonCircuits.rand_conjugate_symplectic_pair(4))
        @test JosephsonCircuits.is_conjugate_symplectic_pair(JosephsonCircuits.rand_conjugate_symplectic_pair(Float64,4))
        @test JosephsonCircuits.is_conjugate_symplectic_pair(JosephsonCircuits.rand_conjugate_symplectic_pair(Complex{Float64},4))

        @test JosephsonCircuits.is_conjugate_symplectic_block(JosephsonCircuits.rand_conjugate_symplectic_block(4))
        @test JosephsonCircuits.is_conjugate_symplectic_block(JosephsonCircuits.rand_conjugate_symplectic_block(Float64,4))
        @test JosephsonCircuits.is_conjugate_symplectic_block(JosephsonCircuits.rand_conjugate_symplectic_block(Complex{Float64},4))

    end

    @testset "random bogoliubov" begin

        @test JosephsonCircuits.is_bogoliubov_pair(JosephsonCircuits.rand_bogoliubov_pair(4))
        @test JosephsonCircuits.is_bogoliubov_pair(JosephsonCircuits.rand_bogoliubov_pair(Float64,4))
        @test JosephsonCircuits.is_bogoliubov_pair(JosephsonCircuits.rand_bogoliubov_pair(Complex{Float64},4))

        @test JosephsonCircuits.is_bogoliubov_block(JosephsonCircuits.rand_bogoliubov_block(4))
        @test JosephsonCircuits.is_bogoliubov_block(JosephsonCircuits.rand_bogoliubov_block(Float64,4))
        @test JosephsonCircuits.is_bogoliubov_block(JosephsonCircuits.rand_bogoliubov_block(Complex{Float64},4))

    end

    @testset "random orthogonal bogoliubov" begin

        @test JosephsonCircuits.is_orthogonal_bogoliubov_pair(JosephsonCircuits.rand_orthogonal_bogoliubov_pair(4))
        @test JosephsonCircuits.is_orthogonal_bogoliubov_pair(JosephsonCircuits.rand_orthogonal_bogoliubov_pair(Float64,4))
        @test JosephsonCircuits.is_orthogonal_bogoliubov_pair(JosephsonCircuits.rand_orthogonal_bogoliubov_pair(Complex{Float64},4))

        @test JosephsonCircuits.is_orthogonal_bogoliubov_block(JosephsonCircuits.rand_orthogonal_bogoliubov_block(4))
        @test JosephsonCircuits.is_orthogonal_bogoliubov_block(JosephsonCircuits.rand_orthogonal_bogoliubov_block(Float64,4))
        @test JosephsonCircuits.is_orthogonal_bogoliubov_block(JosephsonCircuits.rand_orthogonal_bogoliubov_block(Complex{Float64},4))

    end

    @testset "random pseudo-unitary" begin

        @test JosephsonCircuits.is_pseudo_unitary_pair(JosephsonCircuits.rand_pseudo_unitary_pair(4))
        @test JosephsonCircuits.is_pseudo_unitary_pair(JosephsonCircuits.rand_pseudo_unitary_pair(Float64,4))
        @test JosephsonCircuits.is_pseudo_unitary_pair(JosephsonCircuits.rand_pseudo_unitary_pair(Complex{Float64},4))
        
        @test JosephsonCircuits.is_pseudo_unitary_block(JosephsonCircuits.rand_pseudo_unitary_block(4))
        @test JosephsonCircuits.is_pseudo_unitary_block(JosephsonCircuits.rand_pseudo_unitary_block(Float64,4))
        @test JosephsonCircuits.is_pseudo_unitary_block(JosephsonCircuits.rand_pseudo_unitary_block(Complex{Float64},4))

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
            DimensionMismatch("The length of the bogoliubov vector must be double that of the scattering parameter vector."),
            JosephsonCircuits.scattering_to_ladder_pair!(zeros(10),zeros(10),ones(5)),
        )

        @test_throws(
            DimensionMismatch("The length of the bogoliubov vector must be double that of the scattering parameter vector."),
            JosephsonCircuits.scattering_to_ladder_block!(zeros(10),zeros(10),ones(5)),
        )

        @test_throws(
            DimensionMismatch("The length of the symplectic vector must be double that of the scattering parameter vector."),
            JosephsonCircuits.scattering_to_quadrature_pair!(zeros(10),zeros(10),ones(5)),
        )

        @test_throws(
            DimensionMismatch("The length of the symplectic vector must be double that of the scattering parameter vector."),
            JosephsonCircuits.scattering_to_quadrature_block!(zeros(10),zeros(10),ones(5)),
        )

        @test_throws(
            DimensionMismatch("Length of scattering vector must be integer multiples of the number of modes."),
            JosephsonCircuits.scattering_to_ladder_pair!(zeros(10),zeros(5),ones(6)),
        )

        @test_throws(
            DimensionMismatch("Length of scattering vector must be integer multiples of the number of modes."),
            JosephsonCircuits.scattering_to_ladder_block!(zeros(10),zeros(5),ones(6)),
        )

        @test_throws(
            DimensionMismatch("Length of scattering vector must be integer multiples of the number of modes."),
            JosephsonCircuits.scattering_to_quadrature_pair!(zeros(10),zeros(5),ones(6)),
        )

        @test_throws(
            DimensionMismatch("Length of scattering vector must be integer multiples of the number of modes."),
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
                # basis of its rows and of its columns
                @test isapprox(JosephsonCircuits.ladder_to_quadrature_pair(
                    JosephsonCircuits.scattering_to_ladder_pair(S, w)),
                    JosephsonCircuits.scattering_to_quadrature_pair(S, w))
                @test isapprox(JosephsonCircuits.quadrature_to_ladder_block(
                    JosephsonCircuits.scattering_to_quadrature_block(S, w)),
                    JosephsonCircuits.scattering_to_ladder_block(S, w))
            end
        end

        X = randn(10)
        w = sign.(randn(5))
        @test isequal(
            JosephsonCircuits.scattering_to_quadrature_block(X,w),
            JosephsonCircuits.scattering_to_quadrature_block(Complex.(X),w)
        )

        X = randn(10)
        w = sign.(randn(5))
        @test isequal(
            JosephsonCircuits.scattering_to_quadrature_pair(X,w),
            JosephsonCircuits.scattering_to_quadrature_pair(Complex.(X),w)
        )
    end

    @testset "ladder_to_scattering_pair" begin

        # complex floating point input
        X = JosephsonCircuits.rand_unitary(10)
        w=sign.(randn(size(X,1)))
        S = JosephsonCircuits.scattering_to_ladder_pair(X,w)
        @test isapprox(
            JosephsonCircuits.ladder_to_scattering_pair(S,w),
            X,
        )

        S = JosephsonCircuits.rand_symplectic_pair(5)
        @test_throws(
            ErrorException(lazy"Error in Bogoliubov to scattering parameter conversion larger than `atol` and `rtol`."),
            JosephsonCircuits.ladder_to_scattering_pair(S,w),
        )

        # complex floating point input carrying rounding is converted within
        # the default tolerance, and the tolerances can be given
        X = JosephsonCircuits.rand_unitary(6)
        w = [1.0, -1.0, 1.0]
        U = JosephsonCircuits.rand_unitary(12)
        S = JosephsonCircuits.scattering_to_ladder_pair(X,w)*(U*U')
        @test isapprox(JosephsonCircuits.ladder_to_scattering_pair(S,w), X)
        @test isapprox(JosephsonCircuits.ladder_to_scattering_pair(S,w;rtol=1e-10), X)

        # real floating point input
        X = JosephsonCircuits.rand_unitary(Float64,10)
        w=sign.(randn(size(X,1)))
        S = JosephsonCircuits.scattering_to_ladder_pair(X,w)
        @test isapprox(
            JosephsonCircuits.ladder_to_scattering_pair(S,w),
            X,
        )


        # error 1
        @test_throws(
            ErrorException(lazy"Error in Bogoliubov to scattering parameter conversion larger than `atol` and `rtol`."),
            JosephsonCircuits.ladder_to_scattering_pair([1 0 0 0;0 -1 0 0;0 0 1 0;0 0 0 1],[1,1,1,1]),
        )
        # error 2
        @test_throws(
            ErrorException(lazy"Error in Bogoliubov to scattering parameter conversion larger than `atol` and `rtol`."),
            JosephsonCircuits.ladder_to_scattering_pair([1 1 0 0;0 1 0 0;0 0 1 0;0 0 0 1],[1,1,1,1]),
        )
        # error 3
        @test_throws(
            ErrorException(lazy"Error in Bogoliubov to scattering parameter conversion larger than `atol` and `rtol`."),
            JosephsonCircuits.ladder_to_scattering_pair([1 0 0 0;1 1 0 0;0 0 1 0;0 0 0 1],[1,1,1,1]),
        )

    end

    @testset "ladder_to_scattering_block" begin

        # complex floating point input
        X = JosephsonCircuits.rand_unitary(10)
        w=sign.(randn(size(X,1)))
        S = JosephsonCircuits.scattering_to_ladder_block(X,w)
        @test isapprox(
            JosephsonCircuits.ladder_to_scattering_block(S,w),
            X,
        )

        S = JosephsonCircuits.rand_symplectic_block(5)
        @test_throws(
            ErrorException(lazy"Error in Bogoliubov to scattering parameter conversion larger than `atol` and `rtol`."),
            JosephsonCircuits.ladder_to_scattering_block(S,w),
        )

        # complex floating point input carrying rounding is converted within
        # the default tolerance, and the tolerances can be given
        X = JosephsonCircuits.rand_unitary(6)
        w = [1.0, -1.0, 1.0]
        U = JosephsonCircuits.rand_unitary(12)
        S = JosephsonCircuits.scattering_to_ladder_block(X,w)*(U*U')
        @test isapprox(JosephsonCircuits.ladder_to_scattering_block(S,w), X)
        @test isapprox(JosephsonCircuits.ladder_to_scattering_block(S,w;rtol=1e-10), X)

        # real floating point input
        X = JosephsonCircuits.rand_unitary(Float64,10)
        w=sign.(randn(size(X,1)))
        S = JosephsonCircuits.scattering_to_ladder_block(X,w)
        @test isapprox(
            JosephsonCircuits.ladder_to_scattering_block(S,w),
            X,
        )


        # error 1
        @test_throws(
            ErrorException(lazy"Error in Bogoliubov to scattering parameter conversion larger than `atol` and `rtol`."),
            JosephsonCircuits.ladder_to_scattering_block([1 0 0 0;0 -1 0 0;0 0 1 0;0 0 0 1],[1,1,1,1]),
        )
        # error 2
        @test_throws(
            ErrorException(lazy"Error in Bogoliubov to scattering parameter conversion larger than `atol` and `rtol`."),
            JosephsonCircuits.ladder_to_scattering_block([1 0 0 0;0 1 0 1;0 0 1 0;0 0 0 1],[1,1,1,1]),
        )
        # error 3
        @test_throws(
            ErrorException(lazy"Error in Bogoliubov to scattering parameter conversion larger than `atol` and `rtol`."),
            JosephsonCircuits.ladder_to_scattering_block([1 0 0 0;0 1 0 0;0 0 1 0;1 0 0 1],[1,1,1,1]),
        )

    end

    @testset "quadrature_to_scattering_pair" begin

        # complex floating point input
        X = JosephsonCircuits.rand_unitary(10)
        w=sign.(randn(size(X,1)))
        S = JosephsonCircuits.scattering_to_quadrature_pair(X,w)
        @test isapprox(
            JosephsonCircuits.quadrature_to_scattering_pair(S,w),
            X,
        )
        @test isapprox(
            JosephsonCircuits.quadrature_to_scattering_pair(Complex.(S),w),
            X,
        )

        S = JosephsonCircuits.rand_symplectic_pair(5)
        @test_throws(
            ErrorException(lazy"Error in symplectic to scattering parameter conversion larger than `atol` and `rtol`."),
            JosephsonCircuits.quadrature_to_scattering_pair(S,w),
        )

        # complex floating point input carrying rounding is converted within
        # the default tolerance, and the tolerances can be given
        X = JosephsonCircuits.rand_unitary(6)
        w = [1.0, -1.0, 1.0]
        U = JosephsonCircuits.rand_unitary(12)
        S = JosephsonCircuits.scattering_to_quadrature_pair(X,w)*(U*U')
        @test isapprox(JosephsonCircuits.quadrature_to_scattering_pair(S,w), X)
        @test isapprox(JosephsonCircuits.quadrature_to_scattering_pair(S,w;rtol=1e-10), X)

        # real floating point input
        X = JosephsonCircuits.rand_unitary(Float64,10)
        w=sign.(randn(size(X,1)))
        S = JosephsonCircuits.scattering_to_quadrature_pair(X,w)
        @test isapprox(
            JosephsonCircuits.quadrature_to_scattering_pair(S,w),
            X,
        )


    end

    @testset "quadrature_to_scattering_block" begin

        # complex floating point input
        X = JosephsonCircuits.rand_unitary(10)
        w=sign.(randn(size(X,1)))
        S = JosephsonCircuits.scattering_to_quadrature_block(X,w)
        @test isapprox(
            JosephsonCircuits.quadrature_to_scattering_block(S,w),
            X,
        )
        @test isapprox(
            JosephsonCircuits.quadrature_to_scattering_block(Complex.(S),w),
            X,
        )

        S = JosephsonCircuits.rand_symplectic_block(5)
        @test_throws(
            ErrorException(lazy"Error in symplectic to scattering parameter conversion larger than `atol` and `rtol`."),
            JosephsonCircuits.quadrature_to_scattering_block(S,w),
        )

        # complex floating point input carrying rounding is converted within
        # the default tolerance, and the tolerances can be given
        X = JosephsonCircuits.rand_unitary(6)
        w = [1.0, -1.0, 1.0]
        U = JosephsonCircuits.rand_unitary(12)
        S = JosephsonCircuits.scattering_to_quadrature_block(X,w)*(U*U')
        @test isapprox(JosephsonCircuits.quadrature_to_scattering_block(S,w), X)
        @test isapprox(JosephsonCircuits.quadrature_to_scattering_block(S,w;rtol=1e-10), X)

        # real floating point input
        X = JosephsonCircuits.rand_unitary(Float64,10)
        w=sign.(randn(size(X,1)))
        S = JosephsonCircuits.scattering_to_quadrature_block(X,w)
        @test isapprox(
            JosephsonCircuits.quadrature_to_scattering_block(S,w),
            X,
        )


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
            DimensionMismatch("The number of scattering indices 6 of an axis must be a multiple of the number of modes 4."),
            JosephsonCircuits.ports_modes_to_modes_ports_scattering(rand(6, 6), 4),
        )
        @test_throws(
            DimensionMismatch("The number of scattering indices 3 of an axis must be a multiple of the number of modes 2."),
            JosephsonCircuits.modes_ports_to_ports_modes_block(rand(6, 6), 2),
        )
        @test_throws(
            DimensionMismatch("The length 5 of an axis of a pair or block form must be even."),
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

        for M in [JosephsonCircuits.rand_positive_definite(4), JosephsonCircuits.rand_positive_semi_definite(2,2)]
            d, S = JosephsonCircuits.williamson_pair(M)
            @test JosephsonCircuits.is_symplectic_pair(S)
            # serafini convention
            # @test isapprox(transpose(S)*Diagonal(d)*S,M)
            # literature convention
            @test isapprox(S*Diagonal(d)*transpose(S),M)
        end

       # #  # add a test for this matrix for both block and pair
       #  M = Symmetric([
       #     1.166893242623039673e-01  -1.048482006204835837e-02  -3.493446128036903353e-02   1.565614325188508238e-01;
       #    -1.048482006204835837e-02   4.805855935380488053e-01  -1.069468907857842987e+00  -3.770319099460901491e-01;
       #    -3.493446128036903353e-02  -1.069468907857842987e+00   2.409089316529508640e+00   7.648117805967328264e-01;
       #     1.565614325188508238e-01  -3.770319099460901491e-01   7.648117805967328264e-01   4.847266530147685271e-01;
       # ])


        # julia> M = Float64[1 0 0 0;0 1 0 0;0 0 0 0;0 0 0 0]
        # 4×4 Matrix{Float64}:
        #  1.0  0.0  0.0  0.0
        #  0.0  1.0  0.0  0.0
        #  0.0  0.0  0.0  0.0
        #  0.0  0.0  0.0  0.0

        # @test_throws(
        #     ErrorException(lazy"The rank must be even."),
        #     JosephsonCircuits.williamson_pair([1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 0]),
        # )
    end


    @testset "williamson block " begin

        for M in [JosephsonCircuits.rand_positive_definite(4), JosephsonCircuits.rand_positive_semi_definite(2,2)]
            d, S = JosephsonCircuits.williamson_block(M)
            @test JosephsonCircuits.is_symplectic_block(S)
            # serafini convention
            # @test isapprox(transpose(S)*Diagonal(d)*S,M)
            # literature convention
            @test isapprox(S*Diagonal(d)*transpose(S),M)
        end
    end

    @testset "cholesky_williamson" begin

        L1, rankL1 = JosephsonCircuits.cholesky_williamson([1 0;0 1])
        L2, rankL2 = JosephsonCircuits.cholesky_williamson(JosephsonCircuits.SparseArrays.sparse([1 0;0 1]))

        @test isapprox(L1,L2)
        @test isapprox(rankL1,rankL2)

        @test_throws(
            ErrorException(lazy"Cholesky factorization has failed. Input matrix is not positive semi-definite."),
            JosephsonCircuits.cholesky_williamson([1 0 0 0;0 1 0 0;0 0 0 0;0 0 0 -1]),
        )

        @test_throws(
            ErrorException(lazy"Cholesky factorization has failed. Input matrix is not positive semi-definite."),
            JosephsonCircuits.cholesky_williamson(JosephsonCircuits.SparseArrays.sparse([1 0 0 0;0 1 0 0;0 0 0 0;0 0 0 -1])),
        )

        # this fails but with the wrong error message because of the rank
        # reduction kludge in cholesky_williamson.
        @test_throws(
            ErrorException(lazy"Cholesky factorization has failed. Input matrix is not positive semi-definite."),
            JosephsonCircuits.cholesky_williamson([1 0;0 0]),
        )
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

        @test_throws(
            ErrorException(lazy"M must be symmetric."),
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

        @test_throws(
            ErrorException(lazy"M must be symmetric."),
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

        @test_throws(
            ErrorException(lazy"A must be symplectic."),
            JosephsonCircuits.bloch_messiah_block(Float64[1 1;-1 1]),
        )
    end

    @testset "bloch_messiah_pair" begin

        S = JosephsonCircuits.rand_symplectic_pair(Float64, 4)
        @test JosephsonCircuits.is_symplectic_pair(S)
        O, D, Q = JosephsonCircuits.bloch_messiah_pair(S)
        @test isapprox(O * Diagonal(D) * Q, S)
        @test JosephsonCircuits.is_symplectic_pair(O)
        @test JosephsonCircuits.is_symplectic_pair(Diagonal(D))
        @test JosephsonCircuits.is_symplectic_pair(Q)

    end

    @testset "pre_iwasawa_block" begin
        # real
        S = JosephsonCircuits.rand_symplectic_block(Float64, 4)
        E, D, F = JosephsonCircuits.pre_iwasawa_block(S)
        @test isapprox(S, E * D * F)
        @test JosephsonCircuits.is_symplectic_block(F)

        # complex
        S = JosephsonCircuits.rand_symplectic_block(Complex{Float64}, 4)
        E, D, F = JosephsonCircuits.pre_iwasawa_block(S)
        @test isapprox(S, E * D * F)
        @test JosephsonCircuits.is_symplectic_block(F)

        @test_throws(
            ErrorException(lazy"A must be symplectic."),
            JosephsonCircuits.pre_iwasawa_block(Float64[1 1;-1 1]),
        )
    end

    @testset "pre_iwasawa_pair" begin
        # real
        S = JosephsonCircuits.rand_symplectic_pair(Float64, 4)
        E, D, F = JosephsonCircuits.pre_iwasawa_pair(S)
        @test isapprox(S, E * D * F)
        @test JosephsonCircuits.is_symplectic_pair(F)

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

        # complex
        S = JosephsonCircuits.rand_symplectic_block(Complex{Float64}, 2)
        F = JosephsonCircuits.iwasawa_block(S)
        @test isapprox(F.K * F.K', I(4))
        @test isapprox(F.K * F.A * F.N, S)
        @test JosephsonCircuits.is_symplectic_block(F.K)
        @test JosephsonCircuits.is_symplectic_block(F.A)
        @test JosephsonCircuits.is_symplectic_block(F.N)

        @test_throws(
            ErrorException(lazy"A must be symplectic."),
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

        # singular example
        # this test seems brittle. think about how to replace it.
        vals, vecs = eigen(Aa)
        Aa1 = real(vecs * Diagonal([0, 0, vals[3], vals[4]]) * vecs')
        Q1 = JosephsonCircuits.symplectic_normal_form_pair(Aa1)
        @test isapprox(Aa1, Q1 * Omega * Q1')

        # errors
        @test_throws(
            ErrorException(lazy"A must be skew-symmetric."),
            JosephsonCircuits.symplectic_normal_form_pair([1 1;1 1]),
        )
        @test_throws(
            ErrorException(lazy"A must have even dimensions for a symplectic normal form."),
            JosephsonCircuits.symplectic_normal_form_pair([0 0 1;0 0 0;-1 0 0]),
        )

    end

    @testset "symplectic_normal_form_block" begin

        A = randn(Float64, 4, 4)
        Aa = (A - A') / 2
        Q = JosephsonCircuits.symplectic_normal_form_block(Aa)
        Omega = JosephsonCircuits.symplectic_form_block(2)
        @test isapprox(Aa, Q * Omega * Q')

        # # singular example
        # vals, vecs = eigen(Aa)
        # Aa1 = real(vecs * Diagonal([0, 0, vals[3], vals[4]]) * vecs')
        # Q1 = JosephsonCircuits.symplectic_normal_form_pair(Aa1)
        # @test isapprox(Aa1, Q1 * Omega * Q1')

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

    end

    @testset "Ymin_from_X_quadrature_pair and Ymin_from_X_quadrature_block" begin

        X = rand(Float64,4,4)
        for method in 1:3
            @test JosephsonCircuits.is_cptp_quadrature_pair(X,JosephsonCircuits.Ymin_from_X_quadrature_pair(X;method=method))
            @test JosephsonCircuits.is_cptp_quadrature_block(X,JosephsonCircuits.Ymin_from_X_quadrature_block(X;method=method))
        end

        @test_throws(
            ErrorException(lazy"Unknown method"),
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
            ErrorException(lazy"`Y` must be positive semi-definite."),
            JosephsonCircuits.B_from_X_Y_quadrature_block([1 0;0 1],[1 0;0 -1]),
        )
    end

    @testset "X_Y_to_sympletic_pair" begin

        X = rand(Float64,4,4)
        Y = JosephsonCircuits.Ymin_from_X_quadrature_pair(X)
        S = JosephsonCircuits.X_Y_to_sympletic_pair(X,Y)
        @test JosephsonCircuits.is_symplectic_pair(S)

        X, Y = JosephsonCircuits.rand_cptp_quadrature_pair(2)
        S = JosephsonCircuits.X_Y_to_sympletic_pair(X,Y)
        @test JosephsonCircuits.is_symplectic_pair(S)

    end

    @testset "X_Y_to_sympletic_block" begin

        X = rand(Float64,4,4)
        Y = JosephsonCircuits.Ymin_from_X_quadrature_block(X)
        S = JosephsonCircuits.X_Y_to_sympletic_block(X,Y)
        @test JosephsonCircuits.is_symplectic_block(S)

        X, Y = JosephsonCircuits.rand_cptp_quadrature_block(2)
        S = JosephsonCircuits.X_Y_to_sympletic_block(X,Y)
        @test JosephsonCircuits.is_symplectic_block(S)

    end

    @testset "X_Y_to_bogoliubov_pair" begin

        X, Y = JosephsonCircuits.rand_cptp_ladder_pair(2)
        S = JosephsonCircuits.X_Y_to_bogoliubov_pair(X,Y)
        @test JosephsonCircuits.is_bogoliubov_pair(S)

    end

    @testset "X_Y_to_bogoliubov_block" begin

        X, Y = JosephsonCircuits.rand_cptp_ladder_block(2)
        S = JosephsonCircuits.X_Y_to_bogoliubov_block(X,Y)
        @test JosephsonCircuits.is_bogoliubov_block(S)

    end

    @testset "wmatrix" begin
        # the in-place form fills and returns a matrix of the size it checks
        w = zeros(2, 3)
        @test JosephsonCircuits.wmatrix!(w, 0.1:0.1:0.3, (1.0,), [(1,), (-1,)]) === w
        @test w == JosephsonCircuits.wmatrix(0.1:0.1:0.3, (1.0,), [(1,), (-1,)])
        @test_throws(DimensionMismatch,
            JosephsonCircuits.wmatrix!(zeros(3, 3), 0.1:0.1:0.3, (1.0,), [(1,), (-1,)]))
    end

    @testset "interpolate_scattering" begin

        # test with extrapolation
        w = 0.01:0.01:1.0
        S = JosephsonCircuits.ABCD_tline(50,w)
        @test isapprox(S,JosephsonCircuits.interpolate_scattering(w,S,w;extrap=true))

        # test with negative frequencies
        w = 0.01:0.01:1.0
        S = JosephsonCircuits.ABCD_tline(50,w)
        @test isapprox(S,JosephsonCircuits.interpolate_scattering(w,conj.(S),-w))

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
            ErrorException(lazy"`S` must have 3 dimensions. The first two are ports and the third is frequencies."),
            JosephsonCircuits.interpolate_scattering(w,randn(Complex{Float64},2),w),
        )

        @test_throws(
            ErrorException(lazy"The length of the third dimension of `S` must be equal to the number of frequencies."),
            JosephsonCircuits.interpolate_scattering(w,randn(Complex{Float64},2,2,2*length(w)),w),
        )
    end

end
