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

        # @test JosephsonCircuits.is_symplectic_annihilation_creation_block(
        #     JosephsonCircuits.symplectic_form_annihilation_creation_block(4),
        # )

        # @test JosephsonCircuits.is_symplectic_annihilation_creation_pair(
        #     JosephsonCircuits.symplectic_form_annihilation_creation_pair(4),
        # )

        # Serafini B.3, the inverse equals the adjoint
        @test isapprox(
            inv(Matrix(JosephsonCircuits.symplectic_form_block(4))),
            adjoint(JosephsonCircuits.symplectic_form_block(4)),
        )

        @test isapprox(
            inv(Matrix(JosephsonCircuits.symplectic_form_pair(4))),
            adjoint(JosephsonCircuits.symplectic_form_pair(4)),
        )

        # @test isapprox(
        #     inv(Matrix(JosephsonCircuits.symplectic_form_annihilation_creation_block(4))),
        #     adjoint(JosephsonCircuits.symplectic_form_annihilation_creation_block(4)),
        # )

        # @test isapprox(
        #     inv(Matrix(JosephsonCircuits.symplectic_form_annihilation_creation_pair(4))),
        #     adjoint(JosephsonCircuits.symplectic_form_annihilation_creation_pair(4)),
        # )


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

        # @test isapprox(
        #     adjoint(JosephsonCircuits.symplectic_form_annihilation_creation_block(4)),
        #     -JosephsonCircuits.symplectic_form_annihilation_creation_block(4),
        # )

        # @test isapprox(
        #     adjoint(JosephsonCircuits.symplectic_form_annihilation_creation_pair(4)),
        #     -JosephsonCircuits.symplectic_form_annihilation_creation_pair(4),
        # )

        # Serafini pg. 31, symplectic form times itself is minus identity
        @test isapprox(
            JosephsonCircuits.symplectic_form_block(4)^2,
            -I(2 * 4),
        )

        @test isapprox(
            JosephsonCircuits.symplectic_form_pair(4)^2,
            -I(2 * 4),
        )

        # @test isapprox(
        #     JosephsonCircuits.symplectic_form_annihilation_creation_block(4)^2,
        #     -I(2 * 4),
        # )

        # @test isapprox(
        #     JosephsonCircuits.symplectic_form_annihilation_creation_pair(4)^2,
        #     -I(2 * 4),
        # )

        # Serafini pg. 31, symplectic form times adjoint is the identity
        @test isapprox(
            JosephsonCircuits.symplectic_form_block(4) * JosephsonCircuits.symplectic_form_block(4)',
            I(2 * 4),
        )

        @test isapprox(
            JosephsonCircuits.symplectic_form_pair(4) * JosephsonCircuits.symplectic_form_pair(4)',
            I(2 * 4),
        )

        # @test isapprox(
        #     JosephsonCircuits.symplectic_form_annihilation_creation_block(4) * JosephsonCircuits.symplectic_form_annihilation_creation_block(4)',
        #     I(2 * 4),
        # )

        # @test isapprox(
        #     JosephsonCircuits.symplectic_form_annihilation_creation_pair(4) * JosephsonCircuits.symplectic_form_annihilation_creation_pair(4)',
        #     I(2 * 4),
        # )


        # convert between the symplectic forms
        @test isapprox(
            JosephsonCircuits.symplectic_form_block(4),
            JosephsonCircuits.pair_to_block(JosephsonCircuits.symplectic_form_pair(4)),
        )

        @test isapprox(
            JosephsonCircuits.symplectic_form_block(4),
            JosephsonCircuits.pair_to_block2(JosephsonCircuits.symplectic_form_pair(4)),
        )

        # @test isapprox(
        #     JosephsonCircuits.symplectic_form_block(4),
        #     JosephsonCircuits.bogoliubov_to_quadrature_block(JosephsonCircuits.symplectic_form_annihilation_creation_block(4)),
        # )

        @test isapprox(
            JosephsonCircuits.symplectic_form_pair(4),
            JosephsonCircuits.block_to_pair(JosephsonCircuits.symplectic_form_block(4)),
        )

        @test isapprox(
            JosephsonCircuits.symplectic_form_pair(4),
            JosephsonCircuits.block_to_pair2(JosephsonCircuits.symplectic_form_block(4)),
        )

        # @test isapprox(
        #     JosephsonCircuits.symplectic_form_pair(4),
        #     JosephsonCircuits.bogoliubov_to_quadrature_pair(JosephsonCircuits.symplectic_form_annihilation_creation_pair(4)),
        # )

        # @test isapprox(
        #     JosephsonCircuits.quadrature_to_bogoliubov_block(JosephsonCircuits.symplectic_form_block(4)),
        #     JosephsonCircuits.symplectic_form_annihilation_creation_block(4),
        # )

        # @test isapprox(
        #     JosephsonCircuits.quadrature_to_ladder_pair(JosephsonCircuits.symplectic_form_pair(4)),
        #     JosephsonCircuits.symplectic_form_annihilation_creation_pair(4),
        # )

        # test the conversions and their inverses
        @test isapprox(
            JosephsonCircuits.symplectic_form_block(4),
            JosephsonCircuits.block_to_pair(JosephsonCircuits.pair_to_block(JosephsonCircuits.symplectic_form_block(4))),
        )

        # @test isapprox(
        #     JosephsonCircuits.symplectic_form_annihilation_creation_pair(4),
        #     JosephsonCircuits.quadrature_to_ladder_pair(JosephsonCircuits.bogoliubov_to_quadrature_pair(JosephsonCircuits.symplectic_form_annihilation_creation_pair(4))),
        # )

        # @test isapprox(
        #     JosephsonCircuits.symplectic_form_annihilation_creation_block(4),
        #     JosephsonCircuits.quadrature_to_bogoliubov_block(JosephsonCircuits.bogoliubov_to_quadrature_block(JosephsonCircuits.symplectic_form_annihilation_creation_block(4))),
        # )
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

        @test JosephsonCircuits.is_posdef_symplectic_pair(JosephsonCircuits.rand_posdef_symplectic_pair(4))
        @test JosephsonCircuits.is_posdef_symplectic_pair(JosephsonCircuits.rand_posdef_symplectic_pair(Float64,4))
        @test JosephsonCircuits.is_posdef_symplectic_pair(JosephsonCircuits.rand_posdef_symplectic_pair(Complex{Float64},4))

        @test JosephsonCircuits.is_posdef_symplectic_block(JosephsonCircuits.rand_posdef_symplectic_block(4))
        @test JosephsonCircuits.is_posdef_symplectic_block(JosephsonCircuits.rand_posdef_symplectic_block(Float64,4))
        @test JosephsonCircuits.is_posdef_symplectic_block(JosephsonCircuits.rand_posdef_symplectic_block(Complex{Float64},4))

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

    @testset "random cptp symplectic" begin

        # symplectic
        @test JosephsonCircuits.is_cptp_pair(JosephsonCircuits.rand_cptp_pair(4)...)
        @test JosephsonCircuits.is_cptp_pair(JosephsonCircuits.rand_cptp_pair(Float64,4)...)
        # @test JosephsonCircuits.is_cptp_pair(JosephsonCircuits.rand_cptp_pair(Complex{Float64},4)...)
        
        @test JosephsonCircuits.is_cptp_block(JosephsonCircuits.rand_cptp_block(4)...)
        @test JosephsonCircuits.is_cptp_block(JosephsonCircuits.rand_cptp_block(Float64,4)...)
        # @test JosephsonCircuits.is_cptp_block(JosephsonCircuits.rand_cptp_block(Complex{Float64},4)...)

    end

    @testset "random cptp bogoliubov" begin

        # bogoliubov
        @test JosephsonCircuits.is_cptp_bogoliubov_pair(JosephsonCircuits.rand_cptp_bogoliubov_pair(4)...)
        # @test JosephsonCircuits.is_cptp_bogoliubov_pair(JosephsonCircuits.rand_cptp_bogoliubov_pair(Float64,4)...)
        @test JosephsonCircuits.is_cptp_bogoliubov_pair(JosephsonCircuits.rand_cptp_bogoliubov_pair(Complex{Float64},4)...)
        
        @test JosephsonCircuits.is_cptp_bogoliubov_block(JosephsonCircuits.rand_cptp_bogoliubov_block(4)...)
        # @test JosephsonCircuits.is_cptp_bogoliubov_block(JosephsonCircuits.rand_cptp_bogoliubov_block(Float64,4)...)
        @test JosephsonCircuits.is_cptp_bogoliubov_block(JosephsonCircuits.rand_cptp_bogoliubov_block(Complex{Float64},4)...)

    end

    @testset "polar" begin

        A = randn(Complex{Float64}, 4, 4)
        P, Y = JosephsonCircuits.polar(A)
        @test isapprox(P * Y, A)
        @test isapprox(Y * Y', I(4))

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
        vals, vecs = eigen(Aa)
        Aa1 = real(vecs * Diagonal([0, 0, vals[3], vals[4]]) * vecs')
        Q1 = JosephsonCircuits.symplectic_normal_form_pair(Aa1)
        @test isapprox(Aa1, Q1 * Omega * Q1')

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
        B1 = JosephsonCircuits.B_from_X_Y_pair(X, Y)
        S1 = JosephsonCircuits.A_B_to_symplectic_pair(A, B1)
        @test JosephsonCircuits.is_symplectic_pair(S1)

    end

end
