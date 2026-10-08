using Test, JosephsonCircuits, LinearAlgebra, SparseArrays

@testset "Local pole response assembly" begin
    JC = JosephsonCircuits

    @testset "coordinate gathering and common-voltage cancellation" begin
        F = sparse([0.0 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0])
        V = sparse([1.0 0 0 0 0 0; 1 0 0 0 0 0; 1 0 0 0 0 0])
        incidences = [[(1,1.0),(2,-1.0)], [(2,0.5),(3,-0.5)], [(1,1.0),(1,-1.0)]]
        cols, localF, localV = JC.polelocalmaps(sparse(transpose(F)), sparse(transpose(V)), incidences, [5])
        E = [1.0 -1 0; 0 0.5 -0.5; 0 0 0]
        @test cols == [1,2,3,5]
        @test localF == (E*F)[:,cols]
        @test localV == (E*V)[:,cols] == zeros(3,4)
    end

    @testset "overlapping stamps and frequency-dependent zeros" begin
        n = 7
        response(s) = [s-1 0.3im; exp(-s) 2+s]
        terms = [JC.PoleResponseTerm(response, 0.0, [2,5], [1,4,7],
            ComplexF64[1 0.2; -0.5 2], ComplexF64[1 0 2; 0 -1 0],
            ComplexF64[0 0.5 0; 1 0 0.2]),
            JC.PoleResponseTerm(s -> fill(s^2,1,1), 0.4, [2], [4],
                fill(2.0+0im,1,1), fill(1.0+0im,1,1), fill(-1.0+0im,1,1))]
        sys = (; Q0 = spdiagm(0 => ComplexF64.(1:n)),
            Q1 = sparse([1],[7],[0.5+0im],n,n), Q2 = spzeros(ComplexF64,n,n),
            terms, scale = 1.0)
        work = JC.PoleMatrixWorkspace(sys)
        for z in (1.0+0im, -0.7+0.4im, 0.0+0im, 2.0-0.3im)
            expected = Matrix(sys.Q0 + z*sys.Q1 + z^2*sys.Q2)
            for term in terms
                # Independent dense embedding is deliberately confined to
                # this small reference: it exercises all local row/col maps.
                L, R0, R1 = zeros(ComplexF64,n,size(term.left,2)),
                    zeros(ComplexF64,size(term.right0,1),n), zeros(ComplexF64,size(term.right1,1),n)
                L[term.rows,:] = term.left
                R0[:,term.cols], R1[:,term.cols] = term.right0, term.right1
                expected += L*term.response(z+im*term.offset)*(R0+z*R1)
            end
            @test JC.polematrix!(work,sys,z) ≈ expected
        end
        bad = JC.PoleResponseTerm(s -> zeros(3,3), 0.0, [1], [1],
            ones(ComplexF64,1,1), ones(ComplexF64,1,1), zeros(ComplexF64,1,1))
        @test_throws DimensionMismatch JC.poleterm(bad,1.0,1.0)
        bad = JC.PoleResponseTerm(s -> fill(Inf,1,1), 0.0, [1], [1],
            ones(ComplexF64,1,1), ones(ComplexF64,1,1), zeros(ComplexF64,1,1))
        @test_throws ArgumentError JC.poleterm(bad,1.0,1.0)
    end
end
