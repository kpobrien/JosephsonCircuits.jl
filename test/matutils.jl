using JosephsonCircuits
using Test

@testset verbose=true "matutils" begin

    @testset "diagrepeat!" begin
        A = [1 2;3 4]
        out = zeros(eltype(A),4,4)
        @test_throws(
            DimensionMismatch("Sizes not consistent"),
            JosephsonCircuits.diagrepeat!(out,A,1)
        )
    end

    @testset "diagcombine" begin
        @test_throws(
            DimensionMismatch("Sizes are not consistent."),
            JosephsonCircuits.diagcombine([[111 121;211 221],[112 122;212 222],[113 123]])
        )
    end

    @testset "diagcombine!" begin
        @test_throws(
            ArgumentError("`mode_index` = 0 must be greater than zero."),
            JosephsonCircuits.diagcombine!(zeros(4,4),zeros(2,2),0)
        )
        @test_throws(
            ArgumentError("`mode_index` = 10 must be less than or equal to the number of modes, which is 2 from the matrix sizes."),
            JosephsonCircuits.diagcombine!(zeros(4,4),zeros(2,2),10)
        )
    end

    @testset "axis_to_modes" begin
        @test_throws(
            DimensionMismatch("The input array needs 3 or more dimensions (two for ports and one for modes)."),
            JosephsonCircuits.axis_to_modes([111 122],3)
        )
        @test_throws(
            ArgumentError("`modes_axis` must be 3 or more (the first two dimensions are ports."),
            JosephsonCircuits.axis_to_modes([111 121;211 221;;; 112 122;212 222;;; 113 123;213 223],0)
        )
        @test_throws(
            ArgumentError("`modes_axis` must be less than or equal to the number of dimensions in input array."),
            JosephsonCircuits.axis_to_modes([111 121;211 221;;; 112 122;212 222;;; 113 123;213 223],4)
        )
    end

    @testset "axis_to_modes!" begin
        @test_throws(
            DimensionMismatch("The S parameter array `S` must have 3 dimensions (the first two dimensions are ports and the last is the modes)."),
            JosephsonCircuits.axis_to_modes!([1.0 1.0;1.0 1.0],[1.0,1.0])
        )
        @test_throws(
            DimensionMismatch("The first two dimensions of `out` must equal the first two dimensions of `S` times `Nmodes`."),
            JosephsonCircuits.axis_to_modes!(zeros(4,5),zeros(2,2,2))
        )
    end

    @testset "spaddkeepzeros" begin
        A = JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [1,2,0],2,2);
        B = JosephsonCircuits.SparseArrays.sparse([1,2], [1,2], [1,1],3,2);
        @test_throws(
            DimensionMismatch("argument shapes must match"),
            JosephsonCircuits.spaddkeepzeros(A,B)
        )
    end

    @testset "sparseadd!" begin
        begin
            A = JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [1,2,-3],2,2)
            As = JosephsonCircuits.SparseArrays.sparse([1,1], [1,2], [3,4],2,2)
            indexmap = JosephsonCircuits.sparseaddmap(A,As)
            @test_throws(
                DimensionMismatch("As cannot have more nonzero elements than A"),
                JosephsonCircuits.sparseadd!(As,A,indexmap)
            )
        end

        begin
            A = JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [1,2,-3],2,2)
            As = JosephsonCircuits.SparseArrays.sparse([1,1], [1,2], [3,4],2,2)
            indexmap = JosephsonCircuits.sparseaddmap(A,As)
            @test_throws(
                DimensionMismatch("The indexmap must be the same length as As"),
                JosephsonCircuits.sparseadd!(A,As,indexmap[1:end-1])
            )
        end

        begin
            A = JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [1,2,-3],2,2)
            As = JosephsonCircuits.SparseArrays.sparse([1,1], [1,2], [3,4],2,2)
            indexmap = JosephsonCircuits.sparseaddmap(A,As)
            As = JosephsonCircuits.SparseArrays.sparse([1,1], [1,2], [3,4],3,3)
            @test_throws(
                DimensionMismatch("A and As must be the same size."),
                JosephsonCircuits.sparseadd!(A,As,indexmap)
            )
        end

        begin
            A = JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [1,2,-3],2,2)
            As = JosephsonCircuits.SparseArrays.sparse([1,1], [1,2], [3,4],2,2)
            As2 = JosephsonCircuits.SparseArrays.sparse([1,1], [1,2], [3,4],3,3)
            indexmap = JosephsonCircuits.sparseaddmap(A,As)
            @test_throws(
                DimensionMismatch("A and As must be the same size."),
                JosephsonCircuits.sparseadd!(A,2,As2,indexmap)
            )
        end

        begin
            A = JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [1,2,-3],2,2)
            As = JosephsonCircuits.SparseArrays.sparse([1,1], [1,2], [3,4],2,2)
            indexmap = JosephsonCircuits.sparseaddmap(A,As)
            @test_throws(
                DimensionMismatch("The indexmap must be the same length as As"),
                JosephsonCircuits.sparseadd!(A,2,As,indexmap[1:end-1])
            )
        end

        begin
            A = JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [1,2,-3],2,2)
            As = JosephsonCircuits.SparseArrays.sparse([1,1], [1,2], [3,4],2,2)
            indexmap = JosephsonCircuits.sparseaddmap(A,As)
            @test_throws(
                DimensionMismatch("As cannot have more nonzero elements than A"),
                JosephsonCircuits.sparseadd!(As,2,A,indexmap)
            )
        end
     
        begin
            A = JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [1,2,-3],2,2)
            As = JosephsonCircuits.SparseArrays.sparse([1,1], [1,2], [3,4],2,2)
            Ad = JosephsonCircuits.LinearAlgebra.Diagonal([1,-2])
            indexmap = JosephsonCircuits.sparseaddmap(A,As)
            @test_throws(
                DimensionMismatch("As cannot have more nonzero elements than A"),
                JosephsonCircuits.sparseadd!(As,2,A,Ad,indexmap)
            )
        end

        begin
            A = JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [1,2,-3],2,2)
            As = JosephsonCircuits.SparseArrays.sparse([1,1], [1,2], [3,4],2,2)
            Ad = JosephsonCircuits.LinearAlgebra.Diagonal([1,-2])
            indexmap = JosephsonCircuits.sparseaddmap(A,As)
            @test_throws(
                DimensionMismatch("The indexmap must be the same length as As"),
                JosephsonCircuits.sparseadd!(A,2,As,Ad,indexmap[1:end-1])
            )
        end

        begin
            A = JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [1,2,-3],2,2)
            As = JosephsonCircuits.SparseArrays.sparse([1,1], [1,2], [3,4],2,2)
            Ad = JosephsonCircuits.LinearAlgebra.Diagonal([1,-2])
            indexmap = JosephsonCircuits.sparseaddmap(A,As)
            As = JosephsonCircuits.SparseArrays.sparse([1,1], [1,2], [3,4],3,3)
            @test_throws(
                DimensionMismatch("A and As must be the same size."),
                JosephsonCircuits.sparseadd!(A,2,As,Ad,indexmap)
            )
        end

        begin
            A = JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [1,2,-3],2,2)
            As = JosephsonCircuits.SparseArrays.sparse([1,1], [1,2], [3,4],2,2)
            Ad = JosephsonCircuits.LinearAlgebra.Diagonal([1,-2,1])
            indexmap = JosephsonCircuits.sparseaddmap(A,As)
            @test_throws(
                DimensionMismatch("A and Ad must be the same size."),
                JosephsonCircuits.sparseadd!(A,2,As,Ad,indexmap)
            )
        end

        begin
            A = JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [1,2,-3],2,2)
            As = JosephsonCircuits.SparseArrays.sparse([1,1], [1,2], [3,4],2,2)
            Ad = JosephsonCircuits.LinearAlgebra.Diagonal([1,-2])
            indexmap = JosephsonCircuits.sparseaddmap(A,As)
            @test_throws(
                DimensionMismatch("As cannot have more nonzero elements than A"),
                JosephsonCircuits.sparseadd!(As,2,Ad,A,indexmap)
            )
        end

        begin
            A = JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [1,2,-3],2,2)
            As = JosephsonCircuits.SparseArrays.sparse([1,1], [1,2], [3,4],2,2)
            Ad = JosephsonCircuits.LinearAlgebra.Diagonal([1,-2])
            indexmap = JosephsonCircuits.sparseaddmap(A,As)
            @test_throws(
                DimensionMismatch("The indexmap must be the same length as As"),
                JosephsonCircuits.sparseadd!(A,2,Ad,As,indexmap[1:end-1])
            )
        end

        begin
            A = JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [1,2,-3],2,2)
            As = JosephsonCircuits.SparseArrays.sparse([1,1], [1,2], [3,4],2,2)
            Ad = JosephsonCircuits.LinearAlgebra.Diagonal([1,-2])
            indexmap = JosephsonCircuits.sparseaddmap(A,As)
            As = JosephsonCircuits.SparseArrays.sparse([1,1], [1,2], [3,4],3,3)
            @test_throws(
                DimensionMismatch("A and As must be the same size."),
                JosephsonCircuits.sparseadd!(A,2,Ad,As,indexmap)
            )
        end

        begin
            A = JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [1,2,-3],2,2)
            As = JosephsonCircuits.SparseArrays.sparse([1,1], [1,2], [3,4],2,2)
            Ad = JosephsonCircuits.LinearAlgebra.Diagonal([1,-2,1])
            indexmap = JosephsonCircuits.sparseaddmap(A,As)
            @test_throws(
                DimensionMismatch("A and Ad must be the same size."),
                JosephsonCircuits.sparseadd!(A,2,Ad,As,indexmap)
            )
        end
    end

    @testset "sparseaddconjsubst!" begin
        begin
            A = JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [1.0+1.0im,2.0+1.0im,-3.0+0.0im],2,2)
            Ad = JosephsonCircuits.LinearAlgebra.Diagonal([1,-2])
            As = JosephsonCircuits.SparseArrays.sparse([1,1], [1,2], [3.0+2.0im,4.0+3.0im],2,2)
            wmodesm = JosephsonCircuits.LinearAlgebra.Diagonal([-1,1,2])
            indexmap = JosephsonCircuits.sparseaddmap(A,As)
            freqsubstindices  = JosephsonCircuits.symbolicindices(As)
            @test_throws(
                DimensionMismatch("A and conjflag must be the same size."),
                JosephsonCircuits.sparseaddconjsubst!(A,2,As,Ad,indexmap,wmodesm .< 0,
                    wmodesm,freqsubstindices,nothing)
            )
        end

        begin
            A = JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [1.0+1.0im,2.0+1.0im,-3.0+0.0im],2,2)
            Ad = JosephsonCircuits.LinearAlgebra.Diagonal([1,-2])
            As = JosephsonCircuits.SparseArrays.sparse([1,1], [1,2], [3.0+2.0im,4.0+3.0im],2,2)
            wmodesm = JosephsonCircuits.LinearAlgebra.Diagonal([-1,1])
            indexmap = JosephsonCircuits.sparseaddmap(A,As)
            freqsubstindices  = JosephsonCircuits.symbolicindices(As)
            @test_throws(
                DimensionMismatch("As cannot have more nonzero elements than A"),
                JosephsonCircuits.sparseaddconjsubst!(As,2,A,Ad,indexmap,wmodesm .< 0,
                    wmodesm,freqsubstindices,nothing)
            )
        end

        begin
            A = JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [1.0+1.0im,2.0+1.0im,-3.0+0.0im],2,2)
            Ad = JosephsonCircuits.LinearAlgebra.Diagonal([1,-2])
            As = JosephsonCircuits.SparseArrays.sparse([1,1], [1,2], [3.0+2.0im,4.0+3.0im],2,2)
            wmodesm = JosephsonCircuits.LinearAlgebra.Diagonal([-1,1])
            indexmap = JosephsonCircuits.sparseaddmap(A,As)
            freqsubstindices  = JosephsonCircuits.symbolicindices(As)
            @test_throws(
                DimensionMismatch("The indexmap must be the same length as As"),
                JosephsonCircuits.sparseaddconjsubst!(A,2,As,Ad,indexmap[1:end-1],
                    wmodesm .< 0,wmodesm,freqsubstindices,nothing)
            )
        end

        begin
            A = JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [1.0+1.0im,2.0+1.0im,-3.0+0.0im],2,2)
            Ad = JosephsonCircuits.LinearAlgebra.Diagonal([1,-2])
            As = JosephsonCircuits.SparseArrays.sparse([1,1], [1,2], [3.0+2.0im,4.0+3.0im],2,2)
            wmodesm = JosephsonCircuits.LinearAlgebra.Diagonal([-1,1])
            wmodesm2 = JosephsonCircuits.LinearAlgebra.Diagonal([-1,1,2])
            indexmap = JosephsonCircuits.sparseaddmap(A,As)
            freqsubstindices  = JosephsonCircuits.symbolicindices(As)
            @test_throws(
                DimensionMismatch("A and wmodesm must be the same size."),
                JosephsonCircuits.sparseaddconjsubst!(A,2,As,Ad,indexmap,wmodesm .< 0,
                    wmodesm2,freqsubstindices,nothing)
            )
        end

        begin
            A = JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [1.0+1.0im,2.0+1.0im,-3.0+0.0im],2,2)
            Ad = JosephsonCircuits.LinearAlgebra.Diagonal([1,-2,1])
            As = JosephsonCircuits.SparseArrays.sparse([1,1], [1,2], [3.0+2.0im,4.0+3.0im],2,2)
            wmodesm = JosephsonCircuits.LinearAlgebra.Diagonal([-1,1])
            indexmap = JosephsonCircuits.sparseaddmap(A,As)
            freqsubstindices  = JosephsonCircuits.symbolicindices(As)
            @test_throws(
                DimensionMismatch("A and Ad must be the same size."),
                JosephsonCircuits.sparseaddconjsubst!(A,2,As,Ad,indexmap,wmodesm .< 0,
                    wmodesm,freqsubstindices,nothing)
            )
        end

        begin
            A = JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [1.0+1.0im,2.0+1.0im,-3.0+0.0im],2,2)
            Ad = JosephsonCircuits.LinearAlgebra.Diagonal([1,-2])
            As = JosephsonCircuits.SparseArrays.sparse([1,1], [1,2], [3.0+2.0im,4.0+3.0im],2,2)
            As2 = JosephsonCircuits.SparseArrays.sparse([1,1], [1,2], [3.0+2.0im,4.0+3.0im],3,3)
            wmodesm = JosephsonCircuits.LinearAlgebra.Diagonal([-1,1])
            indexmap = JosephsonCircuits.sparseaddmap(A,As)
            freqsubstindices  = JosephsonCircuits.symbolicindices(As)
            @test_throws(
                DimensionMismatch("A and As must be the same size."),
                JosephsonCircuits.sparseaddconjsubst!(A,2,As2,Ad,indexmap,
                    wmodesm .< 0,wmodesm,freqsubstindices,nothing)
            )
        end
    end

    @testset "sparseaddmap" begin
        begin
            As = JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [1,2,-3],2,2)
            A = JosephsonCircuits.SparseArrays.sparse([1,2], [1,2], [4,2],2,2)
            @test_throws(
                ErrorException("Coordinate not found. Are the positions of elements in As a subset of the positions of elements in A?"),
                JosephsonCircuits.sparseaddmap(A,As)
            )
        end

        begin
            A = JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [1,2,-3],4,4)
            As = JosephsonCircuits.SparseArrays.sparse([1,2], [1,2], [4,2],2,2)
            @test_throws(
                DimensionMismatch("A and B must be the same size."),
                JosephsonCircuits.sparseaddmap(A,As)
            )
        end
    end

    @testset "conjnegfreq!" begin
        A = JosephsonCircuits.SparseArrays.sparse([1,2,1,2], [1,1,2,2], [1+1im,1+1im,1+1im,1+1im],2,2);
        @test_throws(
            DimensionMismatch("The dimensions of A must be integer multiples of the length of wmodes."),
            JosephsonCircuits.conjnegfreq!(A,[-1,1,1])
        )
    end

    @testset "freqsubst" begin
        begin
            @variables w
            wmodes = [-1,2];
            A = JosephsonCircuits.diagrepeat(JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [w,2*w,3*w],2,2),2);
            @test_throws(
                ErrorException("Set symfreqvar equal to the symbolic variable representing frequency."),
                JosephsonCircuits.freqsubst(A,wmodes,nothing)
            )
        end

        begin
            @variables w
            wmodes = [-1,1,2];
            A = JosephsonCircuits.diagrepeat(JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [w,2*w,3*w],2,2),2);
            @test_throws(
                DimensionMismatch("The dimensions of A must be integer multiples of the length of wmodes."),
                JosephsonCircuits.freqsubst(A,wmodes,w)
            )
        end

        begin
            wmodes = [-1,2];
            A = JosephsonCircuits.diagrepeat(JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [1,2,3],2,2),2);
            @test_throws(
                ErrorException("symfreqvar must be a symbolic variable (or nothing if no symbolic variables)"),
                JosephsonCircuits.freqsubst(A,wmodes,1)
            )
        end
    end

    @testset "spmatmul!" begin
        begin
            a = JosephsonCircuits.sprand(100,100,0.1);
            b = JosephsonCircuits.sprand(100,100,0.1);
            c = a*b;
            d = copy(c)
            xb = fill(false, size(a,1));
            @test_throws(
                DimensionMismatch("Number of columns in A must equal number of rows in B."),
                JosephsonCircuits.spmatmul!(c,a[:,1:end-1],b,xb)
            )
        end

        begin
            a = JosephsonCircuits.sprand(100,100,0.1)
            b = JosephsonCircuits.sprand(100,100,0.1)
            c = a*b
            d = copy(c)
            xb = fill(false, size(a,1))
            @test_throws(
                DimensionMismatch("Number of rows in C must equal number of rows in A."),
                JosephsonCircuits.spmatmul!(c,a[1:end-1,:],b,xb)
            )
        end

        begin
            a = JosephsonCircuits.sprand(100,100,0.1)
            b = JosephsonCircuits.sprand(100,100,0.1)
            c = a*b
            d = copy(c)
            xb = fill(false, size(a,1))
            @test_throws(
                DimensionMismatch("Length of xb vector must equal number of rows in A."),
                JosephsonCircuits.spmatmul!(c,a,b,xb[1:end-1])
            )
        end

        begin
            a = JosephsonCircuits.sprand(100,100,0.1)
            b = JosephsonCircuits.sprand(100,100,0.1)
            c = a*b
            d = copy(c)
            xb = fill(false, size(a,1))
            @test_throws(
                DimensionMismatch("Number of columns in C must equal number of columns in B."),
                JosephsonCircuits.spmatmul!(c[:,1:end-1],a,b,xb)
            )
        end
    end

    @testset "lu_2x2" begin
        A = StaticArrays.SMatrix{2,2}(rand(Complex{Float64},2,2))
        Afact1 = LinearAlgebra.lu(A)
        Afact2 = JosephsonCircuits.lu_2x2(A)
        @test isapprox(Afact1.L,Afact2.L)
        @test isapprox(Afact1.U,Afact2.U)
        @test isapprox(Afact1.p,Afact2.p)

        A = StaticArrays.SMatrix{2,2}(rand(Complex{Float64},2,2))
        A = StaticArrays.SMatrix{2,2}(A[1,1],0*A[2,1],A[1,2],A[2,2])
        Afact1 = LinearAlgebra.lu(A)
        Afact2 = JosephsonCircuits.lu_2x2(A)
        @test isapprox(Afact1.L,Afact2.L)
        @test isapprox(Afact1.U,Afact2.U)
        @test isapprox(Afact1.p,Afact2.p)

        # LU decomposition of a matrix where A[1,1] = A[2,1] = 0
        A = rand(Complex{Float64},2,2)
        A[1,1] = 0
        A[2,1] = 0
        fact = JosephsonCircuits.lu_2x2(A);
        @test isapprox(fact.L*fact.U,A)
    end

    @testset "ldiv_2x2 errors" begin
        A = StaticArrays.SMatrix{2,2}(rand(Complex{Float64},2,2))
        b = StaticArrays.SVector{2}(rand(Complex{Float64},2))
        fact = JosephsonCircuits.lu_2x2(A)
        
        @test_throws(
            ArgumentError("Unknown pivot."),
            JosephsonCircuits.ldiv_2x2(StaticArrays.LU(fact.L,fact.U,
            StaticArrays.SVector{2}(3,4)),b)
        )

        # test the warning
        u11 = fact.U[1,1]
        u12 = fact.U[1,2]
        u22 = fact.U[2,2]*0
        U = LinearAlgebra.UpperTriangular(StaticArrays.SMatrix{2,2}(u11,zero(u11),u12,u22))
        fact2 = StaticArrays.LU(fact.L,U,fact.p)
        @test_throws(
            ArgumentError("Failed to solve linear system."),
            JosephsonCircuits.ldiv_2x2(fact2,b)
        )
    end

    @testset "ldiv_2x2" begin
        A = StaticArrays.SMatrix{2,2}(rand(Complex{Float64},2,2))
        b = StaticArrays.SVector{2}(rand(Complex{Float64},2))
        @test isapprox(
            JosephsonCircuits.ldiv_2x2(JosephsonCircuits.lu_2x2(A),b),
            LinearAlgebra.lu(A) \ b,
        )

        @test isapprox(
            JosephsonCircuits.ldiv_2x2(LinearAlgebra.lu(A),b),
            JosephsonCircuits.lu_2x2(A) \ b,
        )

        # set A21 equal to zero to test LU without pivoting.
        A = StaticArrays.SMatrix{2,2}(rand(Complex{Float64},2,2))
        B = StaticArrays.SMatrix{2,2}(A[1,1],0,A[1,2],A[2,2])
        @test isapprox(
            JosephsonCircuits.ldiv_2x2(JosephsonCircuits.lu_2x2(B),b),
            LinearAlgebra.lu(B) \ b,
        )

        @test isapprox(
            JosephsonCircuits.ldiv_2x2(LinearAlgebra.lu(B),b),
            JosephsonCircuits.lu_2x2(B) \ b,
        )
    end

end