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

    @testset "spaddkeepzeros" begin
        A = JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [1,2,0],2,2);
        B = JosephsonCircuits.SparseArrays.sparse([1,2], [1,2], [1,1],3,2);
        @test_throws(
            DimensionMismatch("argument shapes must match"),
            JosephsonCircuits.spaddkeepzeros(A,B)
        )
    end

    @testset "sparseadd!" begin
        A = JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [1,2,-3],2,2)
        As = JosephsonCircuits.SparseArrays.sparse([1,1], [1,2], [3,4],2,2)
        indexmap = JosephsonCircuits.sparseaddmap(A,As)
        @test_throws(
            DimensionMismatch("As cannot have more nonzero elements than A"),
            JosephsonCircuits.sparseadd!(As,A,indexmap)
        )
    end

    @testset "sparseadd!" begin
        A = JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [1,2,-3],2,2)
        As = JosephsonCircuits.SparseArrays.sparse([1,1], [1,2], [3,4],2,2)
        indexmap = JosephsonCircuits.sparseaddmap(A,As)
        @test_throws(
            DimensionMismatch("The indexmap must be the same length as As"),
            JosephsonCircuits.sparseadd!(A,As,indexmap[1:end-1])
        )
    end

    @testset "sparseadd!" begin
        A = JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [1,2,-3],2,2)
        As = JosephsonCircuits.SparseArrays.sparse([1,1], [1,2], [3,4],2,2)
        indexmap = JosephsonCircuits.sparseaddmap(A,As)
        As = JosephsonCircuits.SparseArrays.sparse([1,1], [1,2], [3,4],3,3)
        @test_throws(
            DimensionMismatch("A and As must be the same size."),
            JosephsonCircuits.sparseadd!(A,As,indexmap)
        )
    end

    @testset "sparseadd!" begin
        A = JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [1,2,-3],2,2)
        As = JosephsonCircuits.SparseArrays.sparse([1,1], [1,2], [3,4],2,2)
        As2 = JosephsonCircuits.SparseArrays.sparse([1,1], [1,2], [3,4],3,3)
        indexmap = JosephsonCircuits.sparseaddmap(A,As)
        @test_throws(
            DimensionMismatch("A and As must be the same size."),
            JosephsonCircuits.sparseadd!(A,2,As2,indexmap)
        )
    end

    @testset "sparseadd!" begin
        A = JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [1,2,-3],2,2)
        As = JosephsonCircuits.SparseArrays.sparse([1,1], [1,2], [3,4],2,2)
        indexmap = JosephsonCircuits.sparseaddmap(A,As)
        @test_throws(
            DimensionMismatch("The indexmap must be the same length as As"),
            JosephsonCircuits.sparseadd!(A,2,As,indexmap[1:end-1])
        )
    end

    @testset "sparseadd!" begin
        A = JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [1,2,-3],2,2)
        As = JosephsonCircuits.SparseArrays.sparse([1,1], [1,2], [3,4],2,2)
        indexmap = JosephsonCircuits.sparseaddmap(A,As)
        @test_throws(
            DimensionMismatch("As cannot have more nonzero elements than A"),
            JosephsonCircuits.sparseadd!(As,2,A,indexmap)
        )
    end
 
    @testset "sparseadd!" begin
        A = JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [1,2,-3],2,2)
        As = JosephsonCircuits.SparseArrays.sparse([1,1], [1,2], [3,4],2,2)
        Ad = JosephsonCircuits.LinearAlgebra.Diagonal([1,-2])
        indexmap = JosephsonCircuits.sparseaddmap(A,As)
        @test_throws(
            DimensionMismatch("As cannot have more nonzero elements than A"),
            JosephsonCircuits.sparseadd!(As,2,A,Ad,indexmap)
        )
    end

    @testset "sparseadd!" begin
        A = JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [1,2,-3],2,2)
        As = JosephsonCircuits.SparseArrays.sparse([1,1], [1,2], [3,4],2,2)
        Ad = JosephsonCircuits.LinearAlgebra.Diagonal([1,-2])
        indexmap = JosephsonCircuits.sparseaddmap(A,As)
        @test_throws(
            DimensionMismatch("The indexmap must be the same length as As"),
            JosephsonCircuits.sparseadd!(A,2,As,Ad,indexmap[1:end-1])
        )
    end

    @testset "sparseadd!" begin
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

    @testset "sparseadd!" begin
        A = JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [1,2,-3],2,2)
        As = JosephsonCircuits.SparseArrays.sparse([1,1], [1,2], [3,4],2,2)
        Ad = JosephsonCircuits.LinearAlgebra.Diagonal([1,-2,1])
        indexmap = JosephsonCircuits.sparseaddmap(A,As)
        @test_throws(
            DimensionMismatch("A and Ad must be the same size."),
            JosephsonCircuits.sparseadd!(A,2,As,Ad,indexmap)
        )
    end

    @testset "sparseadd!" begin
        A = JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [1,2,-3],2,2)
        As = JosephsonCircuits.SparseArrays.sparse([1,1], [1,2], [3,4],2,2)
        Ad = JosephsonCircuits.LinearAlgebra.Diagonal([1,-2])
        indexmap = JosephsonCircuits.sparseaddmap(A,As)
        @test_throws(
            DimensionMismatch("As cannot have more nonzero elements than A"),
            JosephsonCircuits.sparseadd!(As,2,Ad,A,indexmap)
        )
    end

    @testset "sparseadd!" begin
        A = JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [1,2,-3],2,2)
        As = JosephsonCircuits.SparseArrays.sparse([1,1], [1,2], [3,4],2,2)
        Ad = JosephsonCircuits.LinearAlgebra.Diagonal([1,-2])
        indexmap = JosephsonCircuits.sparseaddmap(A,As)
        @test_throws(
            DimensionMismatch("The indexmap must be the same length as As"),
            JosephsonCircuits.sparseadd!(A,2,Ad,As,indexmap[1:end-1])
        )
    end

    @testset "sparseadd!" begin
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

    @testset "sparseadd!" begin
        A = JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [1,2,-3],2,2)
        As = JosephsonCircuits.SparseArrays.sparse([1,1], [1,2], [3,4],2,2)
        Ad = JosephsonCircuits.LinearAlgebra.Diagonal([1,-2,1])
        indexmap = JosephsonCircuits.sparseaddmap(A,As)
        @test_throws(
            DimensionMismatch("A and Ad must be the same size."),
            JosephsonCircuits.sparseadd!(A,2,Ad,As,indexmap)
        )
    end

    @testset "sparseaddconjsubst!" begin
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

    @testset "sparseaddconjsubst!" begin
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

    @testset "sparseaddconjsubst!" begin
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

    @testset "sparseaddconjsubst!" begin
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

    @testset "sparseaddconjsubst!" begin
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

    @testset "sparseaddconjsubst!" begin
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

    @testset "sparseaddmap" begin
        As = JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [1,2,-3],2,2)
        A = JosephsonCircuits.SparseArrays.sparse([1,2], [1,2], [4,2],2,2)
        @test_throws(
            ErrorException("Coordinate not found. Are the positions of elements in As a subset of the positions of elements in A?"),
            JosephsonCircuits.sparseaddmap(A,As)
        )
    end

    @testset "sparseaddmap" begin
        A = JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [1,2,-3],4,4)
        As = JosephsonCircuits.SparseArrays.sparse([1,2], [1,2], [4,2],2,2)
        @test_throws(
            DimensionMismatch("A and B must be the same size."),
            JosephsonCircuits.sparseaddmap(A,As)
        )
    end

    @testset "conjnegfreq!" begin
        A = JosephsonCircuits.SparseArrays.sparse([1,2,1,2], [1,1,2,2], [1+1im,1+1im,1+1im,1+1im],2,2);
        @test_throws(
            DimensionMismatch("The dimensions of A must be integer multiples of the length of wmodes."),
            JosephsonCircuits.conjnegfreq!(A,[-1,1,1])
        )
    end

    @testset "freqsubst" begin
        @variables w
        wmodes = [-1,2];
        A = JosephsonCircuits.diagrepeat(JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [w,2*w,3*w],2,2),2);
        @test_throws(
            ErrorException("Set symfreqvar equal to the symbolic variable representing frequency."),
            JosephsonCircuits.freqsubst(A,wmodes,nothing)
        )
    end

    @testset "freqsubst" begin
        @variables w
        wmodes = [-1,1,2];
        A = JosephsonCircuits.diagrepeat(JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [w,2*w,3*w],2,2),2);
        @test_throws(
            DimensionMismatch("The dimensions of A must be integer multiples of the length of wmodes."),
            JosephsonCircuits.freqsubst(A,wmodes,w)
        )
    end

    @testset "freqsubst" begin
        wmodes = [-1,2];
        A = JosephsonCircuits.diagrepeat(JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [1,2,3],2,2),2);
        @test_throws(
            ErrorException("symfreqvar must be a symbolic variable (or nothing if no symbolic variables)"),
            JosephsonCircuits.freqsubst(A,wmodes,1)
        )
    end

    @testset "spmatmul!" begin
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

    @testset "spmatmul!" begin
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

    @testset "spmatmul!" begin
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

    @testset "spmatmul!" begin
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