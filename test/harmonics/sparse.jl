using JosephsonCircuits
using LinearAlgebra
using Test

@testset verbose=true "the sparse harmonic matrices" begin

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
     
    end

    @testset "sparseaddmap" begin
        begin
            As = JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [1,2,-3],2,2)
            A = JosephsonCircuits.SparseArrays.sparse([1,2], [1,2], [4,2],2,2)
            @test_throws ArgumentError JosephsonCircuits.sparseaddmap(A,As)
        end

        # the same map whatever the index type of the matrices
        begin
            A = JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [1,2,-3],2,2)
            As = JosephsonCircuits.SparseArrays.sparse([1,2], [1,2], [4,2],2,2)
            A32 = JosephsonCircuits.SparseArrays.sparse(Int32[1,2,1], Int32[1,2,2], [1,2,-3],2,2)
            As32 = JosephsonCircuits.SparseArrays.sparse(Int32[1,2], Int32[1,2], [4,2],2,2)
            @test JosephsonCircuits.sparseaddmap(A32,As32) ==
                JosephsonCircuits.sparseaddmap(A,As) == [1, 3]
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
            JosephsonCircuits.@params w
            wmodes = [-1,2];
            A = JosephsonCircuits.diagrepeat(JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [w,2*w,3*w],2,2),2);
            @test_throws(
                str -> occursin("FrequencyDependent closure", str),
                JosephsonCircuits.freqsubst(A,wmodes)
            )
        end

        begin
            wmodes = [-1,1,2];
            f = JosephsonCircuits.FrequencyDependent
            A = JosephsonCircuits.diagrepeat(JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [f(w->w),f(w->2*w),f(w->3*w)],2,2),2);
            @test_throws(
                DimensionMismatch("The dimensions of A must be integer multiples of the length of wmodes."),
                JosephsonCircuits.freqsubst(A,wmodes)
            )
        end

        # each frequency dependent value is resolved at the magnitude of
        # the mode frequency of its own column
        begin
            wmodes = [-1,2];
            f = JosephsonCircuits.FrequencyDependent
            A = JosephsonCircuits.diagrepeat(JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], Any[f(w->1.0*w),f(w->2.0*w),f(w->3.0*w)],2,2),2);
            B = JosephsonCircuits.freqsubst(A,wmodes)
            @test B[1,1] == 1.0 && B[2,2] == 2.0
            @test B[3,3] == 2.0 && B[4,4] == 4.0
            @test B[1,3] == 3.0 && B[2,4] == 6.0
        end
    end

end