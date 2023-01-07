using JosephsonCircuits
using Test

@testset verbose=true "exportnetlist" begin
    @testset "sumvalues" begin
        @test_throws(
            "unknown component type in sumvalues",
            JosephsonCircuits.sumvalues(:V, 1.0, 4.0)
        )
    end

    @testset "componentdictionaries" begin
        typevector = [:P, :I, :R, :L, :K, :K, :L, :C]
        nodeindexarraysorted = [2 2 2 2 0 0 3 3 3; 1 1 1 1 0 0 1 1 1]
        namedict = Dict("L1" => 4, "I1" => 2, "L2" => 7, "C2" => 8, "K2" => 6, "C3" => 9, "R1" => 3, "P1" => 1, "K1" => 5)
        mutualinductorvector = ["L1", "L2", "L1", "L2"]
        @test_throws(
            "DimensionMismatch: Input arrays must have the same length",
            JosephsonCircuits.componentdictionaries(typevector,
                nodeindexarraysorted,namedict,mutualinductorvector)
        )
    end

    @testset "componentdictionaries" begin
        typevector = [:P, :I, :R, :L, :K, :K, :L, :C, :C]
        nodeindexarraysorted = [2 2 2 2 0 0 3 3 3; 1 1 1 1 0 0 1 1 1; 1 1 1 1 0 0 1 1 1]
        namedict = Dict("L1" => 4, "I1" => 2, "L2" => 7, "C2" => 8, "K2" => 6, "C3" => 9, "R1" => 3, "P1" => 1, "K1" => 5)
        mutualinductorvector = ["L1", "L2", "L1", "L2"]
        @test_throws("DimensionMismatch: The length of the first axis must be 2",
            JosephsonCircuits.componentdictionaries(typevector,
                nodeindexarraysorted,namedict,mutualinductorvector)
        )
    end

    @testset "calcCjIcmean" begin
        typevector = [:P, :I, :R, :L, :K, :K, :L, :C]
        nodeindexarraysorted = [2 2 2 2 0 0 3 3 3; 1 1 1 1 0 0 1 1 1]
        namedict = Dict("L1" => 4, "I1" => 2, "L2" => 7, "C2" => 8, "K2" => 6, "C3" => 9, "R1" => 3, "P1" => 1, "K1" => 5)
        mutualinductorvector = ["L1", "L2", "L1", "L2"]
        @test_throws("DimensionMismatch: Input arrays must have the same length",
            JosephsonCircuits.componentdictionaries(typevector,
                nodeindexarraysorted,namedict,mutualinductorvector)
        )
    end

    @testset "calcCjIcmean" begin
        typevector = [:P, :R, :C, :Lj, :C, :C, :Lj, :C]
        nodeindexarray = [2 2 2 3 3 3 4 4; 1 1 3 1 1 4 1 1]
        valuevector = Real[1, 50.0, 1.0e-13, 1.0e-9, 1.0e-12, 1.0e-13, 100*1.1e-9, 1.2e-12]
        namedict = Dict("R1" => 2, "Cc2" => 6, "Cj2" => 8, "Cj1" => 5, "P1" => 1, "Cc1" => 3, "Lj2" => 7, "Lj1" => 4)
        mutualinductorvector = String[]
        countdict = Dict((:Lj, 1, 4) => 1, (:C, 3, 4) => 1, (:C, 1, 4) => 1, (:Lj, 1, 3) => 1, (:R, 1, 2) => 1, (:P, 1, 2) => 1, (:C, 1, 3) => 1, (:C, 2, 3) => 1)
        indexdict = Dict((:C, 2, 3, 1) => 3, (:Lj, 1, 3, 1) => 4, (:C, 1, 3, 1) => 5, (:R, 1, 2, 1) => 2, (:C, 3, 4, 1) => 6, (:P, 1, 2, 1) => 1, (:C, 1, 4, 1) => 8, (:Lj, 1, 4, 1) => 7)
        @test_throws(
            "Minimum junction too much smaller than average for WRSPICE.",
            JosephsonCircuits.calcCjIcmean(typevector, nodeindexarray,
                valuevector, namedict,mutualinductorvector, countdict, indexdict)
        )
    end

    @testset "calcCjIcmean" begin
        typevector = [:P, :R, :C, :Lj, :C, :C, :Lj, :C, :C, :Lj, :C, :C, :Lj, :C, :C, :Lj, :C, :C, :Lj, :C, :C, :Lj, :C, :C, :Lj, :C, :C, :Lj, :C, :C, :Lj, :C, :C, :Lj, :C, :C, :Lj, :C, :C, :Lj, :C, :C, :Lj, :C, :C, :Lj, :C, :C, :Lj, :C, :C, :Lj, :C, :C, :Lj, :C, :C, :Lj, :C, :C, :Lj, :C, :C, :Lj, :C, :C, :Lj, :C, :C, :Lj, :C, :C, :Lj, :C, :C, :Lj, :C, :C, :Lj, :C, :C, :Lj, :C, :C, :Lj, :C, :C, :Lj, :C, :C, :Lj, :C, :C, :Lj, :C, :C, :Lj, :C, :C, :Lj, :C, :C, :Lj, :C, :C, :Lj, :C, :C, :Lj, :C, :C, :Lj, :C, :C, :Lj, :C, :C, :Lj, :C, :C, :Lj, :C, :C, :Lj, :C, :C, :Lj, :C, :C, :Lj, :C, :C, :Lj, :C, :C, :Lj, :C, :C, :Lj, :C, :C, :Lj, :C, :C, :Lj, :C, :C, :Lj, :C, :C, :Lj, :C, :C, :Lj, :C, :C, :Lj, :C, :C, :Lj, :C, :C, :Lj, :C, :C, :R, :P]
        nodeindexarray = [2 2 2 2 2 3 3 3 4 4 4 5 5 5 6 6 6 7 7 7 8 8 8 9 9 9 10 10 10 11 11 11 12 12 12 13 13 13 14 14 14 15 15 15 16 16 16 17 17 17 18 18 18 19 19 19 20 20 20 21 21 21 22 22 22 23 23 23 24 24 24 25 25 25 26 26 26 27 27 27 28 28 28 29 29 29 30 30 30 31 31 31 32 32 32 33 33 33 34 34 34 35 35 35 36 36 36 37 37 37 38 38 38 39 39 39 40 40 40 41 41 41 42 42 42 43 43 43 44 44 44 45 45 45 46 46 46 47 47 47 48 48 48 49 49 49 50 50 50 51 51 51 52 52 52 53 53 53 54 54 54 55 55 55 56 56 56; 1 1 1 3 3 1 4 4 1 5 5 1 6 6 1 7 7 1 8 8 1 9 9 1 10 10 1 11 11 1 12 12 1 13 13 1 14 14 1 15 15 1 16 16 1 17 17 1 18 18 1 19 19 1 20 20 1 21 21 1 22 22 1 23 23 1 24 24 1 25 25 1 26 26 1 27 27 1 28 28 1 29 29 1 30 30 1 31 31 1 32 32 1 33 33 1 34 34 1 35 35 1 36 36 1 37 37 1 38 38 1 39 39 1 40 40 1 41 41 1 42 42 1 43 43 1 44 44 1 45 45 1 46 46 1 47 47 1 48 48 1 49 49 1 50 50 1 51 51 1 52 52 1 53 53 1 54 54 1 55 55 1 56 56 1 1 1]
        vvn = Real[1, 50.0, 2.25e-14, 9.679587529411765e-11, 5.5e-14, 4.5e-14, 9.679587529411765e-11, 5.5e-14, 4.5e-14, 9.679587529411765e-11, 5.5e-14, 4.5e-14, 9.679587529411765e-11, 5.5e-14, 4.5e-14, 9.679587529411765e-11, 5.5e-14, 4.5e-14, 9.679587529411765e-11, 5.5e-14, 4.5e-14, 9.679587529411765e-11, 5.5e-14, 4.5e-14, 9.679587529411765e-11, 5.5e-14, 4.5e-14, 9.679587529411765e-11, 5.5e-14, 4.5e-14, 9.679587529411765e-11, 5.5e-14, 4.5e-14, 9.679587529411765e-11, 5.5e-14, 4.5e-14, 9.679587529411765e-11, 5.5e-14, 4.5e-14, 9.679587529411765e-11, 5.5e-14, 4.5e-14, 9.679587529411765e-11, 5.5e-14, 4.5e-14, 9.679587529411765e-11, 5.5e-14, 4.5e-14, 9.679587529411765e-11, 5.5e-14, 4.5e-14, 9.679587529411765e-11, 5.5e-14, 4.5e-14, 9.679587529411765e-11, 5.5e-14, 4.5e-14, 9.679587529411765e-11, 5.5e-14, 4.5e-14, 9.679587529411765e-11, 5.5e-14, 4.5e-14, 9.679587529411765e-11, 5.5e-14, 4.5e-14, 9.679587529411765e-11, 5.5e-14, 4.5e-14, 9.679587529411765e-11, 5.5e-14, 4.5e-14, 9.679587529411765e-11, 5.5e-14, 4.5e-14, 9.679587529411765e-11, 5.5e-14, 4.5e-14, 9.679587529411765e-11, 5.5e-14, 4.5e-14, 9.679587529411765e-11, 5.5e-14, 4.5e-14, 9.679587529411765e-11, 5.5e-14, 4.5e-14, 9.679587529411765e-11, 5.5e-14, 4.5e-14, 9.679587529411765e-11, 5.5e-14, 4.5e-14, 9.679587529411765e-11, 5.5e-14, 4.5e-14, 9.679587529411765e-11, 5.5e-14, 4.5e-14, 9.679587529411765e-11, 5.5e-14, 4.5e-14, 9.679587529411765e-11, 5.5e-14, 4.5e-14, 9.679587529411765e-11, 5.5e-14, 4.5e-14, 9.679587529411765e-11, 5.5e-14, 4.5e-14, 9.679587529411765e-11, 5.5e-14, 4.5e-14, 9.679587529411765e-11, 5.5e-14, 4.5e-14, 9.679587529411765e-11, 5.5e-14, 4.5e-14, 9.679587529411765e-11, 5.5e-14, 4.5e-14, 9.679587529411765e-11, 5.5e-14, 4.5e-14, 9.679587529411765e-11, 5.5e-14, 4.5e-14, 9.679587529411765e-11, 5.5e-14, 4.5e-14, 9.679587529411765e-11, 5.5e-14, 4.5e-14, 9.679587529411765e-11, 5.5e-14, 4.5e-14, 9.679587529411765e-11, 5.5e-14, 4.5e-14, 9.679587529411765e-11, 5.5e-14, 4.5e-14, 9.679587529411765e-11, 5.5e-14, 4.5e-14, 9.679587529411765e-11, 5.5e-14, 4.5e-14, 9.679587529411764e-14, 5.5e-14, 4.5e-14, 9.679587529411765e-11, 5.5e-14, 4.5e-14, 9.679587529411765e-11, 5.5e-14, 4.5e-14, 9.679587529411765e-11, 5.5e-14, 4.5e-14, 9.679587529411765e-11, 5.5e-14, 2.25e-14, 50.0, 2]
        namedict = Dict("C21_0" => 63, "Lj35_36" => 106, "C52_53" => 158, "C4_0" => 12, "C9_10" => 29, "C29_0" => 87, "C51_52" => 155, "Lj22_23" => 67, "C9_0" => 27, "Lj47_48" => 142, "C18_0" => 54, "Lj1_2" => 4, "Lj33_34" => 100, "C33_0" => 99, "C13_0" => 39, "C35_0" => 105, "C31_0" => 93, "C40_0" => 120, "Lj32_33" => 97, "Lj46_47" => 139, "P55_0" => 167, "C32_0" => 96, "Lj24_25" => 73, "Lj10_11" => 31, "C17_0" => 51, "Lj9_10" => 28, "C6_7" => 20, "C37_0" => 111, "C26_0" => 78, "Lj3_4" => 10, "Lj8_9" => 25, "C19_0" => 57, "C3_0" => 9, "C35_36" => 107, "R1_0" => 2, "C17_18" => 53, "C38_0" => 114, "C10_0" => 30, "Lj53_54" => 160, "Lj23_24" => 70, "C41_42" => 125, "C21_22" => 65, "Lj40_41" => 121, "C12_0" => 36, "C44_45" => 134, "Lj5_6" => 16, "C54_55" => 164, "C19_20" => 59, "Lj34_35" => 103, "C42_0" => 126, "C46_47" => 140, "Lj4_5" => 13, "C5_6" => 17, "Lj26_27" => 79, "Lj39_40" => 118, "C34_0" => 102, "C8_9" => 26, "C52_0" => 156, "C11_0" => 33, "Lj28_29" => 85, "C2_0" => 6, "C45_46" => 137, "C24_0" => 72, "C7_8" => 23, "C23_24" => 71, "C24_25" => 74, "Lj30_31" => 91, "C16_17" => 50, "Lj50_51" => 151, "C36_37" => 110, "C39_40" => 119, "Lj44_45" => 133, "C47_48" => 143, "C49_0" => 147, "C51_0" => 153, "C53_54" => 161, "Lj54_55" => 163, "C41_0" => 123, "R55_0" => 166, "C43_44" => 131, "C3_4" => 11, "C23_0" => 69, "Lj43_44" => 130, "C38_39" => 116, "C15_16" => 47, "Lj17_18" => 52, "C28_29" => 86, "Lj15_16" => 46, "C1_0" => 3, "Lj25_26" => 76, "C22_0" => 66, "Lj16_17" => 49, "C30_0" => 90, "C14_15" => 44, "C26_27" => 80, "C12_13" => 38, "C22_23" => 68, "C2_3" => 8, "Lj12_13" => 37, "Lj11_12" => 34, "C47_0" => 141, "Lj49_50" => 148, "C5_0" => 15, "C6_0" => 18, "C46_0" => 138, "C50_0" => 150, "C48_49" => 146, "C25_26" => 77, "C34_35" => 104, "C54_0" => 162, "C49_50" => 149, "C1_2" => 5, "Lj27_28" => 82, "C28_0" => 84, "Lj29_30" => 88, "C36_0" => 108, "C27_0" => 81, "C8_0" => 24, "Lj36_37" => 109, "Lj2_3" => 7, "Lj37_38" => 112, "Lj38_39" => 115, "C42_43" => 128, "C50_51" => 152, "C7_0" => 21, "C16_0" => 48, "C18_19" => 56, "Lj51_52" => 154, "C53_0" => 159, "Lj20_21" => 61, "Lj6_7" => 19, "Lj14_15" => 43, "C39_0" => 117, "Lj18_19" => 55, "C44_0" => 132, "C43_0" => 129, "Lj19_20" => 58, "C20_0" => 60, "P1_0" => 1, "Lj7_8" => 22, "C11_12" => 35, "C31_32" => 95, "Lj48_49" => 145, "Lj52_53" => 157, "C48_0" => 144, "C33_34" => 101, "C32_33" => 98, "C27_28" => 83, "C45_0" => 135, "C15_0" => 45, "Lj31_32" => 94, "Lj41_42" => 124, "C29_30" => 89, "C20_21" => 62, "Lj42_43" => 127, "Lj45_46" => 136, "C4_5" => 14, "Lj13_14" => 40, "Lj21_22" => 64, "C40_41" => 122, "C14_0" => 42, "C30_31" => 92, "C37_38" => 113, "C55_0" => 165, "C10_11" => 32, "C13_14" => 41, "C25_0" => 75)
        mutualinductorvector = String[]
        countdict = Dict((:C, 1, 7) => 1, (:Lj, 27, 28) => 1, (:Lj, 25, 26) => 1, (:Lj, 42, 43) => 1, (:R, 1, 56) => 1, (:C, 1, 18) => 1, (:Lj, 23, 24) => 1, (:C, 6, 7) => 1, (:C, 1, 34) => 1, (:C, 1, 51) => 1, (:C, 55, 56) => 1, (:C, 1, 15) => 1, (:C, 52, 53) => 1, (:C, 1, 16) => 1, (:Lj, 18, 19) => 1, (:C, 1, 48) => 1, (:C, 1, 52) => 1, (:Lj, 41, 42) => 1, (:C, 14, 15) => 1, (:C, 20, 21) => 1, (:C, 5, 6) => 1, (:C, 1, 22) => 1, (:C, 22, 23) => 1, (:C, 1, 38) => 1, (:C, 41, 42) => 1, (:C, 1, 5) => 1, (:Lj, 30, 31) => 1, (:Lj, 50, 51) => 1, (:C, 1, 6) => 1, (:C, 1, 26) => 1, (:C, 38, 39) => 1, (:Lj, 46, 47) => 1, (:Lj, 52, 53) => 1, (:C, 30, 31) => 1, (:C, 1, 46) => 1, (:C, 1, 47) => 1, (:C, 1, 37) => 1, (:C, 1, 41) => 1, (:C, 4, 5) => 1, (:Lj, 6, 7) => 1, (:Lj, 15, 16) => 1, (:Lj, 55, 56) => 1, (:C, 39, 40) => 1, (:C, 7, 8) => 1, (:C, 1, 4) => 1, (:Lj, 13, 14) => 1, (:C, 1, 36) => 1, (:C, 1, 40) => 1, (:C, 1, 17) => 1, (:C, 50, 51) => 1, (:Lj, 5, 6) => 1, (:Lj, 19, 20) => 1, (:Lj, 38, 39) => 1, (:C, 26, 27) => 1, (:C, 1, 49) => 1, (:C, 23, 24) => 1, (:Lj, 29, 30) => 1, (:C, 2, 3) => 1, (:Lj, 35, 36) => 1, (:Lj, 49, 50) => 1, (:C, 1, 56) => 1, (:C, 1, 29) => 1, (:C, 29, 30) => 1, (:Lj, 16, 17) => 1, (:C, 12, 13) => 1, (:C, 49, 50) => 1, (:R, 1, 2) => 1, (:C, 42, 43) => 1, (:C, 21, 22) => 1, (:C, 10, 11) => 1, (:C, 1, 24) => 1, (:Lj, 2, 3) => 1, (:Lj, 11, 12) => 1, (:C, 1, 23) => 1, (:Lj, 54, 55) => 1, (:C, 31, 32) => 1, (:Lj, 26, 27) => 1, (:Lj, 28, 29) => 1, (:C, 1, 28) => 1, (:C, 1, 19) => 1, (:Lj, 48, 49) => 1, (:C, 34, 35) => 1, (:C, 1, 25) => 1, (:C, 1, 33) => 1, (:C, 47, 48) => 1, (:C, 32, 33) => 1, (:C, 1, 12) => 1, (:Lj, 7, 8) => 1, (:C, 24, 25) => 1, (:C, 17, 18) => 1, (:C, 1, 20) => 1, (:C, 1, 31) => 1, (:C, 40, 41) => 1, (:C, 1, 32) => 1, (:Lj, 31, 32) => 1, (:Lj, 39, 40) => 1, (:C, 1, 54) => 1, (:C, 1, 3) => 1, (:Lj, 9, 10) => 1, (:C, 46, 47) => 1, (:Lj, 43, 44) => 1, (:Lj, 44, 45) => 1, (:C, 1, 9) => 1, (:C, 13, 14) => 1, (:C, 27, 28) => 1, (:C, 1, 11) => 1, (:Lj, 17, 18) => 1, (:C, 9, 10) => 1, (:C, 1, 35) => 1, (:C, 54, 55) => 1, (:Lj, 32, 33) => 1, (:C, 35, 36) => 1, (:Lj, 14, 15) => 1, (:C, 44, 45) => 1, (:C, 16, 17) => 1, (:C, 15, 16) => 1, (:Lj, 21, 22) => 1, (:Lj, 36, 37) => 1, (:C, 43, 44) => 1, (:C, 33, 34) => 1, (:C, 48, 49) => 1, (:Lj, 37, 38) => 1, (:Lj, 8, 9) => 1, (:C, 1, 42) => 1, (:C, 1, 27) => 1, (:Lj, 51, 52) => 1, (:C, 1, 53) => 1, (:C, 28, 29) => 1, (:C, 37, 38) => 1, (:C, 1, 43) => 1, (:C, 1, 14) => 1, (:Lj, 45, 46) => 1, (:C, 1, 21) => 1, (:Lj, 40, 41) => 1, (:C, 1, 45) => 1, (:C, 1, 50) => 1, (:C, 51, 52) => 1, (:Lj, 33, 34) => 1, (:Lj, 20, 21) => 1, (:Lj, 47, 48) => 1, (:Lj, 53, 54) => 1, (:C, 1, 13) => 1, (:C, 1, 30) => 1, (:C, 53, 54) => 1, (:C, 1, 8) => 1, (:Lj, 34, 35) => 1, (:C, 25, 26) => 1, (:Lj, 3, 4) => 1, (:C, 1, 10) => 1, (:Lj, 4, 5) => 1, (:C, 8, 9) => 1, (:C, 11, 12) => 1, (:C, 36, 37) => 1, (:C, 1, 44) => 1, (:C, 1, 2) => 1, (:C, 3, 4) => 1, (:C, 45, 46) => 1, (:C, 1, 55) => 1, (:Lj, 22, 23) => 1, (:P, 1, 56) => 1, (:Lj, 24, 25) => 1, (:C, 18, 19) => 1, (:C, 19, 20) => 1, (:P, 1, 2) => 1, (:C, 1, 39) => 1, (:Lj, 10, 11) => 1, (:Lj, 12, 13) => 1)
        indexdict = Dict((:C, 12, 13, 1) => 35, (:Lj, 35, 36, 1) => 103, (:Lj, 46, 47, 1) => 136, (:C, 55, 56, 1) => 164, (:Lj, 39, 40, 1) => 115, (:C, 1, 32, 1) => 93, (:Lj, 42, 43, 1) => 124, (:C, 27, 28, 1) => 80, (:C, 35, 36, 1) => 104, (:Lj, 25, 26, 1) => 73, (:C, 49, 50, 1) => 146, (:Lj, 21, 22, 1) => 61, (:C, 1, 22, 1) => 63, (:C, 4, 5, 1) => 11, (:Lj, 9, 10, 1) => 25, (:Lj, 37, 38, 1) => 109, (:Lj, 3, 4, 1) => 7, (:Lj, 26, 27, 1) => 76, (:C, 53, 54, 1) => 158, (:C, 1, 25, 1) => 72, (:C, 1, 5, 1) => 12, (:C, 7, 8, 1) => 20, (:C, 10, 11, 1) => 29, (:C, 1, 18, 1) => 51, (:Lj, 8, 9, 1) => 22, (:C, 1, 42, 1) => 123, (:C, 18, 19, 1) => 53, (:C, 17, 18, 1) => 50, (:C, 1, 35, 1) => 102, (:C, 50, 51, 1) => 149, (:Lj, 7, 8, 1) => 19, (:C, 38, 39, 1) => 113, (:Lj, 47, 48, 1) => 139, (:Lj, 32, 33, 1) => 94, (:C, 46, 47, 1) => 137, (:C, 1, 50, 1) => 147, (:C, 1, 16, 1) => 45, (:C, 40, 41, 1) => 119, (:Lj, 55, 56, 1) => 163, (:C, 1, 2, 1) => 3, (:C, 48, 49, 1) => 143, (:Lj, 15, 16, 1) => 43, (:Lj, 20, 21, 1) => 58, (:C, 26, 27, 1) => 77, (:C, 23, 24, 1) => 68, (:C, 22, 23, 1) => 65, (:Lj, 24, 25, 1) => 70, (:C, 1, 38, 1) => 111, (:Lj, 48, 49, 1) => 142, (:Lj, 4, 5, 1) => 10, (:C, 1, 31, 1) => 90, (:C, 1, 45, 1) => 132, (:Lj, 49, 50, 1) => 145, (:C, 1, 34, 1) => 99, (:C, 1, 51, 1) => 150, (:Lj, 54, 55, 1) => 160, (:Lj, 11, 12, 1) => 31, (:C, 33, 34, 1) => 98, (:C, 6, 7, 1) => 17, (:C, 1, 3, 1) => 6, (:C, 3, 4, 1) => 8, (:Lj, 14, 15, 1) => 40, (:Lj, 6, 7, 1) => 16, (:C, 25, 26, 1) => 74, (:C, 1, 53, 1) => 156, (:C, 44, 45, 1) => 131, (:C, 5, 6, 1) => 14, (:Lj, 18, 19, 1) => 52, (:C, 36, 37, 1) => 107, (:C, 1, 13, 1) => 36, (:Lj, 17, 18, 1) => 49, (:C, 16, 17, 1) => 47, (:Lj, 19, 20, 1) => 55, (:Lj, 27, 28, 1) => 79, (:C, 9, 10, 1) => 26, (:C, 37, 38, 1) => 110, (:C, 1, 12, 1) => 33, (:C, 32, 33, 1) => 95, (:R, 1, 2, 1) => 2, (:C, 1, 9, 1) => 24, (:C, 1, 33, 1) => 96, (:C, 1, 11, 1) => 30, (:C, 1, 27, 1) => 78, (:C, 21, 22, 1) => 62, (:Lj, 38, 39, 1) => 112, (:Lj, 44, 45, 1) => 130, (:C, 1, 46, 1) => 135, (:C, 1, 49, 1) => 144, (:Lj, 30, 31, 1) => 88, (:Lj, 28, 29, 1) => 82, (:C, 1, 44, 1) => 129, (:C, 1, 56, 1) => 165, (:C, 39, 40, 1) => 116, (:C, 2, 3, 1) => 5, (:C, 1, 8, 1) => 21, (:C, 47, 48, 1) => 140, (:C, 1, 28, 1) => 81, (:P, 1, 56, 1) => 167, (:C, 52, 53, 1) => 155, (:R, 1, 56, 1) => 166, (:Lj, 2, 3, 1) => 4, (:C, 1, 4, 1) => 9, (:C, 1, 6, 1) => 15, (:C, 1, 10, 1) => 27, (:C, 1, 19, 1) => 54, (:Lj, 16, 17, 1) => 46, (:Lj, 22, 23, 1) => 64, (:C, 31, 32, 1) => 92, (:C, 14, 15, 1) => 41, (:Lj, 40, 41, 1) => 118, (:C, 1, 29, 1) => 84, (:Lj, 50, 51, 1) => 148, (:C, 1, 23, 1) => 66, (:C, 1, 24, 1) => 69, (:Lj, 33, 34, 1) => 97, (:C, 1, 40, 1) => 117, (:C, 20, 21, 1) => 59, (:C, 1, 47, 1) => 138, (:Lj, 52, 53, 1) => 154, (:Lj, 36, 37, 1) => 106, (:C, 1, 21, 1) => 60, (:C, 8, 9, 1) => 23, (:C, 1, 7, 1) => 18, (:C, 19, 20, 1) => 56, (:Lj, 31, 32, 1) => 91, (:Lj, 45, 46, 1) => 133, (:C, 45, 46, 1) => 134, (:C, 43, 44, 1) => 128, (:P, 1, 2, 1) => 1, (:Lj, 10, 11, 1) => 28, (:C, 1, 14, 1) => 39, (:C, 1, 30, 1) => 87, (:Lj, 34, 35, 1) => 100, (:C, 1, 36, 1) => 105, (:C, 1, 37, 1) => 108, (:C, 15, 16, 1) => 44, (:C, 24, 25, 1) => 71, (:C, 13, 14, 1) => 38, (:C, 1, 39, 1) => 114, (:C, 42, 43, 1) => 125, (:C, 1, 26, 1) => 75, (:Lj, 41, 42, 1) => 121, (:C, 1, 48, 1) => 141, (:C, 51, 52, 1) => 152, (:C, 1, 54, 1) => 159, (:C, 54, 55, 1) => 161, (:C, 1, 52, 1) => 153, (:Lj, 5, 6, 1) => 13, (:Lj, 13, 14, 1) => 37, (:Lj, 29, 30, 1) => 85, (:Lj, 12, 13, 1) => 34, (:C, 1, 15, 1) => 42, (:C, 28, 29, 1) => 83, (:C, 1, 20, 1) => 57, (:C, 1, 17, 1) => 48, (:C, 1, 43, 1) => 126, (:Lj, 53, 54, 1) => 157, (:C, 1, 55, 1) => 162, (:Lj, 43, 44, 1) => 127, (:C, 30, 31, 1) => 89, (:Lj, 51, 52, 1) => 151, (:C, 41, 42, 1) => 122, (:C, 34, 35, 1) => 101, (:C, 29, 30, 1) => 86, (:C, 1, 41, 1) => 120, (:Lj, 23, 24, 1) => 67, (:C, 11, 12, 1) => 32)

        @test_throws(
            "Maximum junction too much larger than average for WRSPICE.",
            JosephsonCircuits.calcCjIcmean(typevector, nodeindexarray,
                vvn, namedict,mutualinductorvector, countdict, indexdict)
        )
    end

end