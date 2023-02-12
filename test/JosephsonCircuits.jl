using JosephsonCircuits
using SparseArrays
using Test

@testset verbose=true "JosephsonCircuits" begin

    @testset verbose=true "warmup" begin
        @variables Ipump Rleft Cc Lj Cj w L1
        out1 = JosephsonCircuits.HB(JosephsonCircuits.NonlinearHB((2.9845193040956104e10,), JosephsonCircuits.Frequencies{1}((4,), (5,), (8,), CartesianIndex{1}[CartesianIndex(2,), CartesianIndex(4,)], [(1,), (3,)]), ComplexF64[-0.013172681534609886 - 0.008620192636062398im, 3.387748361622929e-6 + 8.323112211960799e-6im, 0.12179774887407359 + 0.07965319542032055im, 2.1979490981702754e-5 + 7.557330877939172e-7im], sparse([1, 2, 3, 4], [1, 2, 3, 4], [1, 1, 1, 1], 4, 4), sparsevec([2], ComplexF64[1.0e-9 + 0.0im], 2), sparsevec(Int64[], Nothing[], 2), sparsevec([3, 4], ComplexF64[1.0e-9 + 0.0im, 1.0e-9 + 0.0im], 4), 2, 2, ComplexF64[0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im]), JosephsonCircuits.LinearHB([0.9999636971626665 + 0.02657759782873958im 0.02477628090844665 + 0.004460970065011846im; 0.020063516524967025 - 0.015205908530159075im 0.8793137132801236 - 0.4769079135295377im;;; 0.8793041537769256 + 0.4769257467844534im 0.020066345817058493 + 0.015208700168468686im; 0.024780251022882258 - 0.004461160668156066im 0.9999643476849712 - 0.02655684759594325im], Array{ComplexF64, 3}(undef, 0, 0, 0), [0.9993670379479911 0.0006329620520087953; 0.0006329620520087943 0.9993670379479911;;; 0.9993668400044972 0.0006331599955029636; 0.000633159995502963 0.9993668400044972], [0.999367037947991 1.0; 1.0 0.9993670379479901;;; 0.9993668400044973 1.0; 1.0 0.9993668400044973], [-1.0 -0.9999999999999997; 1.000000000000001 0.9999999999999999], ComplexF64[], ComplexF64[], ComplexF64[], 2, 3, 2, 2, 2.827433388230814e10:3.141592653589793e9:3.141592653589793e10, [(-2,), (0,)]))
        out2 = JosephsonCircuits.warmup()
        @test JosephsonCircuits.compare(out1,out2)
    end

    @testset verbose=true "warmupsyms" begin
        @variables Rleft Cc Lj Cj w L1
        out1 = JosephsonCircuits.HB(JosephsonCircuits.NonlinearHB((2.9845193040956104e10,), JosephsonCircuits.Frequencies{1}((4,), (5,), (8,), CartesianIndex{1}[CartesianIndex(2,), CartesianIndex(4,)], [(1,), (3,)]), ComplexF64[-0.01317268153461008 - 0.008620192636063851im, 3.3877483616232297e-6 + 8.323112211967887e-6im, 0.12179774887406367 + 0.07965319542032037im, 2.1979490981718878e-5 + 7.557330878003312e-7im], sparse([1, 2, 3, 4], [1, 2, 3, 4], [1, 1, 1, 1], 4, 4), sparsevec([2], [1.0e-9], 2), sparsevec(Int64[], Nothing[], 2), sparsevec([3, 4], [1.0e-9, 1.0e-9], 4), 2, 2, ComplexF64[0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im]), JosephsonCircuits.LinearHB([0.9999636971626665 + 0.026577597828732157im 0.024776280908444095 + 0.004460970065009537im; 0.02006351652496576 - 0.01520590853015591im 0.8793137132801259 - 0.4769079135295333im;;; 0.8793041537769276 + 0.4769257467844491im 0.020066345817057213 + 0.015208700168465511im; 0.0247802510228797 - 0.004461160668153754im 0.9999643476849712 - 0.02655684759593593im], Array{ComplexF64, 3}(undef, 0, 0, 0), [0.9993670379479913 0.0006329620520086488; 0.0006329620520086476 0.9993670379479913;;; 0.9993668400044972 0.0006331599955028162; 0.0006331599955028159 0.9993668400044972], [0.9993670379479914 1.0; 1.0 0.9993670379479901;;; 0.9993668400044982 1.0; 1.0 0.9993668400044973], [-0.9999999999999998 -0.9999999999999991; 1.000000000000001 0.9999999999999998], ComplexF64[], ComplexF64[], ComplexF64[], 2, 3, 2, 2, 2.827433388230814e10:3.141592653589793e9:3.141592653589793e10, [(-2,), (0,)]))
        out2 = JosephsonCircuits.warmupsyms()
        @test JosephsonCircuits.compare(out1,out2)
    end

    @testset verbose=true "warmupsyms2" begin
        @variables Rleft Cc Lj Cj w L1
        out1 = JosephsonCircuits.HB(JosephsonCircuits.NonlinearHB((2.9845193040956104e10,), JosephsonCircuits.Frequencies{1}((4,), (5,), (8,), CartesianIndex{1}[CartesianIndex(2,), CartesianIndex(4,)], [(1,), (3,)]), ComplexF64[-0.013189713633104328 - 0.008651147208076224im, 2.638535275923149e-5 - 5.730931597816161e-6im, 0.1215732825534493 + 0.07973637714135108im, 1.358389069474495e-5 - 6.46691838941983e-5im], sparse([1, 2, 3, 4], [1, 2, 3, 4], [1, 1, 1, 1], 4, 4), sparsevec([2], [1.0e-9], 2), sparsevec(Int64[], Nothing[], 2), sparsevec([3, 4], [1.0e-9, 1.0e-9], 4), 2, 2, ComplexF64[0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im]), JosephsonCircuits.LinearHB([0.8793569504109443 - 0.47681967169273876im 0.001512328226489033 - 0.00032760895856871557im 0.020025158424520146 - 0.015067615722927218im; -0.0007934007343287591 + 0.0013286644026411442im 0.7304798277512885 - 0.6829325423534629im -0.0001705218151428818 - 4.3738391467828484e-5im; 0.024679629873200038 + 0.0043538778518115665im -6.37823775245524e-5 + 0.00016295097095725243im 0.9999644639585391 + 0.026441290141526786im;;; 0.9999647914129438 - 0.02638429641209193im -0.0016915420142723293 - 1.2595516192949867e-5im 0.024690396695749902 - 0.004354416922277141im; 0.0011629933572913824 - 0.0012285222260809im 0.7132140455028935 - 0.7009442824244202im -0.00010956779790372759 - 0.00010761689252906938im; 0.020032903114943077 + 0.01507513767399213im -0.00010855760658972386 + 0.00010702508091378949im 0.8793339176946257 + 0.47686524450209766im], Array{ComplexF64, 3}(undef, 0, 0, 0), [0.9993703564049484 2.391460428690098e-6 0.0006272521346227394; 2.394833671643844e-6 0.9999975741755938 3.09907344069648e-8; 0.000627252464855619 3.0582794285956745e-8 0.9993727169523501;;; 0.9993693547090311 2.857880239584042e-6 0.0006277874107293358; 2.861820274078016e-6 0.9999971145932292 2.358649678243965e-8; 0.0006277877280917366 2.320994234326287e-8 0.9993721890619659], [0.9993751363371318 1.0 1.0; 1.0 1.0 1.0; 1.0 1.0 0.9993727169523507;;; 0.9993750668975425 1.0 1.0; 1.0 1.0 1.0; 1.0 1.0 0.9993721890619665], [1.0000000000000004 1.0; 1.0000000000000002 1.0; -0.9999999999999996 -0.9999999999999996], ComplexF64[], ComplexF64[], ComplexF64[], 3, 3, 2, 2, 2.827433388230814e10:3.141592653589793e9:3.141592653589793e10, [(0,), (2,), (-2,)]))
        out2 = JosephsonCircuits.warmupsyms2()
        @test JosephsonCircuits.compare(out1,out2)
    end

    @testset verbose=true "warmupparse" begin
        @variables Ipump Rleft Cc Lj Cj w L1
        out1 = JosephsonCircuits.ParsedCircuit([1, 2, 1, 2, 1, 2, 1, 3, 3, 2, 3, 2], ["1", "0", "2"], String[], ["P1", "I1", "R1", "C1", "Lj1", "C2"], [:P, :I, :R, :C, :Lj, :C], Num[1, Ipump, Rleft, Cc, Lj, Cj], Dict("I1" => 2, "C1" => 4, "C2" => 6, "R1" => 3, "P1" => 1, "Lj1" => 5), 3)
        out2 = JosephsonCircuits.warmupparse()
        @test JosephsonCircuits.compare(out1,out2)
    end

    @testset verbose=true "warmupparsesort" begin
        @variables Ipump Rleft Cc Lj Cj w L1
        out1 = JosephsonCircuits.ParsedSortedCircuit([2 2 2 2 3 3; 1 1 1 3 1 1], ["0", "1", "2"], String[], ["P1", "I1", "R1", "C1", "Lj1", "C2"], [:P, :I, :R, :C, :Lj, :C], Num[1, Ipump, Rleft, Cc, Lj, Cj], Dict("I1" => 2, "C1" => 4, "C2" => 6, "R1" => 3, "P1" => 1, "Lj1" => 5), 3)
        out2 = JosephsonCircuits.warmupparsesort()
        @test JosephsonCircuits.compare(out1,out2)
    end

    @testset verbose=true "warmupnumericmatrices" begin
        out1 = JosephsonCircuits.CircuitMatrices(sparse([1, 2, 1, 2], [1, 1, 2, 2], [1.0e-13, -1.0e-13, -1.0e-13, 1.1e-12], 2, 2), sparse([1], [1], [0.02], 2, 2), sparsevec(Int64[], Nothing[], 2), sparsevec(Int64[], Nothing[], 2), sparsevec([2], [1.0e-9], 2), sparsevec([2], [1.0e-9], 2), sparse(Int64[], Int64[], Nothing[], 2, 2), sparse(Int64[], Int64[], Nothing[], 2, 2), sparse([1, 2], [1, 2], [1, 1], 2, 2), [1], [1], [3], Int64[], 1.0e-9, Real[1, 1.0e-8, 50.0, 1.0e-13, 1.0e-9, 1.0e-12])
        out2 = JosephsonCircuits.warmupnumericmatrices()
        @test JosephsonCircuits.compare(out1,out2)

    end

    @testset verbose=true "warmuphblinsolve" begin
        out1 = JosephsonCircuits.LinearHB([0.8952708641229391 - 0.44552225517090227im;;; 0.8415115570832487 - 0.5402391130743189im;;; 0.6457820691998714 - 0.7635217869189669im;;; -0.9968560060568034 + 0.07923448231975308im;;; 0.9316787544566122 + 0.36328322077158454im;;; 0.9988570509555925 + 0.04779740323801577im], Array{ComplexF64, 3}(undef, 0, 0, 0), [1.0;;; 1.0;;; 1.0;;; 1.0;;; 1.0;;; 1.0], [0.9999999999999993;;; 0.9999999999999996;;; 1.0;;; 0.9999999999999991;;; 1.0;;; 0.9999999999999993], [1.0000000000000007 1.0000000000000004 0.9999999999999997 1.0000000000000009 1.0 1.0000000000000007], ComplexF64[], ComplexF64[], ComplexF64[], 1, 3, 2, 1, 2.827433388230814e10:6.283185307179586e8:3.141592653589793e10, [(0,)])
        out2 = JosephsonCircuits.warmuphblinsolve()
        @test JosephsonCircuits.compare(out1,out2)
    end

    @testset verbose=true "warmupvvn" begin
        out1 = Real[1, 1.0e-8, 50.0, 1.0e-13, 1.0e-9, 1.0e-12]
        out2 = JosephsonCircuits.warmupvvn()
        @test JosephsonCircuits.compare(out1,out2)
    end

end