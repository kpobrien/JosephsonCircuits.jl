using JosephsonCircuits
using SparseArrays
using Test

@testset verbose=true "JosephsonCircuits" begin

    @testset verbose=true "warmup" begin
        @variables Ipump Rleft Cc Lj Cj w L1
        out1 = JosephsonCircuits.HB(JosephsonCircuits.NonlinearHB((2.9845193040956104e10,), JosephsonCircuits.Frequencies{1}((4,), (5,), (8,), CartesianIndex{1}[CartesianIndex(2,), CartesianIndex(4,)], [(1,), (3,)]), ComplexF64[-0.013189713633103771 - 0.00865114720807541im, 2.6385352759232326e-5 - 5.730931597819341e-6im, 0.12157328255345531 + 0.07973637714134815im, 1.3583890694738675e-5 - 6.466918389420335e-5im], sparse([1, 2, 3, 4], [1, 2, 3, 4], [1, 1, 1, 1], 4, 4), sparsevec([2], ComplexF64[1.0e-9 + 0.0im], 2), sparsevec(Int64[], Nothing[], 2), sparsevec([3, 4], ComplexF64[1.0e-9 + 0.0im, 1.0e-9 + 0.0im], 4), 2, 2, ["1", "2"], [1], [(1,), (3,)], ComplexF64[-0.3984172062272202 - 0.9171852686842317im -0.0 - 0.0im; 0.0006902518137315259 + 0.0031779366560466767im 0.0 - 0.0im]), JosephsonCircuits.LinearHB([0.8793432544356942 - 0.47684762507575995im 0.020027766663628266 - 0.015070010763669602im; 0.02468311994730562 + 0.00435433710417341im 0.9999645263305277 + 0.026441685697269325im;;; 0.9999651737600814 - 0.026420913217696473im 0.024687073668703135 - 0.004354512827247238im; 0.02003059333090143 + 0.015072780735441067im 0.8793336998975008 + 0.47686545021930254im], Array{ComplexF64, 3}(undef, 0, 0, 0), [0.999372571659925 0.0006274283400750277; 0.0006274283400750285 0.999372571659925;;; 0.9993723754270823 0.0006276245729177577; 0.0006276245729177582 0.9993723754270823], [0.9993725716599243 1.0; 1.0 0.999372571659925;;; 0.9993723754270819 1.0; 1.0 0.9993723754270828], [1.0000000000000007 1.0000000000000004; -0.9999999999999999 -0.9999999999999994], ComplexF64[], ComplexF64[], ComplexF64[], 2, 3, 2, ["1", "2"], [1], 1, 2.827433388230814e10:3.141592653589793e9:3.141592653589793e10, [(0,), (-2,)]))
        out2 = JosephsonCircuits.warmup()
        @test JosephsonCircuits.compare(out1,out2)
    end

    @testset verbose=true "warmupsymsold" begin
        @variables Rleft Cc Lj Cj w L1
        out1 = JosephsonCircuits.HB(JosephsonCircuits.NonlinearHB((2.9845193040956104e10,), JosephsonCircuits.Frequencies{1}((4,), (5,), (8,), CartesianIndex{1}[CartesianIndex(2,), CartesianIndex(4,)], [(1,), (3,)]), ComplexF64[-0.01317268153461008 - 0.008620192636063851im, 3.3877483616232297e-6 + 8.323112211967887e-6im, 0.12179774887406367 + 0.07965319542032037im, 2.1979490981718878e-5 + 7.557330878003312e-7im], sparse([1, 2, 3, 4], [1, 2, 3, 4], [1, 1, 1, 1], 4, 4), sparsevec([2], [1.0e-9], 2), sparsevec(Int64[], Nothing[], 2), sparsevec([3, 4], [1.0e-9, 1.0e-9], 4), 2, 2, ["1", "2"], [1], [(1,), (3,)], ComplexF64[-0.4005697228198604 - 0.9160008919594795im -0.0 - 0.0im; -0.0010024623749457897 + 0.0004080312966859029im 0.0 + 0.0im]), JosephsonCircuits.LinearHB([0.9999636971626665 + 0.026577597828732157im 0.024776280908444095 + 0.004460970065009537im; 0.02006351652496576 - 0.01520590853015591im 0.8793137132801259 - 0.4769079135295333im;;; 0.8793041537769276 + 0.4769257467844491im 0.020066345817057213 + 0.015208700168465511im; 0.0247802510228797 - 0.004461160668153754im 0.9999643476849712 - 0.02655684759593593im], Array{ComplexF64, 3}(undef, 0, 0, 0), [0.9993670379479913 0.0006329620520086488; 0.0006329620520086476 0.9993670379479913;;; 0.9993668400044972 0.0006331599955028162; 0.0006331599955028159 0.9993668400044972], [0.9993670379479914 1.0; 1.0 0.9993670379479901;;; 0.9993668400044982 1.0; 1.0 0.9993668400044973], [-0.9999999999999998 -0.9999999999999991; 1.000000000000001 0.9999999999999998], ComplexF64[], ComplexF64[], ComplexF64[], 2, 3, 2, ["1", "2"], [1], 2, 2.827433388230814e10:3.141592653589793e9:3.141592653589793e10, [(-2,), (0,)]))
        out2 = JosephsonCircuits.warmupsymsold()
        @test JosephsonCircuits.compare(out1,out2)
    end

    @testset verbose=true "warmupsyms" begin
        @variables Rleft Cc Lj Cj w L1
        out1 = JosephsonCircuits.HB(JosephsonCircuits.NonlinearHB((2.9845193040956104e10,), JosephsonCircuits.Frequencies{1}((4,), (5,), (8,), CartesianIndex{1}[CartesianIndex(2,), CartesianIndex(4,)], [(1,), (3,)]), ComplexF64[-0.013189713633104328 - 0.008651147208076224im, 2.638535275923149e-5 - 5.730931597816161e-6im, 0.1215732825534493 + 0.07973637714135108im, 1.358389069474495e-5 - 6.46691838941983e-5im], sparse([1, 2, 3, 4], [1, 2, 3, 4], [1, 1, 1, 1], 4, 4), sparsevec([2], [1.0e-9], 2), sparsevec(Int64[], Nothing[], 2), sparsevec([3, 4], [1.0e-9, 1.0e-9], 4), 2, 2, ["1", "2"], [1], [(1,), (3,)], ComplexF64[-0.3984172062271636 - 0.9171852686842704im -0.0 - 0.0im; 0.0006902518137311428 + 0.0031779366560465757im 0.0 - 0.0im]), JosephsonCircuits.LinearHB([0.8793432544356953 - 0.4768476250757582im 0.02002776666362855 - 0.015070010763667296im; 0.024683119947304828 + 0.004354337104171201im 0.9999645263305277 + 0.026441685697266244im;;; 0.9999651737600814 - 0.0264209132176933im 0.02468707366870235 - 0.004354512827245029im; 0.02003059333090171 + 0.015072780735438763im 0.8793336998975017 + 0.4768654502193009im], Array{ComplexF64, 3}(undef, 0, 0, 0), [0.999372571659925 0.0006274283400749695; 0.0006274283400749705 0.999372571659925;;; 0.9993723754270823 0.0006276245729177001; 0.0006276245729177002 0.9993723754270823], [0.9993725716599241 1.0; 1.0 0.9993725716599252;;; 0.9993723754270819 1.0; 1.0 0.9993723754270828], [1.0000000000000009 1.0000000000000002; -0.9999999999999998 -0.9999999999999993], ComplexF64[], ComplexF64[], ComplexF64[], 2, 3, 2, ["1", "2"], [1], 1, 2.827433388230814e10:3.141592653589793e9:3.141592653589793e10, [(0,), (-2,)]))
        out2 = JosephsonCircuits.warmupsyms()
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
        out1 = JosephsonCircuits.LinearHB([0.8952708641229391 - 0.44552225517090227im;;; 0.8415115570832487 - 0.5402391130743189im;;; 0.6457820691998714 - 0.7635217869189669im;;; -0.9968560060568034 + 0.07923448231975308im;;; 0.9316787544566122 + 0.36328322077158454im;;; 0.9988570509555925 + 0.04779740323801577im], Array{ComplexF64, 3}(undef, 0, 0, 0), [1.0;;; 1.0;;; 1.0;;; 1.0;;; 1.0;;; 1.0], [0.9999999999999993;;; 0.9999999999999996;;; 1.0;;; 0.9999999999999991;;; 1.0;;; 0.9999999999999993], [1.0000000000000007 1.0000000000000004 0.9999999999999997 1.0000000000000009 1.0 1.0000000000000007], ComplexF64[], ComplexF64[], ComplexF64[], 1, 3, 2, ["1", "2"], [1], 1, 2.827433388230814e10:6.283185307179586e8:3.141592653589793e10, [(0,)])
        out2 = JosephsonCircuits.warmuphblinsolve()
        @test JosephsonCircuits.compare(out1,out2)
    end

    @testset verbose=true "warmupvvn" begin
        out1 = Real[1, 1.0e-8, 50.0, 1.0e-13, 1.0e-9, 1.0e-12]
        out2 = JosephsonCircuits.warmupvvn()
        @test JosephsonCircuits.compare(out1,out2)
    end

end