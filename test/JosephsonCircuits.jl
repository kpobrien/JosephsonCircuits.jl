using JosephsonCircuits
using SparseArrays
using Test

@testset verbose=true "JosephsonCircuits" begin

    @testset verbose=true "warmup" begin
        @variables Ipump Rleft Cc Lj Cj w L1
        out1 = JosephsonCircuits.HB(JosephsonCircuits.NonlinearHB(ComplexF64[-0.013172681534609886 - 0.008620192636062398im, 3.387748361622929e-6 + 8.323112211960799e-6im, 0.12179774887407359 + 0.07965319542032055im, 2.1979490981702754e-5 + 7.557330877939172e-7im], sparse([1, 2, 3, 4], [1, 2, 3, 4], [1, 1, 1, 1], 4, 4), sparsevec([2], ComplexF64[1.0e-9 + 0.0im], 2), sparsevec(Int64[], Nothing[], 2), sparsevec([3, 4], ComplexF64[1.0e-9 + 0.0im, 1.0e-9 + 0.0im], 4), 2, 2, ComplexF64[0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im]), ComplexF64[0.9789325751107415 + 0.0im; 0.1205118781611217 + 0.07881274759778024im; -0.0042304188751325356 - 0.009645207525001455im; 0.00019110015334428132 + 0.0im;;], JosephsonCircuits.LinearHB([0.9999636971626665 + 0.02657759782873958im 0.02477628090844665 + 0.004460970065011846im; 0.020063516524967025 - 0.015205908530159075im 0.8793137132801236 - 0.4769079135295377im;;; 0.8793041537769256 + 0.4769257467844534im 0.020066345817058493 + 0.015208700168468686im; 0.024780251022882258 - 0.004461160668156066im 0.9999643476849712 - 0.02655684759594325im], Array{ComplexF64, 3}(undef, 0, 0, 0), [0.9993670379479911 0.0006329620520087953; 0.0006329620520087943 0.9993670379479911;;; 0.9993668400044972 0.0006331599955029636; 0.000633159995502963 0.9993668400044972], [0.999367037947991 1.0; 1.0 0.9993670379479901;;; 0.9993668400044973 1.0; 1.0 0.9993668400044973], [-1.0 -0.9999999999999997; 1.000000000000001 0.9999999999999999], ComplexF64[], ComplexF64[], 2, 3, 2, 2, StepRangeLen(Base.TwicePrecision{Float64}(2.827433388230814e10, -1.7570437194081023e-6),Base.TwicePrecision{Float64}(3.141592653589793e9, -8.926326700020581e-8),2)))
        out2 = JosephsonCircuits.warmup()
        @test JosephsonCircuits.compare(out1,out2)
    end

    # @testset verbose=true "warmupsyms" begin
    #     # @variables Ipump Rleft Cc Lj Cj w L1
    #     out1 = JosephsonCircuits.HB(JosephsonCircuits.NonlinearHB(ComplexF64[-0.01317268153461008 - 0.008620192636063851im, 3.3877483616232297e-6 + 8.323112211967887e-6im, 0.12179774887406367 + 0.07965319542032037im, 2.1979490981718878e-5 + 7.557330878003312e-7im], sparse([1, 2, 3, 4], [1, 2, 3, 4], [1, 1, 1, 1], 4, 4), sparsevec([2], [1.0e-9], 2), sparsevec(Int64[], Nothing[], 2), sparsevec([3, 4], [1.0e-9, 1.0e-9], 4), 2, 2, ComplexF64[0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im]), ComplexF64[0.9789325751107439 + 0.0im; 0.12051187816111208 + 0.07881274759778019im; -0.004230418875131352 - 0.00964520752500062im; 0.00019110015334431832 + 0.0im;;], JosephsonCircuits.LinearHB([0.9999636971626665 + 0.026577597828732157im 0.024776280908444095 + 0.004460970065009537im; 0.02006351652496576 - 0.01520590853015591im 0.8793137132801259 - 0.4769079135295333im;;; 0.8793041537769276 + 0.4769257467844491im 0.020066345817057213 + 0.015208700168465511im; 0.0247802510228797 - 0.004461160668153754im 0.9999643476849712 - 0.02655684759593593im], Array{ComplexF64, 3}(undef, 0, 0, 0), [0.9993670379479913 0.0006329620520086488; 0.0006329620520086476 0.9993670379479913;;; 0.9993668400044972 0.0006331599955028162; 0.0006331599955028159 0.9993668400044972], [0.9993670379479914 1.0; 1.0 0.9993670379479901;;; 0.9993668400044982 1.0; 1.0 0.9993668400044973], [-0.9999999999999998 -0.9999999999999991; 1.000000000000001 0.9999999999999998], ComplexF64[], ComplexF64[], 2, 3, 2, 2, StepRangeLen(Base.TwicePrecision{Float64}(2.827433388230814e10, -1.7570437194081023e-6),Base.TwicePrecision{Float64}(3.141592653589793e9, -8.926326700020581e-8),2)))
    #     out2 = JosephsonCircuits.warmupsyms()
    #     @test JosephsonCircuits.compare(out1,out2)
    # end

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
        out1 = JosephsonCircuits.LinearHB([0.8952708641229391 - 0.44552225517090227im;;; 0.8415115570832487 - 0.5402391130743189im;;; 0.6457820691998714 - 0.7635217869189669im;;; -0.9968560060568034 + 0.07923448231975308im;;; 0.9316787544566122 + 0.36328322077158454im;;; 0.9988570509555925 + 0.04779740323801577im], Array{Complex{Float64}, 3}(undef, 0, 0, 0), [1.0;;; 1.0;;; 1.0;;; 1.0;;; 1.0;;; 1.0], [0.9999999999999993;;; 0.9999999999999996;;; 1.0;;; 0.9999999999999991;;; 1.0;;; 0.9999999999999993], [1.0000000000000007 1.0000000000000004 0.9999999999999997 1.0000000000000009 1.0 1.0000000000000007], ComplexF64[], ComplexF64[], 1, 3, 2, 1, StepRangeLen(Base.TwicePrecision{Float64}(2.827433388230814e10, -1.7570437194081023e-6),Base.TwicePrecision{Float64}(6.283185307179585e8, 1.7288220988120884e-7),6))
        out2 = JosephsonCircuits.warmuphblinsolve()
        @test JosephsonCircuits.compare(out1,out2)
    end

    @testset verbose=true "warmupvvn" begin
        out1 = Real[1, 1.0e-8, 50.0, 1.0e-13, 1.0e-9, 1.0e-12]
        out2 = JosephsonCircuits.warmupvvn()
        @test JosephsonCircuits.compare(out1,out2)
    end

end