using JosephsonCircuits
using SparseArrays
import AxisKeys
import AxisKeys.NamedDims
import AxisKeys.NamedDims: NamedDimsArray

using Test

# The fixtures of the warmup circuit which are not part of the precompile
# workload: the same single junction amplifier as `warmup`, taken through
# the compile, the numeric matrices, a pump free sweep, the value
# resolution and every network conversion.
# Elaboration and compilation of the typed circuit alone.
function warmupcompile()

    JosephsonCircuits.@params Rleft Cc Lj Cj
    return JosephsonCircuits.compile(JosephsonCircuits.warmupcircuit(Rleft, Cc, Lj, Cj))
end

# The capacitance and inverse inductance matrices at numeric values.
function warmupnumericmatrices()

    JosephsonCircuits.@params Rleft Cc Lj Cj
    circuit = JosephsonCircuits.warmupcircuit(Rleft, Cc, Lj, Cj)
    return JosephsonCircuits.numericmatrices(circuit, JosephsonCircuits.warmupdefs(Rleft, Cc, Lj, Cj))
end

# A linear (no pump) frequency sweep.
function warmuphblinsolve()

    JosephsonCircuits.@params Rleft Cc Lj Cj
    circuit = JosephsonCircuits.warmupcircuit(Rleft, Cc, Lj, Cj)
    return JosephsonCircuits.hblinsolve(2*pi*(4.5:0.1:5.0)*1e9, circuit,
        JosephsonCircuits.warmupdefs(Rleft, Cc, Lj, Cj))
end

# Resolving the compiled component values to numbers.
function warmupvvn()

    JosephsonCircuits.@params Rleft Cc Lj Cj
    psc = JosephsonCircuits.compile(JosephsonCircuits.warmupcircuit(Rleft, Cc, Lj, Cj); sorting = :number)

    return JosephsonCircuits.componentvaluestonumber(psc.componentvalues,
        JosephsonCircuits.warmupdefs(Rleft, Cc, Lj, Cj))
end



@testset verbose=true "JosephsonCircuits" begin

    @testset verbose=true "warmup and warmupsyms" begin
        # the precompile workloads solve the warmup circuit with symbol
        # valued and with @params valued components; both must give the
        # numbers a plain numeric circuit gives
        circuit = Circuit(
            ["P1" => Port(1; Z0 = 50.0), "C1" => Capacitor(100.0e-15),
             "Lj1" => JosephsonJunction(1000.0e-12), "C2" => Capacitor(1000.0e-15)],
            [Net("1", [("P1",1), ("C1",1)]),
             Net("2", [("C1",2), ("Lj1",1), ("C2",1)]),
             Net("0", [("P1",2), ("Lj1",2), ("C2",2), Ground])])
        ref = hbsolve(2*pi*(4.5:0.5:5.0)*1e9, (2*pi*4.75001*1e9,),
            [(mode=(1,),port=1,current=0.00565e-6)], (2,), (4,), circuit;
            ftol = 1e-12)
        @test JosephsonCircuits.compare(ref, JosephsonCircuits.warmup())
        @test JosephsonCircuits.compare(ref, JosephsonCircuits.warmupsyms())
    end

    @testset verbose=true "warmupcompile" begin
        out = warmupcompile()
        @test out.componentnames ==
            ["P1", "P1/termination", "C1", "Lj1", "C2"]
        @test out.componenttypes == [:P, :R, :C, :Lj, :C]
        @test out.nodenames == ["0", "1", "2"]
        @test out.nodeindices == [2 2 2 3 3; 1 1 3 1 1]
        @test out.Nnodes == 3
        # the port owns the reference impedance the legacy netlist wrote as a
        # separate resistor, so the node layout is the legacy netlist's
        @test out.componentnames[only(out.ports).environment] ==
            "P1/termination"
    end

    @testset verbose=true "warmupnumericmatrices" begin
        out1 = JosephsonCircuits.CircuitMatrices(sparse([1, 2, 1, 2], [1, 1, 2, 2], [1.0e-13, -1.0e-13, -1.0e-13, 1.1e-12], 2, 2), sparse([1], [1], [0.02], 2, 2), sparsevec(Int64[], Float64[], 2), sparsevec(Int64[], Float64[], 2), sparsevec([2], [1.0e-9], 2), sparsevec([2], [1.0e-9], 2), sparse(Int64[], Int64[], Float64[], 2, 2), sparse(Int64[], Int64[], Float64[], 2, 2), sparse([1, 2], [1, 2], [1, 1], 2, 2), [1], [1], [50.0], [2], Int64[], 1.0e-9, [50.0, 50.0, 1.0e-13, 1.0e-9, 1.0e-12])
        out2 = warmupnumericmatrices()
        @test JosephsonCircuits.compare(out1,out2)

    end

    @testset verbose=true "warmuphblinsolve" begin
        # JosephsonCircuits.testshow(stdout,warmuphblinsolve())
        out1 = JosephsonCircuits.LinearizedHB(collect(2*pi*(4.5:0.1:5.0)*1e9), [(0,)], AxisKeys.KeyedArray(NamedDimsArray(ComplexF64[0.895270864122939 - 0.4455222551709022im;;;;; 0.8415115570832487 - 0.5402391130743189im;;;;; 0.6457820691998714 - 0.7635217869189669im;;;;; -0.9968560060568034 + 0.07923448231975308im;;;;; 0.9316787544566122 + 0.36328322077158454im;;;;; 0.9988570509555925 + 0.04779740323801577im], (:outputmode, :outputport, :inputmode, :inputport, :freqindex)), ([(0,)], [1], [(0,)], [1], 1:6)), Array{ComplexF64, 3}(undef, 0, 0, 0), Array{ComplexF64, 3}(undef, 0, 0, 0), Array{ComplexF64, 4}(undef, 0, 0, 0, 0), AxisKeys.KeyedArray(NamedDimsArray([1.0;;;;; 1.0;;;;; 1.0;;;;; 1.0;;;;; 1.0;;;;; 1.0], (:outputmode, :outputport, :inputmode, :inputport, :freqindex)), ([(0,)], [1], [(0,)], [1], 1:6)), AxisKeys.KeyedArray(NamedDimsArray([0.9999999999999996;;;;; 0.9999999999999996;;;;; 1.0;;;;; 0.9999999999999991;;;;; 1.0;;;;; 0.9999999999999993], (:outputmode, :outputport, :inputmode, :inputport, :freqindex)), ([(0,)], [1], [(0,)], [1], 1:6)), AxisKeys.KeyedArray(NamedDimsArray([1.0000000000000004;;; 1.0000000000000004;;; 0.9999999999999997;;; 1.0000000000000009;;; 1.0;;; 1.0000000000000007], (:outputmode, :outputport, :freqindex)), ([(0,)], [1], 1:6)), Array{ComplexF64, 3}(undef, 0, 0, 0), Array{ComplexF64, 3}(undef, 0, 0, 0), Array{ComplexF64, 3}(undef, 0, 0, 0), Array{ComplexF64, 3}(undef, 0, 0, 0), ["0", "1", "2"], [2 2 2 3 3; 1 1 3 1 1], ["P1", "P1/termination", "C1", "Lj1", "C2"], [:P, :R, :C, :Lj, :C], Dict("C1" => 3, "C2" => 5, "P1/termination" => 2, "P1" => 1, "Lj1" => 4), String[], [1], [1], [50.0], Int64[], String[], Int64[], 1, 3, 2, 1, 1)
        out2 = warmuphblinsolve()
        @test JosephsonCircuits.compare(out1,out2)
    end

    @testset verbose=true "warmupvvn" begin
        # the port's slot holds its reference impedance, so every entry is
        # a quantity and the vector is concretely typed
        out1 = [50.0, 50.0, 1.0e-13, 1.0e-9, 1.0e-12]
        out2 = warmupvvn()
        @test JosephsonCircuits.compare(out1,out2)
    end

    @testset verbose=true "warmupconnect" begin
        @test JosephsonCircuits.warmupconnect()
    end

end