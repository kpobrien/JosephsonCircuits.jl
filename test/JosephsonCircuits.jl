using JosephsonCircuits
using SparseArrays

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
    # the occupations of this lossless circuit in its vacuum are rounding,
    # which a relative comparison cannot hold, so they are not asked for
    return JosephsonCircuits.hblinsolve(2*pi*(4.5:0.1:5.0)*1e9, circuit,
        JosephsonCircuits.warmupdefs(Rleft, Cc, Lj, Cj); returnnbar = false)
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
            atol = 1e-12)
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
        # the port owns its reference impedance through the termination
        # generated for it
        @test out.componentnames[only(out.ports).environment] ==
            "P1/termination"
    end

    @testset verbose=true "warmupnumericmatrices" begin
        out1 = JosephsonCircuits.CircuitMatrices(sparse([1, 2, 1, 2], [1, 1, 2, 2], [1.0e-13, -1.0e-13, -1.0e-13, 1.1e-12], 2, 2), sparse([1], [1], [0.02], 2, 2), sparsevec(Int64[], Float64[], 2), sparsevec([2], [1.0e-9], 2), sparsevec([2], [1.0e-9], 2), sparse(Int64[], Int64[], Float64[], 2, 2), sparse(Int64[], Int64[], Float64[], 2, 2), sparse([1, 2], [1, 2], [1, 1], 2, 2), [1], [1], [50.0], [2], 1.0e-9, [50.0, 50.0, 1.0e-13, 1.0e-9, 1.0e-12])
        out2 = warmupnumericmatrices()
        @test JosephsonCircuits.compare(out1,out2)

    end

    @testset verbose=true "warmuphblinsolve" begin
        # unpumped, the amplifier is its coupling capacitor in series with
        # the inductance and the capacitance of its junction in parallel,
        # which the matched port sees as the reflection of that impedance,
        # lossless: the quantum efficiency and the commutation relations
        # of a lossless reflection are one
        out = warmuphblinsolve()
        Cc, Lj, Cj, Z0 = 100e-15, 1000e-12, 1000e-15, 50.0
        Z(w) = 1/(im*w*Cc) + 1/(im*w*Cj + 1/(im*w*Lj))
        @test vec(Array(out.S)) ≈ [(Z(w) - Z0)/(Z(w) + Z0) for w in out.w] rtol = 1e-12
        @test all(x -> isapprox(x, 1; atol = 1e-12), out.QE)
        @test all(x -> isapprox(x, 1; atol = 1e-12), out.CM)
    end

    @testset verbose=true "warmupvvn" begin
        # the port's slot holds its reference impedance, so every entry is
        # a quantity, compared by value; the vector is a `Vector{Any}`, as
        # componentvaluestonumber returns it
        out1 = [50.0, 50.0, 1.0e-13, 1.0e-9, 1.0e-12]
        out2 = warmupvvn()
        @test JosephsonCircuits.compare(out1,out2)
    end

    @testset verbose=true "warmupconnect" begin
        @test JosephsonCircuits.warmupconnect()
    end

end