using JosephsonCircuits
using LinearAlgebra
using SparseArrays
using Test

# A component value given as a closure of the mode frequency: the
# provider the circuit value carries, through a solve.
@testset verbose=true "circuit values" begin

    wp = (2*pi*4.75001*1e9,)
    src = [(mode=(1,), port=1, current=0.00565e-6)]
    ws = 2*pi*(4.5:0.1:5.0)*1e9

    @testset "FrequencyDependent provider" begin
        # a lossy frequency dependent resistor with an arbitrary law (trig
        # inside the closure -- the generality a closed expression set
        # cannot offer)
        law = w -> 50.0*(1 + 0.1*sin(w/2e10)) + im*abs(w)*1e-9
        circuit = Tuple{String,String,String,Any}[
            ("P1","1","0",1), ("R1","1","0", FrequencyDependent(law)),
            ("C1","1","2",100e-15), ("Lj1","2","0",1000e-12),
            ("C2","2","0",1000e-15)]
        out = hbnlsolve(wp, (8,), src, circuit;
            keyedarrays = false)
        @test out.solverinfo.converged
        o2 = hbsolve(ws, wp, src, (2,), (8,), circuit)
        @test all(isfinite, Array(o2.linearized.S))

        # a CONSTANT law must agree exactly with a plain numeric value
        circa = Tuple{String,String,String,Any}[
            ("P1","1","0",1), ("R1","1","0", FrequencyDependent(w -> 50.0)),
            ("C1","1","2",100e-15), ("Lj1","2","0",1000e-12),
            ("C2","2","0",1000e-15)]
        circb = Tuple{String,String,String,Any}[
            ("P1","1","0",1), ("R1","1","0", 50.0),
            ("C1","1","2",100e-15), ("Lj1","2","0",1000e-12),
            ("C2","2","0",1000e-15)]
        Sa = hbsolve(ws, wp, src, (2,), (8,), circa).linearized.S
        Sb = hbsolve(ws, wp, src, (2,), (8,), circb).linearized.S
        @test isapprox(Array(Sa), Array(Sb), rtol = 1e-12)

        # providers combine with numbers through the value arithmetic
        circc = Tuple{String,String,String,Any}[
            ("P1","1","0",1), ("R1","1","0", 2*FrequencyDependent(w -> 25.0)),
            ("C1","1","2",100e-15), ("Lj1","2","0",1000e-12),
            ("C2","2","0",1000e-15)]
        Sc = hbsolve(ws, wp, src, (2,), (8,), circc).linearized.S
        @test isapprox(Array(Sc), Array(Sb), rtol = 1e-12)

        # the linear-only entry point also accepts a circuit alone
        lin = hblinsolve(ws, circb; Nmodulationharmonics = (2,))
        @test all(isfinite, Array(lin.S))
    end
end
