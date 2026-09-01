using Symbolics
using JosephsonCircuits
using LinearAlgebra
using Test

@testset verbose=true "circuit values" begin

    @testset "Symbolics.jl and complex numbers" begin
    # this functionality is currently broken. this test will error once it is
    # fixed and then we can revise valuetonumber

        @variables Lj1
        @test(Symbolics.evaluate(Lj1,Dict(Lj1=>3.0e-12))==3.0e-12)

        @variables Lj1::Complex
        @test(Symbolics.evaluate(Lj1,Dict(Lj1=>3.0e-12))==3.0e-12,broken=true)

        @variables Lj1 Lj2
        @test(Symbolics.evaluate(Lj1+Lj2,Dict(Lj1=>3.0e-12,Lj2=>1.0e-12)) == 4.0e-12)

        @variables Lj1::Complex Lj2::Complex
        @test(Symbolics.evaluate(Lj1+Lj2,Dict(Lj1=>3.0e-12,Lj2=>1.0e-12)) == 4.0e-12,broken=true)
    end
end
