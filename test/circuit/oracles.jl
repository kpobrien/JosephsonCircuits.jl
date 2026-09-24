using JosephsonCircuits
using LinearAlgebra
using Random
using Test

# An independent reference for the linear solve, and the invariance of a
# circuit's response to the way its description is written.
#
# The reference is the nodal admittance of a two port written out by hand,
# Y(w) = i w C + G + inv(L)/(i w), and its scattering matrix at a common
# reference impedance, (I - Z0 Y)/(I + Z0 Y). It shares no code with the
# package: no compilation, no branch map, no stamps, no augmentation. The
# same circuit is then written five ways -- in the order it was built, in a
# shuffled order, with the nodes renamed and the two terminal branches
# reversed (a reversed inductor reverses the sign of its mutual coupling),
# as a subcircuit of a hierarchy, and already compiled -- and every one of
# them must give the reference.
@testset "an independent nodal oracle, and the frontends agree with it" begin
    rng = MersenneTwister(72841)
    ws = 2*pi*[1.7e9, 4.3e9, 7.1e9]
    Z0 = 50.0
    for _ in 1:12
        L1, L2 = (0.5 .+ rand(rng, 2))*1e-9
        # |k| < 1 keeps the inductance matrix positive definite
        k = 1.4*rand(rng) - 0.7
        C1, C2, Cc = (0.1 .+ rand(rng, 3))*1e-12
        R1, R2 = 100 .+ 300*rand(rng, 2)
        tand = 0.01*rand(rng)

        C = [C1/(1 + im*tand) + Cc  -Cc; -Cc  C2 + Cc]
        G = Diagonal([1/R1, 1/R2])
        L = [L1  k*sqrt(L1*L2); k*sqrt(L1*L2)  L2]
        reference = cat([(Y = im*w*C + G + inv(L)/(im*w);
            (I - Z0*Y)/(I + Z0*Y)) for w in ws]...; dims = 3)

        entries = Any[
            (:p1, 1, 0, Port(1; Z0 = Z0)), (:p2, 2, 0, Port(2; Z0 = Z0)),
            (:c1, 1, 0, Capacitor(C1/(1 + im*tand))),
            (:c2, 2, 0, Capacitor(C2)), (:cc, 1, 2, Capacitor(Cc)),
            (:r1, 1, 0, Resistor(R1)), (:r2, 2, 0, Resistor(R2)),
            (:l1, 1, 0, Inductor(L1)), (:l2, 2, 0, Inductor(L2)),
            (:k, :l1, :l2, MutualInductor(k))]
        # every branch but the second inductor and the two ports is written
        # the other way round, and the first inductor's reversal reverses
        # the sign of the coupling between the two
        reversed = map(entries) do (name, n1, n2, component)
            name === :k && return (name, n1, n2, MutualInductor(-k))
            node(n) = n == 0 ? "0" : n == 1 ? "z_input" : "a_output"
            name in (:l1, :c1, :c2, :cc, :r1, :r2) ?
                (name, node(n2), node(n1), component) :
                (name, node(n1), node(n2), component)
        end
        nested = Circuit([:device => Circuit(reversed; pins = [1 => (:p1, 1)])],
            [])
        original = Circuit(entries)
        for c in (original, Circuit(shuffle(rng, entries)), Circuit(reversed),
                nested, compile(original))
            S = Array(hblinsolve(ws, c; keyedarrays = false).S)
            @test reshape(S, 2, 2, :) ≈ reference rtol = 1e-11 atol = 1e-13
        end
    end
end
