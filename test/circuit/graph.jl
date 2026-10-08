using Random
using JosephsonCircuits
using Test

@testset verbose=true "graphproc" begin
    # the branches the incidence matrix is built from
    @testset "extractbranches" begin
        @test_throws(
            DimensionMismatch("componenttypes must have the same length as the number of node indices"),
            JosephsonCircuits.extractbranches(
                [:P,:I,:R,:C,:Lj,:C],
                [2 2 2 2 3; 1 1 1 3 1]
            )
        )

        @test_throws(
            DimensionMismatch("the length of the first axis must be 2"),
            JosephsonCircuits.extractbranches(
                [:P,:I,:R,:C,:Lj,:C],
                [2 2 2 2 3 3; 1 1 1 3 1 1; 0 0 0 0 0 0],
            )
        )
    end
end

@testset "the incidence matrix of the branches" begin
    # written from its definition, without the graph library: the distinct
    # branches in ascending order of their endpoints, each row -1 at its
    # lower node and +1 at its higher, ground (node 1) having no column
    # and a self loop no entries, and either orientation of a branch naming
    # its row
    rng = Random.MersenneTwister(190926)
    JC = JosephsonCircuits
    for trial in 1:20
        n = rand(rng, 1:40)
        edges = [(rand(rng,1:n), rand(rng,1:n)) for _ in 1:rand(rng,0:80)]
        t = JC.circuittopology(edges, n)
        branches = sort(unique([minmax(e...) for e in edges]))
        Rbn = zeros(Int, length(branches), n - 1)
        for (i, (a, b)) in enumerate(branches)
            a == b && continue
            a > 1 && (Rbn[i, a-1] = -1)
            b > 1 && (Rbn[i, b-1] = 1)
        end
        row(e) = findfirst(==(minmax(e...)), branches)
        @test t.Nbranches == length(branches) && Matrix(t.Rbn) == Rbn &&
            all(t.edge2indexdict[e] == t.edge2indexdict[reverse(e)] == row(e)
                for e in edges)
    end
end

@testset "components with both terminals on one node" begin
    # such a component carries no current: the response of the circuit, its
    # noise included, is the one without it, and a port so placed is refused
    ws = 2pi*(4.5:0.5:5.5)*1e9
    sol(c) = hblinsolve(ws, Circuit(c); keyedarrays = false, returnSnoise = true)
    # node 2 has no inductive branch, so a self loop there is its only one
    base = [(:p1, 1, 0, Port(1)), (:c1, 1, 2, Capacitor(100e-15)),
        (:c2, 2, 0, Capacitor(1e-12)), (:r2, 2, 0, Resistor(1e4))]
    ref = sol(base)
    # and on node 5, which nothing else touches
    for extra in ((:cx, 0, 0, Capacitor(1e-12)), (:rx, 0, 0, Resistor(50.0)),
            (:rx, 2, 2, Resistor(50.0)), (:lx, 0, 0, Inductor(1e-9)),
            (:lx, 2, 2, Inductor(1e-9)), (:jx, 2, 2, JosephsonJunction(1e-9)),
            (:cx, 5, 5, Capacitor(1e-12)), (:lx, 5, 5, Inductor(1e-9)),
            (:jx, 5, 5, JosephsonJunction(1e-9)))
        s = sol(vcat(base, [extra]))
        @test s.S ≈ ref.S && s.QE ≈ ref.QE && size(s.Snoise) == size(ref.Snoise)
    end
    @test_throws ArgumentError sol([(:p1, 1, 1, Port(1)), (:c1, 1, 0, Capacitor(1e-12))])
    # on a node of its own, a secondary coupled to a primary shorts it, which
    # leaves the primary L1(1 - K^2), and a line whose far port is shorted is
    # a stub of input impedance i Z0 tan(w l/vp)
    L1, L2, K = 1e-9, 2e-9, 0.6
    shorted = sol([(:p1, 1, 0, Port(1)), (:l1, 1, 0, Inductor(L1)),
        (:l2, 3, 3, Inductor(L2)), (:k, :l1, :l2, MutualInductor(K)),
        (:c1, 1, 0, Capacitor(1e-12))])
    reduced = sol([(:p1, 1, 0, Port(1)), (:l1, 1, 0, Inductor(L1*(1 - K^2))),
        (:c1, 1, 0, Capacitor(1e-12))])
    @test shorted.S ≈ reduced.S
    Z0, len = 50.0, 0.004
    stub = sol([(:p1, 1, 0, Port(1)),
        (:line, 1, 0, 2, 2, TransmissionLine(Z0, len; grounded = false))])
    zin = [im*Z0*tan(w*len/JosephsonCircuits.speed_of_light) for w in ws]
    @test vec(stub.S) ≈ (zin .- 50)./(zin .+ 50)
end

@testset "a subnetwork no element connects to ground" begin
    # harmonic balance cannot determine its potential, so both solvers
    # refuse it by name; the transient fixes that potential and solves it
    # (transient/solve.jl)
    c = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:c1, 1, 0, Capacitor(1e-12)),
        (:c2, 2, 3, Capacitor(1e-12)), (:l2, 2, 3, Inductor(1e-9))])
    psc = compile(c)
    @test [psc.nodenames[i] for i in only(JosephsonCircuits.isolatedsubnetworks(psc))] ==
        ["2", "3"]
    @test_throws ArgumentError hblinsolve(2pi*[4e9, 5e9], c)
    @test_throws ArgumentError hbnlsolve((2pi*5e9,), (4,),
        [(mode = (1,), port = 1, current = 1e-8)], c)
end
