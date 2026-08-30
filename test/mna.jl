using Symbolics
using JosephsonCircuits
using Symbolics: @variables, Num
using Test
isdefined(Main, :testjpacircuit) || include(joinpath(@__DIR__, "testcircuits.jl"))

@testset verbose=true "mna" begin

    # a JPA-like circuit used by several testsets below
    # the canonical JPA lives in testcircuits.jl, shared across test files
    jpacircuit() = testjpacircuit()

    @testset "hbnlsolve nodal reference values" begin

        # reference values computed with the purely nodal formulation
        # (which this formulation replaces) at ftol = 1e-12; the two
        # formulations were verified to agree to ~1e-11 before the nodal
        # path was removed.
        circuit, circuitdefs = jpacircuit()
        wp = (2*pi*4.75001*1e9,)
        sources = [(mode=(1,),port=1,current=0.00565e-6)]

        out = hbnlsolve(wp, (16,), sources, circuit, circuitdefs;
            ftol = 1e-12)
        @test out.solverinfo.converged
        @test isapprox(Vector(out.nodeflux(outputmode=(1,))),
            ComplexF64[-0.013189575461243642 - 0.008650771565475453im,
                0.12157593799903155 + 0.07973582686023277im], atol = 1e-8)
        @test isapprox(
            out.S(outputmode=(1,),outputport=1,inputmode=(1,),inputport=1),
            -0.3984433276327756 - 0.9171756605002952im, atol = 1e-8)
    end

    @testset "hbnlsolve nodal reference values two tone" begin

        # nodal formulation reference values, as above
        circuit, circuitdefs = jpacircuit()
        wp = (2*pi*4.75001*1e9, 2*pi*4.35001*1e9)
        sources = [
            (mode=(1,0),port=1,current=0.00565e-6),
            (mode=(0,1),port=1,current=0.00265e-6),
        ]

        out = hbnlsolve(wp, (8,8), sources, circuit, circuitdefs;
            ftol = 1e-12)
        @test out.solverinfo.converged
        @test isapprox(Vector(out.nodeflux(outputmode=(1,0))),
            ComplexF64[-0.013075460477876147 - 0.008395813901112565im,
                0.1233985871775909 + 0.07922607187414707im], atol = 1e-8)
        @test isapprox(Vector(out.nodeflux(outputmode=(0,1))),
            ComplexF64[-0.0028520505572584252 - 0.014156019369081027im,
                0.00134982177864443 + 0.006713723189472766im], atol = 1e-8)
    end

    @testset "hbnlsolve mna dc gauge fixing isolated node" begin

        # node 1 is connected only by a resistor and a capacitor, so with
        # dc = true its DC node flux is unconstrained and the nodal system
        # matrix is structurally singular. the mna formulation adds a gauge
        # fixing equation instead of requiring a workaround inductor.
        circuit, circuitdefs = jpacircuit()
        wp = (2*pi*4.75001*1e9,)
        sources = [(mode=(1,),port=1,current=0.00565e-6)]

        mnadc = hbnlsolve(wp, (16,), sources, circuit, circuitdefs;
            ftol = 1e-12, dc = true, odd = true)
        @test mnadc.solverinfo.converged

        # the DC node flux of the gauge fixed node should be zero
        @test isapprox(abs(mnadc.nodeflux(outputmode=(0,),node="1")), 0.0,
            atol = 1e-12)

        # with no DC source, including the DC mode should not change the
        # nonzero frequency modes
        mnanodc = hbnlsolve(wp, (16,), sources, circuit, circuitdefs;
            ftol = 1e-12)
        @test mnanodc.solverinfo.converged
        for m in mnanodc.modes
            for node in ["1", "2"]
                @test isapprox(
                    mnadc.nodeflux(outputmode=m, node=node),
                    mnanodc.nodeflux(outputmode=m, node=node),
                    atol = 1e-9)
            end
        end
    end

    @testset "hbnlsolve mna floating inductor island" begin

        # nodes 2 and 3 are joined by an inductor but have no inductive
        # path to ground, so their common DC flux is a gauge degree of
        # freedom even though every column of the inverse inductance matrix
        # is structurally nonempty. with dc = true the nodal system matrix
        # is singular; the mna formulation adds one gauge equation for the
        # component and solves exactly.
        circuit = Tuple{String,String,String,Union{Complex{Float64},Symbol,Int64}}[]
        push!(circuit,("P1","1","0",1))
        push!(circuit,("R1","1","0",:Rleft))
        push!(circuit,("Lj1","1","0",:Lj))
        push!(circuit,("C1","1","2",:Cc))
        push!(circuit,("L1","2","3",:L1))
        push!(circuit,("C2","3","0",:Cj))
        circuitdefs = Dict(
            :Lj =>1000.0e-12,
            :Cc => 100.0e-15,
            :Cj => 1000.0e-15,
            :Rleft => 50.0,
            :L1 => 300.0e-12,
        )
        wp = (2*pi*4.75001*1e9,)
        sources = [(mode=(1,),port=1,current=0.00565e-6)]

        mnadc = hbnlsolve(wp, (8,), sources, circuit, circuitdefs;
            ftol = 1e-12, dc = true, odd = true)
        @test mnadc.solverinfo.converged

        # the DC flux of the reference node of the floating component is
        # gauge fixed to zero
        @test isapprox(abs(mnadc.nodeflux(outputmode=(0,),node="2")), 0.0,
            atol = 1e-12)

        # with no DC source, the nonzero frequency modes are unchanged by
        # including the DC mode
        mnanodc = hbnlsolve(wp, (8,), sources, circuit, circuitdefs;
            ftol = 1e-12)
        @test mnanodc.solverinfo.converged
        for m in mnanodc.modes
            for node in ["1", "2", "3"]
                @test isapprox(
                    mnadc.nodeflux(outputmode=m, node=node),
                    mnanodc.nodeflux(outputmode=m, node=node),
                    atol = 1e-9)
            end
        end
    end

    @testset "hbnlsolve mna floating Josephson island" begin

        # a junction joins the floating nodes 2 and 3: the junction
        # constrains only the difference of their DC fluxes, so the common
        # mode is a gauge degree of freedom, and node 1 forms a second
        # floating component on its own. one gauge equation per component
        # is added.
        circuit = Tuple{String,String,String,Union{Complex{Float64},Symbol,Int64}}[]
        push!(circuit,("P1","1","0",1))
        push!(circuit,("R1","1","0",:Rleft))
        push!(circuit,("C1","1","2",:Cc))
        push!(circuit,("Lj1","2","3",:Lj))
        push!(circuit,("C2","3","0",:Cj))
        circuitdefs = Dict(
            :Lj =>1000.0e-12,
            :Cc => 100.0e-15,
            :Cj => 1000.0e-15,
            :Rleft => 50.0,
        )
        wp = (2*pi*4.75001*1e9,)
        sources = [(mode=(1,),port=1,current=0.00565e-6)]

        mnadc = hbnlsolve(wp, (8,), sources, circuit, circuitdefs;
            ftol = 1e-12, dc = true, odd = true)
        @test mnadc.solverinfo.converged

        # both floating components are gauge fixed at their reference nodes
        @test isapprox(abs(mnadc.nodeflux(outputmode=(0,),node="1")), 0.0,
            atol = 1e-12)
        @test isapprox(abs(mnadc.nodeflux(outputmode=(0,),node="2")), 0.0,
            atol = 1e-12)

        mnanodc = hbnlsolve(wp, (8,), sources, circuit, circuitdefs;
            ftol = 1e-12)
        @test mnanodc.solverinfo.converged
        for m in mnanodc.modes
            for node in ["1", "2", "3"]
                @test isapprox(
                    mnadc.nodeflux(outputmode=m, node=node),
                    mnanodc.nodeflux(outputmode=m, node=node),
                    atol = 1e-9)
            end
        end
    end

    @testset "hbnlsolve mna balanced dc excitation" begin

        # a DC current driven between the two nodes of a floating inductive
        # island is balanced (zero net injection), so a periodic solution
        # exists: the DC current flows through the island inductor and its
        # branch flux is L1*Idc.
        circuit = Tuple{String,String,String,Union{Complex{Float64},Symbol,Int64}}[]
        push!(circuit,("P1","1","0",1))
        push!(circuit,("R1","1","0",:Rleft))
        push!(circuit,("Lj1","1","0",:Lj))
        push!(circuit,("C1","1","2",:Cc))
        push!(circuit,("L1","2","3",:L1))
        push!(circuit,("P2","2","3",2))
        push!(circuit,("R2","2","3",:Rright))
        push!(circuit,("C2","3","0",:Cj))
        circuitdefs = Dict(
            :Lj =>1000.0e-12,
            :Cc => 100.0e-15,
            :Cj => 1000.0e-15,
            :Rleft => 50.0,
            :Rright => 50.0,
            :L1 => 1000.0e-12,
        )
        wp = (2*pi*4.75001*1e9,)
        Idc = 1.0e-5
        sources = [
            (mode=(1,),port=1,current=0.00565e-6),
            (mode=(0,),port=2,current=Idc),
        ]

        mnadc = hbnlsolve(wp, (8,), sources, circuit, circuitdefs;
            ftol = 1e-12, dc = true, odd = true)
        @test mnadc.solverinfo.converged

        # the DC branch flux of the island inductor equals L1*Idc. the
        # node fluxes are in units of the reduced flux quantum phi0.
        dphi = abs(mnadc.nodeflux(outputmode=(0,),node="2") -
            mnadc.nodeflux(outputmode=(0,),node="3"))
        @test isapprox(dphi, circuitdefs[:L1]*Idc/JosephsonCircuits.phi0,
            rtol = 1e-6)
    end

    @testset "hbnlsolve mna unbalanced dc excitation throws" begin

        # a nonzero net DC current into a floating component has no
        # periodic solution. the compatibility check must throw instead of
        # letting the gauge equation absorb the source into the flux
        # reference and report an apparently converged nonphysical result.
        circuit, circuitdefs = jpacircuit()
        wp = (2*pi*4.75001*1e9,)
        sources = [
            (mode=(0,),port=1,current=1.0e-7),
            (mode=(1,),port=1,current=0.00565e-6),
        ]

        @test_throws ArgumentError hbnlsolve(wp, (8,), sources, circuit,
            circuitdefs; ftol = 1e-12, dc = true, odd = true)
    end

    @testset "hbnlsolve mna complex storage of real resistors" begin

        # the package's own examples store circuit definitions in
        # containers such as Dict{Symbol,Complex{Float64}}, so a real
        # resistance can arrive as 50.0 + 0.0im. it must still be promoted.
        circuit, _ = jpacircuit()
        circuitdefs = Dict{Symbol,Complex{Float64}}(
            :Lj =>1000.0e-12,
            :Cc => 100.0e-15,
            :Cj => 1000.0e-15,
            :Rleft => 50.0,
        )
        wp = (2*pi*4.75001*1e9,)
        sources = [(mode=(1,),port=1,current=0.00565e-6)]

        out = hbnlsolve(wp, (16,), sources, circuit, circuitdefs;
            ftol = 1e-12)
        @test out.solverinfo.converged
        # the result must match the real-typed definitions (the nodal
        # formulation reference values from the first testset)
        @test isapprox(Vector(out.nodeflux(outputmode=(1,))),
            ComplexF64[-0.013189575461243642 - 0.008650771565475453im,
                0.12157593799903155 + 0.07973582686023277im], atol = 1e-8)

        # and the DC gauge case with complex-typed definitions
        outdc = hbnlsolve(wp, (16,), sources, circuit, circuitdefs;
            ftol = 1e-12, dc = true, odd = true)
        @test outdc.solverinfo.converged
    end

    @testset "hbnlsolve mna mixed numeric and symbolic resistors" begin

        # the numeric port resistor is promoted while the symbolic
        # frequency dependent interior resistor is left in the conductance
        # matrix. this exercises the partial promotion path where the
        # promoted subset is subtracted from Gnm.
        @variables w
        circuit = Tuple{String,String,String,Union{Complex{Float64},Symbol,Int64,Num}}[]
        push!(circuit,("P1","1","0",1))
        push!(circuit,("R1","1","0",:Rleft))
        push!(circuit,("C1","1","2",:Cc))
        push!(circuit,("Lj1","2","0",:Lj))
        push!(circuit,("C2","2","0",:Cj))
        push!(circuit,("R2","2","0",1.0e6 + 1e-3*w))
        circuitdefs = Dict(
            :Lj =>1000.0e-12,
            :Cc => 100.0e-15,
            :Cj => 1000.0e-15,
            :Rleft => 50.0,
        )
        wp = (2*pi*4.75001*1e9,)
        sources = [(mode=(1,),port=1,current=0.00565e-6)]

        out = hbnlsolve(wp, (16,), sources, circuit, circuitdefs;
            ftol = 1e-12, symfreqvar = w)
        @test out.solverinfo.converged
        # nodal formulation reference values, as in the first testset
        @test isapprox(Vector(out.nodeflux(outputmode=(1,))),
            ComplexF64[-0.013190622855691949 - 0.00865617332363941im,
                0.12153869209076193 + 0.07973744395061082im], atol = 1e-8)
    end

    @testset "hbnlsolve mna warm start" begin

        circuit, circuitdefs = jpacircuit()
        wp = (2*pi*4.75001*1e9,)
        sources = [(mode=(1,),port=1,current=0.00565e-6)]

        out = hbnlsolve(wp, (16,), sources, circuit, circuitdefs;
            ftol = 1e-12)
        @test out.solverinfo.converged

        # the converged node fluxes supplied as a plain vector; the
        # auxiliary currents are initialized from the constitutive
        # relations
        outv = hbnlsolve(wp, (16,), sources, circuit, circuitdefs;
            ftol = 1e-12, x0 = Vector(out.nodeflux[:]))
        @test outv.solverinfo.converged
        @test isapprox(out.nodeflux, outv.nodeflux, atol = 1e-9)

        # the converged node fluxes supplied as a keyed array
        outk = hbnlsolve(wp, (16,), sources, circuit, circuitdefs;
            ftol = 1e-12, x0 = out.nodeflux)
        @test outk.solverinfo.converged
        @test isapprox(out.nodeflux, outk.nodeflux, atol = 1e-9)

        # an initial value with an invalid length throws
        @test_throws DimensionMismatch hbnlsolve(wp, (16,), sources,
            circuit, circuitdefs; ftol = 1e-12,
            x0 = zeros(Complex{Float64}, 3))
    end

    @testset "hbnlsolve mna gauge normalized dc warm start" begin

        # a common DC shift of every node flux of a floating component is
        # exactly the gauge degree of freedom: it changes no branch flux
        # and no physical quantity, so a shifted converged solution must be
        # an equally good initial guess. the guess is transformed into the
        # selected gauge before its residual is first evaluated.
        circuit, circuitdefs = jpacircuit()
        wp = (2*pi*4.75001*1e9,)
        sources = [(mode=(1,),port=1,current=0.00565e-6)]

        ref = hbnlsolve(wp, (8,), sources, circuit, circuitdefs;
            ftol = 1e-12, dc = true, odd = true)
        @test ref.solverinfo.converged
        Nmodes = length(ref.modes)
        dcindex = findfirst(==((0,)), collect(ref.modes))

        # shift the DC flux of node 1 (its own floating component) by a
        # large offset
        x0 = Vector(ref.nodeflux[:])
        x0[dcindex] += 5.0
        shifted = hbnlsolve(wp, (8,), sources, circuit, circuitdefs;
            ftol = 1e-12, dc = true, odd = true, x0 = x0)
        @test shifted.solverinfo.converged
        @test isapprox(Vector(ref.nodeflux[:]),
            Vector(shifted.nodeflux[:]), atol = 1e-9)
        # the gauge fixed flux is zero in the returned solution
        @test abs(shifted.nodeflux[dcindex]) <= 1e-12

    end

    @testset "hbnlsolve mna near balanced incompatible dc source" begin

        # an imbalance far below sqrt(eps) but far above summation
        # roundoff cannot correspond to a periodic excitation and must be
        # rejected rather than absorbed into the flux reference. the two
        # sources drive opposite DC currents into the same floating island
        # through two ports, mismatched by a relative 1e-10.
        circuit = Tuple{String,String,String,Union{Complex{Float64},Symbol,Int64}}[]
        push!(circuit,("P1","1","0",1))
        push!(circuit,("R1","1","0",:Rleft))
        push!(circuit,("Lj1","1","0",:Lj))
        push!(circuit,("C1","1","2",:Cc))
        push!(circuit,("P2","2","0",2))
        push!(circuit,("R2","2","0",:Rleft))
        push!(circuit,("L1","2","3",:L1))
        push!(circuit,("P3","3","0",3))
        push!(circuit,("R3","3","0",:Rleft))
        push!(circuit,("C3","3","0",:Cj))
        circuitdefs = Dict(:Lj =>1000.0e-12, :Cc => 100.0e-15,
            :Cj => 1000.0e-15, :Rleft => 50.0, :L1 => 500.0e-12)
        wp = (2*pi*4.75001*1e9,)
        Idc = 1.0e-6
        sources = [
            (mode=(1,),port=1,current=0.00565e-6),
            (mode=(0,),port=2,current=Idc),
            (mode=(0,),port=3,current=-Idc*(1 - 1e-10)),
        ]
        @test_throws ArgumentError hbnlsolve(wp, (8,), sources, circuit,
            circuitdefs; ftol = 1e-12, dc = true, odd = true)

        # exactly balanced sources through the same two ports are accepted
        sources2 = [
            (mode=(1,),port=1,current=0.00565e-6),
            (mode=(0,),port=2,current=Idc),
            (mode=(0,),port=3,current=-Idc),
        ]
        out = hbnlsolve(wp, (8,), sources2, circuit, circuitdefs;
            ftol = 1e-12, dc = true, odd = true)
        @test out.solverinfo.converged
    end

    @testset "hbnlsolve mna gauge only path" begin

        # a floating node whose only resistor is symbolic exercises the
        # gauge fixing equations with no promoted resistors (Naux == 0),
        # so gauge correctness is not entangled with resistor promotion.
        @variables w
        circuit = Tuple{String,String,String,Union{Complex{Float64},Symbol,Int64,Num}}[]
        push!(circuit,("P1","1","0",1))
        push!(circuit,("R1","1","0",50.0 + 1e-12*w))
        push!(circuit,("C1","1","2",:Cc))
        push!(circuit,("Lj1","2","0",:Lj))
        push!(circuit,("C2","2","0",:Cj))
        circuitdefs = Dict(:Lj =>1000.0e-12, :Cc => 100.0e-15,
            :Cj => 1000.0e-15)
        wp = (2*pi*4.75001*1e9,)
        sources = [(mode=(1,),port=1,current=0.00565e-6)]
        out = hbnlsolve(wp, (8,), sources, circuit, circuitdefs;
            ftol = 1e-12, dc = true, odd = true, symfreqvar = w)
        @test out.solverinfo.converged
        @test isapprox(abs(out.nodeflux(outputmode=(0,),node="1")), 0.0,
            atol = 1e-12)
    end

    @testset "hbnlsolve mna parallel interior resistors" begin

        # several parallel resistors on the same interior node pair are
        # physically identical to a single resistor of the parallel value.
        # interior resistors are not promoted, so these exercise the nodal
        # conductance stamps alongside the promoted port resistor.
        function circuitwith(rs)
            circuit = Tuple{String,String,String,Union{Complex{Float64},Symbol,Int64}}[]
            push!(circuit,("P1","1","0",1))
            push!(circuit,("R1","1","0",:Rleft))
            push!(circuit,("C1","1","2",:Cc))
            push!(circuit,("Lj1","2","0",:Lj))
            push!(circuit,("C2","2","0",:Cj))
            # the parallel bank sits on a non-port node pair
            for (k, r) in enumerate(rs)
                push!(circuit,("R$(k+1)","2","0",complex(r)))
            end
            return circuit
        end
        circuitdefs = Dict(:Lj =>1000.0e-12, :Cc => 100.0e-15,
            :Cj => 1000.0e-15, :Rleft => 50.0)
        wp = (2*pi*4.75001*1e9,)
        sources = [(mode=(1,),port=1,current=0.00565e-6)]
        many = hbnlsolve(wp, (16,), sources,
            circuitwith([30000.0, 30000.0, 30000.0]), circuitdefs;
            ftol = 1e-12)
        one_ = hbnlsolve(wp, (16,), sources, circuitwith([10000.0]),
            circuitdefs; ftol = 1e-12)
        @test many.solverinfo.converged
        @test one_.solverinfo.converged
        @test isapprox(Vector(many.nodeflux[:]), Vector(one_.nodeflux[:]),
            atol = 1e-9)

        # a symbolic resistor in parallel with a numeric resistor on an
        # internal node pair converges, with both in the conductance matrix
        # and only the port resistor promoted
        @variables w
        cmix = Tuple{String,String,String,Union{Complex{Float64},Symbol,Int64,Num}}[]
        push!(cmix,("P1","1","0",1))
        push!(cmix,("R1","1","0",complex(50.0)))
        push!(cmix,("C1","1","2",:Cc))
        push!(cmix,("Lj1","2","0",:Lj))
        push!(cmix,("C2","2","0",:Cj))
        push!(cmix,("R2","2","3",complex(1.0e6)))
        push!(cmix,("R3","2","3",2.0e6 + 1e-3*w))
        push!(cmix,("C3","3","0",:Cj))
        mix = hbnlsolve(wp, (8,), sources, cmix, circuitdefs;
            ftol = 1e-12, symfreqvar = w)
        @test mix.solverinfo.converged
    end

    @testset "mnagaugenormalize! and mnaungaugedkcl" begin
        # normalization subtracts the reference node DC flux from every
        # node of the component and leaves other modes untouched
        Nmodes = 2
        wmodes = [0.0, 2pi*4e9]
        x = Complex{Float64}[1+2im, 3, 4+2im, 5, 9, 9]
        JosephsonCircuits.mnagaugenormalize!(x, [[2,3]], wmodes, Nmodes)
        @test x == Complex{Float64}[0, 3, 3, 5, 9, 9]

        # the reconstruction removes exactly the gauge contribution x[g]
        # from the gauge rows of the node block and nothing else
        F = Complex{Float64}[0.5, 2.0, 0.0, 7.0, 9.0]
        xs = Complex{Float64}[0.5, 1.0, 0.25, 1.0, 1.0]
        Fkcl = JosephsonCircuits.mnaungaugedkcl(F, xs, [1, 3], 4)
        @test Fkcl == Complex{Float64}[0.0, 2.0, -0.25, 7.0]
        # a case a bare gauge-coordinate surrogate would wrongly accept:
        # the gauge flux is tiny but the augmented residual at the gauge
        # row is large, so the physical KCL residual F[g] - x[g] is large
        F2 = Complex{Float64}[0.3, 0.0]
        x2 = Complex{Float64}[1e-14, 0.0]
        @test abs(JosephsonCircuits.mnaungaugedkcl(F2, x2, [1], 2)[1]) > 0.29

        # the acceptance policy: block-relative infinity norms on both
        # sides, independent of the auxiliary state and of the number of
        # driven rows, and rejecting non-finite residuals
        b = Complex{Float64}[1.0, 1.0, 0.0]
        ftol = 1e-10
        # a compatible residual well inside the tolerance is accepted
        Fok = Complex{Float64}[1e-11, 0.0, 0.0]
        xok = Complex{Float64}[0.0, 0.0, 0.0, 1e9]
        ok, nk, tol = JosephsonCircuits.mnavalidatekcl(Fok, xok, [1], 3,
            b, ftol)
        @test ok
        @test nk == 1e-11
        # the applied tolerance is returned for diagnostics
        @test tol == 10*ftol*(1 + 1.0)
        # the surrogate-defeating case is rejected by the policy
        okbad, nkbad = JosephsonCircuits.mnavalidatekcl(
            Complex{Float64}[0.3, 0.0, 0.0], Complex{Float64}[1e-14, 0, 0, 0],
            [1], 3, b, ftol)
        @test !okbad
        @test nkbad > 0.29
        # duplicating driven rows does not loosen the per-row tolerance
        # (infinity-norm source scale), unlike a Euclidean scale
        bmany = Complex{Float64}[fill(1.0, 100); 0.0]
        Fedge = zeros(Complex{Float64}, 101)
        Fedge[1] = 25*ftol
        xz = zeros(Complex{Float64}, 102)
        okedge, = JosephsonCircuits.mnavalidatekcl(Fedge, xz, [1], 101,
            bmany, ftol)
        @test !okedge  # 25*ftol > 10*ftol*(1 + 1) regardless of row count
        # a non-finite reconstructed residual fails
        oknan, = JosephsonCircuits.mnavalidatekcl(
            Complex{Float64}[NaN, 0.0, 0.0], xok, [1], 3, b, ftol)
        @test !oknan
        # a non-finite source scale would make the tolerance infinite and
        # accept anything, so it fails the validation
        okinf, = JosephsonCircuits.mnavalidatekcl(Fok, xok, [1], 3,
            Complex{Float64}[Inf, 0.0, 0.0], ftol)
        @test !okinf
    end

    @testset "isnumericallyzero" begin
        # two pump tones one hertz apart near 5 GHz: the beat mode
        # frequency is a small number formed from ~3e10 scale terms, so
        # its rounding uncertainty is set by the pump frequencies. a 1 Hz
        # signal nominally cancelling the beat leaves a ~1e-6 remainder
        # which a collapsed-value bound misses by eight orders of
        # magnitude.
        w1 = 2*pi*5e9
        w2 = 2*pi*(5e9 - 1)
        wi = 2*pi*1.0
        total = wi + (w2 - w1)
        @test abs(total) > 1e-7                    # far above eps(2*pi)
        @test JosephsonCircuits.isnumericallyzero(total, (wi, -w1, w2))
        # a genuine 1 Hz signal with no pump contribution is not rejected
        @test !JosephsonCircuits.isnumericallyzero(wi, (wi, 0.0, 0.0))
        # a genuinely nonzero beat on the accepting side of the bound: a
        # 1 kHz signal against the 1 Hz beat gives a ~999 Hz total, far
        # above the roundoff of the ~6e10 contributing terms. the value
        # is the floating point combination of the terms, per the helper
        # contract.
        terms1k = (2*pi*1e3, -w1, w2)
        total1k = (terms1k[1] + terms1k[2]) + terms1k[3]
        @test abs(total1k) > 2*pi*900
        @test !JosephsonCircuits.isnumericallyzero(total1k, terms1k)
        # a nominally zero combination whose floating point evaluation
        # leaves a demonstrably nonzero remainder of several ulps of the
        # term scale is still classified as zero
        terms = (1e16, 1.0, -1e16, -1.0)
        acc = (((terms[1] + terms[2]) + terms[3]) + terms[4])
        @test !iszero(acc)
        @test JosephsonCircuits.isnumericallyzero(acc, terms)
        # non-finite inputs are rejected rather than misclassified
        @test_throws ArgumentError JosephsonCircuits.isnumericallyzero(
            Inf, (Inf, -w1))
        @test_throws ArgumentError JosephsonCircuits.isnumericallyzero(
            0.0, (NaN, -w1))
    end

    @testset "checkstaticstiffnessvalues and element contract" begin
        # zero and NaN inductances are rejected with informative errors
        @test_throws ArgumentError JosephsonCircuits.checkstaticstiffnessvalues(
            [:P,:R,:L], Any[1, 50.0, 0.0])
        @test_throws ArgumentError JosephsonCircuits.checkstaticstiffnessvalues(
            [:P,:R,:Lj], Any[1, 50.0, NaN])
        # an infinite inductance would make Lmean infinite; an open
        # circuit is represented by omitting the branch
        @test_throws ArgumentError JosephsonCircuits.checkstaticstiffnessvalues(
            [:P,:R,:L], Any[1, 50.0, Inf])
        # finite and symbolic values are accepted
        @variables Lsym
        @test isnothing(JosephsonCircuits.checkstaticstiffnessvalues(
            [:P,:R,:L,:L], Any[1, 50.0, 1e-9, Lsym]))

        # end-to-end: a zero inductance in a solve is rejected
        circuit = Tuple{String,String,String,Union{Complex{Float64},Symbol,Int64}}[]
        push!(circuit,("P1","1","0",1)); push!(circuit,("R1","1","0",:Rleft))
        push!(circuit,("C1","1","2",:Cc)); push!(circuit,("Lj1","2","0",:Lj))
        push!(circuit,("L1","2","0",:Lz)); push!(circuit,("C2","2","0",:Cj))
        circuitdefs = Dict(:Lj =>1000.0e-12, :Cc => 100.0e-15,
            :Cj => 1000.0e-15, :Rleft => 50.0, :Lz => 0.0)
        wp = (2*pi*4.75001*1e9,)
        sources = [(mode=(1,),port=1,current=0.00565e-6)]
        @test_throws ArgumentError hbnlsolve(wp, (8,), sources, circuit,
            circuitdefs; ftol = 1e-12)

        # end-to-end: an infinite inductance in a solve is rejected with
        # an informative error (an open circuit is represented by omitting
        # the branch, as in the floating island tests above)
        circuit2 = Tuple{String,String,String,Union{Complex{Float64},Symbol,Int64}}[]
        push!(circuit2,("P1","1","0",1)); push!(circuit2,("R1","1","0",:Rleft))
        push!(circuit2,("Lj1","1","0",:Lj))
        push!(circuit2,("C1","1","2",:Cc)); push!(circuit2,("L1","2","3",:L1))
        push!(circuit2,("L2","3","0",:Linf)); push!(circuit2,("C2","3","0",:Cj))
        circuitdefs2 = Dict(:Lj =>1000.0e-12, :Cc => 100.0e-15,
            :Cj => 1000.0e-15, :Rleft => 50.0, :L1 => 300.0e-12,
            :Linf => Inf)
        @test_throws ArgumentError hbnlsolve(wp, (8,), sources, circuit2,
            circuitdefs2; ftol = 1e-12, dc = true, odd = true)
    end

    @testset "non-finite frequency and source inputs" begin
        # non-finite drive frequencies and source currents are rejected
        # before mode construction and source assembly, in both solvers
        circuit, circuitdefs = jpacircuit()
        srcs(I) = [(mode=(1,),port=1,current=I)]
        for wbad in (Inf, NaN)
            @test_throws ArgumentError hbnlsolve((wbad,), (4,),
                srcs(1e-6), circuit, circuitdefs)
        end
        for Ibad in (Inf, NaN, complex(1e-6, Inf))
            @test_throws ArgumentError hbnlsolve((2*pi*5.0001e9,), (4,),
                srcs(Ibad), circuit, circuitdefs)
        end
        lincir = Tuple{String,String,String,Union{Complex{Float64},Symbol,Int64}}[]
        push!(lincir,("P1","1","0",1)); push!(lincir,("R1","1","0",:R))
        push!(lincir,("C1","1","2",:C)); push!(lincir,("L1","2","0",:L))
        push!(lincir,("P2","2","0",2)); push!(lincir,("R2","2","0",:R))
        lind = Dict(:R=>50.0,:C=>1e-13,:L=>1e-9)
        for wbad in (Inf, NaN)
            @test_throws ArgumentError hblinsolve([2*pi*5e9, wbad],
                lincir, lind)
        end
    end

    @testset "invalid inductances through hblinsolve" begin
        # the element value contract fires with the informative message in
        # the linearized solver as well: matrix assembly builds Inf or NaN
        # entries silently, so the contract check is empirically the first
        # failure for all of these values (verified for both solvers).
        cir = Tuple{String,String,String,Union{Complex{Float64},Symbol,Int64}}[]
        push!(cir,("P1","1","0",1)); push!(cir,("R1","1","0",:R))
        push!(cir,("C1","1","2",:C)); push!(cir,("Lj1","2","0",:Lj))
        push!(cir,("L1","2","0",:Lbad)); push!(cir,("C2","2","0",:C))
        for (v, fragment) in [(0.0, "A zero value"),
                              (Inf, "An infinite value"),
                              (NaN, "has a NaN value")]
            d = Dict(:R=>50.0,:C=>1e-13,:Lj=>1e-9,:Lbad=>v)
            err = try
                hblinsolve([2*pi*5.0001e9], cir, d)
                nothing
            catch e
                e
            end
            @test err isa ArgumentError
            @test occursin(fragment, sprint(showerror, err))
            errnl = try
                hbnlsolve((2*pi*5.0001e9,), (4,),
                    [(mode=(1,),port=1,current=1e-6)], cir, d)
                nothing
            catch e
                e
            end
            @test errnl isa ArgumentError
            @test occursin(fragment, sprint(showerror, errnl))
        end
    end

    @testset "debugJacobian with the MNA formulation" begin
        # the debug early return provides the augmented residual, Jacobian,
        # state, and assembly closure. with a DC gauge row present the
        # constant MNA structure is in place after one evaluation; no
        # components are promoted, so there are no auxiliary coordinates.
        circuit, circuitdefs = jpacircuit()
        wp = (2*pi*4.75001*1e9,)
        sources = [(mode=(1,),port=1,current=0.00565e-6)]
        d = hbnlsolve(wp, (8,), sources, circuit, circuitdefs;
            dc = true, odd = true, debugJacobian = true)
        F, J, x, fj! = d.F, d.Jx, d.x, d.fj
        Ntot = length(x)
        @test length(F) == Ntot
        @test size(J) == (Ntot, Ntot)
        # two non-ground nodes and no promoted components: the system is
        # purely nodal (port resistors stay as node conductances)
        @test Ntot % 2 == 0
        Nmodes = Ntot ÷ 2
        Nnodal = 2*Nmodes
        @test d.Nnodal == Nnodal
        @test isempty(d.mnaindices)
        @test length(d.gaugeindices) == 1
        # evaluate the closure at the zero state to fill the true residual
        # and Jacobian
        x .= 0
        fj!(F, J, x)
        @test all(isfinite, F)
        # the gauge fixing equation of the floating port node sits at the
        # DC coordinate of node 1: an exact one on the diagonal, since all
        # other contributions to that entry vanish at zero frequency
        @test J[1,1] == 1
        # the real Jacobian of the equivalent real system carries the same
        # augmented structure through the shared plan machinery
        Fr = copy(d.Fr)
        xr = copy(d.xr)
        xr .= 0
        d.fjreal(Fr, d.Jr, xr)
        @test all(isfinite, Fr)
    end

    @testset "hbnlsolve commensurate drive frequencies" begin
        # a nonzero mode tuple whose physical frequency cancels to zero
        # (or within roundoff of the combined magnitudes) would duplicate
        # the DC coordinate and produce vanishing capacitor and resistor
        # stamps. it must be rejected with an explanation.
        circuit, circuitdefs = jpacircuit()
        wp = (2*pi*5.0*1e9, 2*pi*10.0*1e9)
        sources = [
            (mode=(1,0),port=1,current=0.00565e-6),
            (mode=(0,1),port=1,current=0.00265e-6),
        ]
        # the mode (2,-1) has frequency 2*w1 - w2 = 0
        @test_throws ArgumentError hbnlsolve(wp, (2,2), sources, circuit,
            circuitdefs; ftol = 1e-12)
        # incommensurate drives are unaffected
        wp2 = (2*pi*5.0*1e9, 2*pi*10.00001*1e9)
        out = hbnlsolve(wp2, (2,2), sources, circuit, circuitdefs;
            ftol = 1e-12)
        @test out.solverinfo.converged
    end

    @testset "hblinsolve beat frequency cancellation" begin
        # two pump tones one hertz apart: the retained four wave mixing
        # signal mode (1,-1) has the -1 Hz beat frequency, formed from
        # ~5 GHz terms. a +1 Hz signal nominally cancels it, leaving a
        # remainder of order 1e-6 rad/s - far above eps of the beat but
        # far below the rounding scale of the contributing pump terms -
        # which must be rejected by the term-aware bound.
        circuit, circuitdefs = jpacircuit()
        wp = (2*pi*(5.0e9 - 1), 2*pi*5.0e9)
        sources = [
            (mode=(1,0),port=1,current=0.00565e-6),
            (mode=(0,1),port=1,current=0.00265e-6),
        ]
        @test_throws ArgumentError hbsolve([2*pi*1.0], wp, sources, (1,1),
            (2,2), circuit, circuitdefs; ftol = 1e-12,
            threewavemixing = false, fourwavemixing = true)
        # a signal well away from the beat is accepted
        out = hbsolve([2*pi*1.0e3], wp, sources, (1,1), (2,2), circuit,
            circuitdefs; ftol = 1e-12, threewavemixing = false,
            fourwavemixing = true)
        @test out.nonlinear.solverinfo.converged
    end

    @testset "perfectly coupled inductors" begin
        # a perfectly coupled inductor pair (|k| = 1) makes the branch
        # inductance matrix singular, so its inverse does not exist. the
        # modified nodal analysis promotion keeps the branch inductance
        # matrix un-inverted, so the system remains well posed whenever
        # the surrounding circuit determines the branch currents.
        cir = Tuple{String,String,String,Union{Complex{Float64},Symbol,Int64}}[]
        push!(cir,("P1","1","0",1)); push!(cir,("R1","1","0",:R))
        push!(cir,("L1","1","0",:L)); push!(cir,("L2","2","0",:L))
        push!(cir,("K1","L1","L2",:K)); push!(cir,("C1","2","0",:C))
        push!(cir,("Lj1","2","0",:Lj))
        d = Dict(:R=>50.0,:L=>1e-9,:K=>1.0,:C=>1e-12,:Lj=>1e-9)
        # this previously threw a SingularException from the inverse of the
        # singular branch inductance matrix in numericmatrices. with the
        # coupled branches excluded from the inverse inductance matrix and
        # represented by auxiliary branch currents with the un-inverted
        # branch inductance matrix, the perfectly coupled pair is well
        # posed and solves.
        outk = hbnlsolve((2*pi*5.0001e9,),(4,),
            [(mode=(1,),port=1,current=1e-6)],cir,d)
        @test outk.solverinfo.converged
    end

    @testset "hbnlsolve large auxiliary scale acceptance" begin
        # in an inductor-free circuit the Lmean = 1 fallback makes the
        # auxiliary currents u = Lmean*i/phi0 of order 1e9 while the node
        # fluxes are of order 1. the physical KCL acceptance must not be
        # loosened by these large auxiliary entries: the solve converges,
        # the gauge fixed DC flux is exactly zero, and the pump mode
        # matches the analytic value.
        circuit = [("P1","1","0",1),("R1","1","0",50.0)]
        wpz = 2*pi*5.0*1e9
        Ip = 1.0e-6
        out = hbnlsolve((wpz,),(1,),[(mode=(1,),port=1,current=Ip)],
            circuit, Dict(); dc = true, odd = true)
        @test out.solverinfo.converged
        @test isapprox(abs(out.nodeflux[1]), 0.0, atol = 1e-15)
        @test isapprox(im*out.nodeflux[2]*wpz*JosephsonCircuits.phi0/50, Ip)
    end

    @testset "hbnlsolve mna keyedarrays false" begin

        circuit, circuitdefs = jpacircuit()
        wp = (2*pi*4.75001*1e9,)
        sources = [(mode=(1,),port=1,current=0.00565e-6)]

        keyed = hbnlsolve(wp, (8,), sources, circuit, circuitdefs;
            ftol = 1e-12)
        plain = hbnlsolve(wp, (8,), sources, circuit, circuitdefs;
            ftol = 1e-12, keyedarrays = false)

        # the output length is the number of node fluxes, without the
        # auxiliary variables, and the values match the keyed output
        @test length(plain.nodeflux) == 2*length(keyed.modes)
        @test length(plain.nodeflux) == length(keyed.nodeflux)
        @test isapprox(Vector(keyed.nodeflux[:]), plain.nodeflux,
            atol = 1e-12)
    end

    @testset "hbsolve nodal reference values" begin

        # full pipeline reference values computed with the purely nodal
        # formulation prior to its removal; the linearized gain amplifies
        # small operating point differences, hence the looser tolerance.
        circuit, circuitdefs = jpacircuit()
        ws = 2*pi*(4.6:0.1:4.9)*1e9
        wp = (2*pi*4.75001*1e9,)
        sources = [(mode=(1,),port=1,current=0.00565e-6)]

        out = hbsolve(ws, wp, sources, (4,), (8,), circuit, circuitdefs;
            ftol = 1e-12)
        @test out.nonlinear.solverinfo.converged
        @test isapprox(Vector(out.linearized.S(outputmode=(0,),
                outputport=1, inputmode=(0,), inputport=1)),
            ComplexF64[0.7946384428948396 - 0.6107482382057416im,
                0.38523539452862193 - 1.0146026970089426im,
                0.8028844946933809 + 0.7303408200854153im,
                0.9926807608401819 + 0.13802125272539348im], atol = 1e-6)
        @test isapprox(Vector(out.linearized.QE(outputmode=(0,),
                outputport=1, inputmode=(0,), inputport=1)),
            [0.9955631849311455, 0.8687783963970223, 0.8686637478305292,
                0.9955599960987674], atol = 1e-6)
    end

    @testset "calcstaticfluxcomponents and calcdcgaugeindices" begin
        # node 2 is isolated (no static edges); nodes 3 and 4 are joined by
        # a junction with no path to ground
        comps = JosephsonCircuits.calcstaticfluxcomponents(
            [:P,:R,:C,:Lj,:C],[2 2 2 3 4;1 1 3 4 1],
            Any[1,50.0,1e-13,1e-9,1e-12],4)
        @test comps == [[2],[3,4]]
        # an inductor to ground removes the gauge freedom
        @test JosephsonCircuits.calcstaticfluxcomponents(
            [:P,:R,:L],[2 2 2;1 1 1],Any[1,50.0,1e-9],2) == Vector{Int64}[]
        # defense in depth on inputs the public solvers reject before
        # classification (checkstaticstiffnessvalues): a zero value, were
        # it to reach the classifier, would be treated as a finite
        # stiffness edge, and a NaN, like an infinity, contributes no edge
        @test JosephsonCircuits.calcstaticfluxcomponents(
            [:P,:R,:L],[2 2 2;1 1 1],Any[1,50.0,0.0],2) == Vector{Int64}[]
        @test JosephsonCircuits.calcstaticfluxcomponents(
            [:P,:R,:L],[2 2 2;1 1 1],Any[1,50.0,NaN],2) == [[2]]
        # an infinite inductance is an open circuit at zero frequency and
        # provides no static stiffness, so the node remains floating
        @test JosephsonCircuits.calcstaticfluxcomponents(
            [:P,:R,:L],[2 2 2;1 1 1],Any[1,50.0,Inf],2) == [[2]]
        # symbolic values are conservatively assumed to provide stiffness
        @variables Lsym
        @test JosephsonCircuits.calcstaticfluxcomponents(
            [:P,:R,:L],[2 2 2;1 1 1],Any[1,50.0,Lsym],2) == Vector{Int64}[]
        # one gauge per floating component and zero-frequency mode
        @test JosephsonCircuits.calcdcgaugeindices(
            comps,[0.0,2pi*4e9],2) == [1,3]
        @test JosephsonCircuits.calcdcgaugeindices(
            comps,[2pi*4e9,2pi*8e9],2) == Int[]
    end

    @testset "calcAmna stamp identity" begin

        # eliminating the auxiliary variables from the mna stamp must
        # reproduce the nodal conductance contribution im*Gnm*wmodesm
        # exactly, for positive, negative, and zero mode frequencies.
        componenttypes = [:R]
        nodeindices = reshape([2,3],2,1)
        vvn = [50.0]
        Nnodes = 4
        Nmodes = 3
        w0 = 2pi*4e9
        wmodes = [-w0, 0.0, w0]
        Lmean = 1e-9
        Nnodal = (Nnodes-1)*Nmodes

        # promote the resistor by explicit index: calcAmna stamps whatever
        # index set it is given, independent of the promotion policy
        mnaindices = [1]
        Amna = JosephsonCircuits.calcAmna(mnaindices, nodeindices, vvn,
            Int[], wmodes, Nmodes, Nnodes, Lmean)

        Gnm = JosephsonCircuits.calcGn(componenttypes, nodeindices, vvn,
            Nmodes, Nnodes)
        JosephsonCircuits.conjnegfreq!(Gnm, wmodes)
        wmodesm = JosephsonCircuits.LinearAlgebra.Diagonal(
            repeat(wmodes, outer = Nnodes-1))

        # a deterministic test vector keeps this exact algebraic identity
        # independent of the global random number generator state
        phi = Complex{Float64}[complex(sin(k), cos(2k)) for k in 1:Nnodal]
        x = vcat(phi, zeros(Complex{Float64}, Nmodes))
        JosephsonCircuits.mnainitialaux!(x, Amna, Nnodal)
        r = Amna*x

        # the constitutive rows are exactly zero after consistent
        # initialization of the auxiliary currents
        @test all(iszero, r[Nnodal+1:end])
        # the Kirchhoff current law rows reproduce the nodal stamp
        @test isapprox(r[1:Nnodal], im*Lmean*Gnm*wmodesm*phi,
            rtol = 1e-14, atol = 1e-30)
    end

end

@testset verbose=true "mna hblinsolve" begin

    # analytic two-port network parameters from a chain of ABCD matrices,
    # for validating the linearized solver against exact circuit theory.
    # elements: (:series, Z) or (:shunt, Y).
    function abcdchain(elements, w)
        A = ComplexF64[1 0; 0 1]
        for (kind, val) in elements
            if kind == :series
                A = A*ComplexF64[1 val(w); 0 1]
            else
                A = A*ComplexF64[1 0; val(w) 1]
            end
        end
        return A
    end
    function abcdtoS(A, Z0)
        a, b, c, d = A[1,1], A[1,2], A[2,1], A[2,2]
        den = a + b/Z0 + c*Z0 + d
        S11 = (a + b/Z0 - c*Z0 - d)/den
        S21 = 2/den
        S12 = 2*(a*d - b*c)/den
        S22 = (-a + b/Z0 - c*Z0 + d)/den
        return S11, S12, S21, S22
    end
    @testset "hblinsolve analytic linear network" begin
        circuit = Tuple{String,String,String,Union{Complex{Float64},Symbol,Int64}}[]
        push!(circuit,("P1","1","0",1)); push!(circuit,("R1","1","0",:Rleft))
        push!(circuit,("C1","1","2",:Cc)); push!(circuit,("L2","2","0",:L2))
        push!(circuit,("C2","2","0",:C2)); push!(circuit,("C3","2","3",:Cc))
        push!(circuit,("P2","3","0",2)); push!(circuit,("R2","3","0",:Rleft))
        circuitdefs = Dict(:Rleft => 50.0, :Cc => 30.0e-15, :L2 => 1.2e-9,
            :C2 => 0.9e-12)
        ws = collect(2*pi*(3.0:0.1:7.0)*1e9)
        out = hblinsolve(ws, circuit, circuitdefs)
        elements = [
            (:series, w -> 1/(im*w*circuitdefs[:Cc])),
            (:shunt, w -> 1/(im*w*circuitdefs[:L2]) + im*w*circuitdefs[:C2]),
            (:series, w -> 1/(im*w*circuitdefs[:Cc])),
        ]
        for (i, w) in enumerate(ws)
            S11, S12, S21, S22 = abcdtoS(abcdchain(elements, w), 50.0)
            @test isapprox(out.S(outputmode=(0,),outputport=1,
                inputmode=(0,),inputport=1,freqindex=i), S11, atol = 1e-10)
            @test isapprox(out.S(outputmode=(0,),outputport=2,
                inputmode=(0,),inputport=1,freqindex=i), S21, atol = 1e-10)
        end
    end

    @testset "hblinsolve zero frequency error" begin
        circuit = Tuple{String,String,String,Union{Complex{Float64},Symbol,Int64}}[]
        push!(circuit,("P1","1","0",1)); push!(circuit,("R1","1","0",:Rleft))
        push!(circuit,("R2","1","2",:Rmid))
        push!(circuit,("P2","2","0",2)); push!(circuit,("R3","2","0",:Rleft))
        circuitdefs = Dict(:Rleft => 50.0, :Rmid => 100.0)
        @test_throws ArgumentError hblinsolve([0.0], circuit, circuitdefs)
        # this network is frequency independent, so the result at 1 Hz
        # equals its DC limit exactly
        S = hblinsolve([2*pi*1.0], circuit, circuitdefs)
        @test isapprox(S.S(outputmode=(0,),outputport=2,inputmode=(0,),
            inputport=1,freqindex=1), 0.5, atol = 1e-10)
    end

    @testset "hblinsolve mna internal resistor" begin
        # an internal damping resistor between two non-port nodes is
        # promoted along with the port resistors, exercising the auxiliary
        # current Kirchhoff couplings away from the ports.
        circuit = Tuple{String,String,String,Union{Complex{Float64},Symbol,Int64}}[]
        push!(circuit,("P1","1","0",1)); push!(circuit,("R1","1","0",:Rleft))
        push!(circuit,("C1","1","2",:Cc)); push!(circuit,("L2","2","0",:L2))
        push!(circuit,("C2","2","0",:C2)); push!(circuit,("R4","2","3",:Rint))
        push!(circuit,("C3","3","0",:C2)); push!(circuit,("C4","3","4",:Cc))
        push!(circuit,("P2","4","0",2)); push!(circuit,("R2","4","0",:Rleft))
        circuitdefs = Dict(:Rleft => 50.0, :Cc => 30.0e-15, :L2 => 1.2e-9,
            :C2 => 0.9e-12, :Rint => 700.0)
        ws = collect(2*pi*(3.0:0.2:7.0)*1e9)
        out = hblinsolve(ws, circuit, circuitdefs; returnSnoise = true)
        elements = [
            (:series, w -> 1/(im*w*circuitdefs[:Cc])),
            (:shunt, w -> 1/(im*w*circuitdefs[:L2]) + im*w*circuitdefs[:C2]),
            (:series, w -> complex(circuitdefs[:Rint])),
            (:shunt, w -> im*w*circuitdefs[:C2]),
            (:series, w -> 1/(im*w*circuitdefs[:Cc])),
        ]
        for (i, w) in enumerate(ws)
            A = abcdchain(elements, w)
            S11, S12, S21, S22 = abcdtoS(A, 50.0)
            @test isapprox(out.S(outputmode=(0,),outputport=1,
                inputmode=(0,),inputport=1,freqindex=i), S11, atol = 1e-10)
            @test isapprox(out.S(outputmode=(0,),outputport=2,
                inputmode=(0,),inputport=1,freqindex=i), S21, atol = 1e-10)
        end
        # the internal resistor contributes the only noise component; the
        # port resistors are the scattering parameter references. the
        # keyed axes are (inputmode, component, outputmode, outputport,
        # freqindex).
        @test size(out.Snoise, 2) == 1
        @test all(isfinite, out.Snoise)
    end

    @testset "hbsolve mna noise and output slicing" begin
        # exercise the adjoint solve (noise, quantum efficiency, and
        # commutation relations) and the node flux and voltage output
        # slicing with the taller augmented solution matrix.
        circuit = Tuple{String,String,String,Union{Complex{Float64},Symbol,Int64}}[]
        push!(circuit,("P1","1","0",1)); push!(circuit,("R1","1","0",:Rleft))
        push!(circuit,("C1","1","2",:Cc)); push!(circuit,("Lj1","2","0",:Lj))
        push!(circuit,("C2","2","0",:Cj)); push!(circuit,("R3","2","0",:Rint))
        circuitdefs = Dict(:Lj =>1000.0e-12, :Cc => 100.0e-15,
            :Cj => 1000.0e-15, :Rleft => 50.0, :Rint => 20000.0)
        ws = 2*pi*(4.6:0.1:4.9)*1e9
        wp = (2*pi*4.75001*1e9,)
        sources = [(mode=(1,),port=1,current=0.00565e-6)]
        kwargs = (ftol = 1e-12, returnSnoise = true, returnnodeflux = true,
            returnvoltage = true, returnnodefluxadjoint = true,
            keyedarrays = false)
        out = hbsolve(ws, wp, sources, (4,), (8,), circuit, circuitdefs;
            kwargs...)
        @test out.nonlinear.solverinfo.converged
        # quantum efficiency is a physical efficiency
        @test all(q -> 0 < q <= 1 + 1e-10, out.linearized.QE)
        @test all(isfinite, out.linearized.S)
        @test all(isfinite, out.linearized.Snoise)
        # the returned node quantities contain only the node rows: two
        # resistors are promoted internally, but the outputs are nodal
        Nsignalmodes = size(out.linearized.S, 1)  # one port
        Nnodal = 2*Nsignalmodes
        @test size(out.linearized.nodeflux, 1) == Nnodal
        @test size(out.linearized.voltage, 1) == Nnodal
        @test size(out.linearized.nodefluxadjoint, 1) == Nnodal
        # the voltage output is exactly im*w*phi row by row, which
        # validates the slicing of both arrays from the augmented solution
        wpumpmodes = [m[1]*wp[1] for m in out.linearized.modes]
        for i in eachindex(ws), col in axes(out.linearized.nodeflux, 2)
            for n in 0:1, (k, wm) in enumerate(wpumpmodes)
                row = n*Nsignalmodes + k
                @test isapprox(out.linearized.voltage[row, col, i],
                    im*(ws[i] + wm)*out.linearized.nodeflux[row, col, i],
                    rtol = 1e-12, atol = 1e-30)
            end
        end
    end

    @testset "calcAmnasplit stamp identity" begin
        # eliminating the auxiliary currents from the per-frequency
        # linearized assembly Amna0 + im*AmnaG*wmodesm (with the same
        # negative-mode conjugation applied as for the conductance matrix)
        # must reproduce the removed nodal stamp exactly. two promoted
        # resistors share one node pair and one is stored as a zero
        # imaginary part ComplexF64.
        componenttypes = [:R, :R]
        nodeindices = [2 2; 3 3]
        vvn = Any[100.0, 300.0 + 0.0im]
        Nnodes = 3
        Nmodes = 3
        w0 = 2pi*5e9
        wmodes = [-w0 + 2pi*4e9, 2pi*1e5, w0 + 2pi*4e9]
        Nnodal = (Nnodes-1)*Nmodes
        # promote both resistors by explicit index, independent of the
        # promotion policy
        mnaindices = [1, 2]
        Amna0, AmnaG = JosephsonCircuits.calcAmnasplit(mnaindices,
            nodeindices, vvn, Nmodes, Nnodes)
        Ntot = size(Amna0, 1)
        Naux = Ntot - Nnodal
        wmodesmfull = JosephsonCircuits.LinearAlgebra.Diagonal(
            repeat(wmodes, outer = Nnodes-1+length(mnaindices)))
        conjmask = real.(wmodesmfull) .< 0

        # assemble the augmented matrix with the production machinery,
        # zeroing the values while keeping the sparsity structure
        X = JosephsonCircuits.spaddkeepzeros(Amna0, AmnaG)
        fill!(X.nzval, 0)
        m0 = JosephsonCircuits.sparseaddmap(X, Amna0)
        mG = JosephsonCircuits.sparseaddmap(X, AmnaG)
        JosephsonCircuits.sparseadd!(X, 1, Amna0, m0)
        JosephsonCircuits.sparseaddconjsubst!(X, im, AmnaG, wmodesmfull,
            mG, conjmask, wmodesmfull, JosephsonCircuits.symbolicindices(AmnaG),
            nothing)

        # assemble the nodal stamp of the same resistors identically
        Gp = JosephsonCircuits.calcGn(componenttypes, nodeindices, vvn,
            Nmodes, Nnodes)
        Gp = JosephsonCircuits.mnapad(Gp, Naux)
        Y = copy(Gp)
        fill!(Y.nzval, 0)
        mY = JosephsonCircuits.sparseaddmap(Y, Gp)
        JosephsonCircuits.sparseaddconjsubst!(Y, im, Gp, wmodesmfull,
            mY, conjmask, wmodesmfull, JosephsonCircuits.symbolicindices(Gp),
            nothing)

        # Schur elimination of the auxiliary block of X onto the nodes
        Xd = Matrix(X)
        A11 = Xd[1:Nnodal, 1:Nnodal]
        A12 = Xd[1:Nnodal, Nnodal+1:end]
        A21 = Xd[Nnodal+1:end, 1:Nnodal]
        A22 = Xd[Nnodal+1:end, Nnodal+1:end]
        schur = A11 - A12*(A22\A21)
        @test isapprox(schur, Matrix(Y)[1:Nnodal, 1:Nnodal],
            rtol = 1e-14, atol = 1e-30)
    end

    @testset "hblinsolve zero total sideband cancellation" begin
        # a signal frequency intended to cancel a pump mode frequency can
        # leave an ulp-level remainder when the two are computed through
        # separate arithmetic paths. this must be detected, not passed
        # through to a catastrophically ill conditioned flux-basis matrix.
        circuit = Tuple{String,String,String,Union{Complex{Float64},Symbol,Int64}}[]
        push!(circuit,("P1","1","0",1)); push!(circuit,("R1","1","0",:Rleft))
        push!(circuit,("C1","1","2",:Cc)); push!(circuit,("Lj1","2","0",:Lj))
        push!(circuit,("C2","2","0",:Cj))
        circuitdefs = Dict(:Lj =>1000.0e-12, :Cc => 100.0e-15,
            :Cj => 1000.0e-15, :Rleft => 50.0)
        wp = (2*pi*4.75001*1e9,)
        sources = [(mode=(1,),port=1,current=0.00565e-6)]
        # three wave mixing includes the (-1,) signal mode, so a signal at
        # (numerically) the pump frequency produces a zero total frequency
        ws = [wp[1]*(1 + 2*eps())]
        @test_throws ArgumentError hbsolve(ws, wp, sources, (4,), (8,),
            circuit, circuitdefs; ftol = 1e-12, threewavemixing = true)
    end

    @testset "hblinsolve complex storage and symbolic resistor" begin
        # the promotion predicate paths already unit tested for the
        # nonlinear solver are exercised through the linearized assembly:
        # real resistances stored as ComplexF64 are promoted and match the
        # analytic network exactly, and a symbolic frequency dependent
        # resistor is retained in the conductance matrix.
        circuit = Tuple{String,String,String,Union{Complex{Float64},Symbol,Int64}}[]
        push!(circuit,("P1","1","0",1)); push!(circuit,("R1","1","0",:Rleft))
        push!(circuit,("C1","1","2",:Cc)); push!(circuit,("L2","2","0",:L2))
        push!(circuit,("C2","2","0",:C2)); push!(circuit,("C3","2","3",:Cc))
        push!(circuit,("P2","3","0",2)); push!(circuit,("R2","3","0",:Rleft))
        circuitdefs = Dict{Symbol,Complex{Float64}}(:Rleft => 50.0,
            :Cc => 30.0e-15, :L2 => 1.2e-9, :C2 => 0.9e-12)
        ws = collect(2*pi*(3.0:0.5:7.0)*1e9)
        out = hblinsolve(ws, circuit, circuitdefs)
        elements = [
            (:series, w -> 1/(im*w*real(circuitdefs[:Cc]))),
            (:shunt, w -> 1/(im*w*real(circuitdefs[:L2])) +
                im*w*real(circuitdefs[:C2])),
            (:series, w -> 1/(im*w*real(circuitdefs[:Cc]))),
        ]
        for (i, w) in enumerate(ws)
            S11, S12, S21, S22 = abcdtoS(abcdchain(elements, w), 50.0)
            @test isapprox(out.S(outputmode=(0,),outputport=2,
                inputmode=(0,),inputport=1,freqindex=i), S21, atol = 1e-10)
        end

        # a symbolic frequency dependent resistor across the resonator
        @variables w
        circuit2 = Tuple{String,String,String,Union{Complex{Float64},Symbol,Int64,Num}}[]
        push!(circuit2,("P1","1","0",1)); push!(circuit2,("R1","1","0",complex(50.0)))
        push!(circuit2,("C1","1","2",:Cc)); push!(circuit2,("L2","2","0",:L2))
        push!(circuit2,("C2","2","0",:C2))
        push!(circuit2,("R3","2","0",1.0e4 + 1e-6*w))
        push!(circuit2,("C3","2","3",:Cc))
        push!(circuit2,("P2","3","0",2)); push!(circuit2,("R2","3","0",complex(50.0)))
        circuitdefs2 = Dict(:Cc => 30.0e-15, :L2 => 1.2e-9, :C2 => 0.9e-12)
        out2 = hblinsolve(ws, circuit2, circuitdefs2; symfreqvar = w)
        elements2 = [
            (:series, wv -> 1/(im*wv*30.0e-15)),
            (:shunt, wv -> 1/(im*wv*1.2e-9) + im*wv*0.9e-12 +
                1/(1.0e4 + 1e-6*wv)),
            (:series, wv -> 1/(im*wv*30.0e-15)),
        ]
        for (i, wv) in enumerate(ws)
            S11, S12, S21, S22 = abcdtoS(abcdchain(elements2, wv), 50.0)
            @test isapprox(out2.S(outputmode=(0,),outputport=2,
                inputmode=(0,),inputport=1,freqindex=i), S21, atol = 1e-10)
        end
    end

    @testset "hblinsolve batches and output forms" begin
        # threaded batch splitting must not change results, keyed and
        # plain outputs must agree, and the adjoint output flags must
        # produce finite node-sized arrays with a promoted internal
        # resistor present. the sensitivity flags are exercised separately
        # below.
        circuit = Tuple{String,String,String,Union{Complex{Float64},Symbol,Int64}}[]
        push!(circuit,("P1","1","0",1)); push!(circuit,("R1","1","0",:Rleft))
        push!(circuit,("C1","1","2",:Cc)); push!(circuit,("L2","2","0",:L2))
        push!(circuit,("C2","2","0",:C2)); push!(circuit,("R4","2","3",:Rint))
        push!(circuit,("C3","3","0",:C2)); push!(circuit,("C4","3","4",:Cc))
        push!(circuit,("P2","4","0",2)); push!(circuit,("R2","4","0",:Rleft))
        circuitdefs = Dict(:Rleft => 50.0, :Cc => 30.0e-15, :L2 => 1.2e-9,
            :C2 => 0.9e-12, :Rint => 700.0)
        ws = collect(2*pi*(3.0:0.25:7.0)*1e9)
        b1 = hblinsolve(ws, circuit, circuitdefs; nbatches = 1)
        b3 = hblinsolve(ws, circuit, circuitdefs; nbatches = 3)
        @test isapprox(b1.S, b3.S, atol = 1e-14)

        keyed = hblinsolve(ws, circuit, circuitdefs;
            returnnodeflux = true, returnvoltage = true)
        plain = hblinsolve(ws, circuit, circuitdefs;
            returnnodeflux = true, returnvoltage = true,
            keyedarrays = false)
        # the keyed outputs split the mode/node and mode/port axes; the
        # column-major reshape recovers the plain layout
        @test isapprox(reshape(Array(keyed.nodeflux),
            size(plain.nodeflux)), plain.nodeflux, atol = 1e-14)
        @test isapprox(reshape(Array(keyed.voltage),
            size(plain.voltage)), plain.voltage, atol = 1e-14)
        # node-sized outputs: 4 non-ground nodes, one signal mode
        @test size(plain.nodeflux, 1) == 4
        @test size(plain.voltage, 1) == 4

        adj = hblinsolve(ws, circuit, circuitdefs; keyedarrays = false,
            returnnodefluxadjoint = true, returnvoltageadjoint = true)
        @test all(isfinite, adj.nodefluxadjoint)
        @test all(isfinite, adj.voltageadjoint)
        @test size(adj.nodefluxadjoint, 1) == 4
        @test size(adj.voltageadjoint, 1) == 4

        # the sensitivity output flags, with the promoted internal
        # resistor itself as the sensitivity component: finite outputs of
        # the expected leading dimension. these are structural checks of
        # the output paths through the augmented solution matrix; the
        # sensitivity semantics themselves are upstream functionality.
        sens = hblinsolve(ws, circuit, circuitdefs; keyedarrays = false,
            sensitivitynames = ["R4"], returnSsensitivity = true)
        @test all(isfinite, sens.Ssensitivity)
        @test size(sens.Ssensitivity)[1:2] == size(sens.S)[1:2]
        @test size(sens.Ssensitivity, 3) == 1
    end

    @testset "hblinsolve noise port symmetry" begin
        # an internal resistor at the symmetry center of a mirror
        # symmetric two-port network must couple with equal magnitude to
        # both ports. this checks the noise wave transfer through the
        # augmented solve without depending on a normalization convention.
        circuit = Tuple{String,String,String,Union{Complex{Float64},Symbol,Int64}}[]
        push!(circuit,("P1","1","0",1)); push!(circuit,("R1","1","0",:Rleft))
        push!(circuit,("C1","1","2",:Cc)); push!(circuit,("R3","2","0",:Rint))
        push!(circuit,("L2","2","0",:L2))
        push!(circuit,("C3","2","3",:Cc))
        push!(circuit,("P2","3","0",2)); push!(circuit,("R2","3","0",:Rleft))
        circuitdefs = Dict(:Rleft => 50.0, :Cc => 30.0e-15, :L2 => 1.2e-9,
            :Rint => 5000.0)
        ws = collect(2*pi*(3.0:0.5:7.0)*1e9)
        out = hblinsolve(ws, circuit, circuitdefs; returnSnoise = true,
            keyedarrays = false)
        # one noise component (the internal resistor; port resistors are
        # the scattering references), two output ports
        @test size(out.Snoise, 1) == 1
        @test size(out.Snoise, 2) == 2
        for i in eachindex(ws)
            @test isapprox(abs(out.Snoise[1,1,i]), abs(out.Snoise[1,2,i]),
                rtol = 1e-10)
        end
    end

    @testset "hbsolve mna dc singular end to end" begin
        # the full pipeline on a circuit with no inductive path to ground
        # from the port node: the pump solve with dc = true requires the
        # gauge equations, and the linearized solve then runs on the MNA
        # operating point. the nodal formulation cannot solve this circuit
        # with a DC pump mode at all.
        circuit = Tuple{String,String,String,Union{Complex{Float64},Symbol,Int64}}[]
        push!(circuit,("P1","1","0",1)); push!(circuit,("R1","1","0",:Rleft))
        push!(circuit,("C1","1","2",:Cc)); push!(circuit,("Lj1","2","0",:Lj))
        push!(circuit,("C2","2","0",:Cj))
        circuitdefs = Dict(:Lj =>1000.0e-12, :Cc => 100.0e-15,
            :Cj => 1000.0e-15, :Rleft => 50.0)
        ws = 2*pi*(4.6:0.1:4.9)*1e9
        wp = (2*pi*4.75001*1e9,)
        sources = [(mode=(1,),port=1,current=0.00565e-6)]
        withdc = hbsolve(ws, wp, sources, (4,), (8,), circuit,
            circuitdefs; ftol = 1e-12, dc = true, threewavemixing = true)
        @test withdc.nonlinear.solverinfo.converged
        @test all(isfinite, withdc.linearized.S)
        # with no DC drive, including the DC pump mode must not change the
        # scattering parameters
        plain = hbsolve(ws, wp, sources, (4,), (8,), circuit, circuitdefs;
            ftol = 1e-12)
        @test isapprox(
            plain.linearized.S(outputmode=(0,),outputport=1,
                inputmode=(0,),inputport=1),
            withdc.linearized.S(outputmode=(0,),outputport=1,
                inputmode=(0,),inputport=1), atol = 1e-6)
    end

    @testset "solver scale" begin

        # Z0/w0 with geometric means over port impedances and nonzero
        # drive frequencies
        @test isapprox(JosephsonCircuits.calcsolverscale((2*pi*5e9,),
            [:P,:R], Any[1,50.0], [2], 1e-9), 50.0/(2*pi*5e9))
        # zero frequencies are ignored; geomean(4,9) = 6, geomean(25,100) = 50
        @test isapprox(JosephsonCircuits.calcsolverscale(
            (2*pi*4e9, 0.0, 2*pi*9e9), [:R,:R], Any[25.0,100.0], [1,2],
            1.0), 50.0/(2*pi*6e9))
        # fallback to all constant real resistors when no port impedances
        @test isapprox(JosephsonCircuits.calcsolverscale((2*pi*1e9,),
            [:R], Any[200.0], Int[], 1.0), 200.0/(2*pi*1e9))
        # fallback to 50 ohms when no resistors at all
        @test isapprox(JosephsonCircuits.calcsolverscale((2*pi*1e9,),
            [:C], Any[1e-15], Int[], 1.0), 50.0/(2*pi*1e9))
        # purely zero-frequency drive falls back to the mean inductance
        @test JosephsonCircuits.calcsolverscale((0.0,), [:R], Any[50.0],
            [1], 3.3e-9) == 3.3e-9
        # complex storage with zero imaginary part is accepted
        @test isapprox(JosephsonCircuits.calcsolverscale((2*pi*5e9,),
            [:R], Any[50.0+0.0im], [1], 1e-9), 50.0/(2*pi*5e9))
    end

    @testset "coupled inductor MNA stamp Schur identity" begin

        # eliminating the auxiliary branch currents from the coupled
        # inductor stamps must reproduce exactly the coupled part of the
        # scaled nodal inverse inductance matrix, Lscale*Rbn'*inv(L)*Rbn,
        # assembled by the production conjugation-free machinery.
        Lb = JosephsonCircuits.SparseArrays.sparsevec([1,2],
            [3e-10, 5e-10], 3)
        k = 0.9
        M = k*sqrt(3e-10*5e-10)
        Mb = JosephsonCircuits.SparseArrays.sparse([1,2],[2,1],[M,M],3,3)
        Rbn = JosephsonCircuits.SparseArrays.sparse([1,1,2,2,3],
            [1,2,2,3,3],[1,-1,1,-1,1],3,3)
        Nmodes = 2
        Nnodal = 3*Nmodes
        Lscale = 2e-10
        cb = JosephsonCircuits.mnacoupledbranches(Mb)
        @test cb == [1,2]
        Ntot = Nnodal + length(cb)*Nmodes
        A = Matrix(JosephsonCircuits.calcAmnaind(cb, Lb, Mb, Rbn, Nmodes,
            Nnodal, Ntot, Lscale))
        B = A[1:Nnodal, Nnodal+1:end]
        C = A[Nnodal+1:end, 1:Nnodal]
        D = A[Nnodal+1:end, Nnodal+1:end]
        schur = -B*(D\C)
        ref = Lscale .* Matrix(JosephsonCircuits.calcinvLn(Lb, Mb, Rbn,
            Nmodes))
        @test isapprox(schur, ref, rtol = 1e-12)

        # a coupled branch without a self inductance is rejected
        Lbbad = JosephsonCircuits.SparseArrays.sparsevec([1], [3e-10], 3)
        @test_throws ArgumentError JosephsonCircuits.calcAmnaind(cb, Lbbad,
            Mb, Rbn, Nmodes, Nnodal, Ntot, Lscale)
    end

    @testset "mutually coupled inductors near unity coupling" begin

        # the promoted formulation keeps the system matrix entries bounded
        # as the coupling coefficient approaches one, where the nodal
        # inverse inductance entries diverge as 1/(1-k^2), and both solvers
        # remain well behaved.
        function coupled(k)
            circuit = Tuple{String,String,String,Union{Complex{Float64},Symbol,Int64}}[]
            push!(circuit,("P1","1","0",1))
            push!(circuit,("R1","1","0",:Rl))
            push!(circuit,("C1","1","2",:Cc))
            push!(circuit,("Lj1","2","0",:Lj))
            push!(circuit,("C2","2","0",:Cj))
            push!(circuit,("L1","2","0",:Ll))
            push!(circuit,("L2","3","0",:Lf))
            push!(circuit,("P2","3","0",2))
            push!(circuit,("R2","3","0",:Rl))
            push!(circuit,("K1","L1","L2",:K1))
            circuitdefs = Dict(:Lj=>500.0e-12, :Cc=>100.0e-15,
                :Cj=>1000.0e-15, :Ll=>300.0e-12, :Lf=>300.0e-12,
                :Rl=>50.0, :K1=>k)
            return circuit, circuitdefs
        end
        wp = (2*pi*4.75001*1e9,)
        sources = [(mode=(1,),port=1,current=1.0e-6),
            (mode=(1,),port=2,current=20.0e-6)]
        Sprev = nothing
        for k in (0.9, 0.9999)
            circuit, circuitdefs = coupled(k)
            d = JosephsonCircuits.hbnlsolve(wp, (4,), sources, circuit,
                circuitdefs; debugJacobian=true)
            @test length(d.coupledbranches) == 2
            # bounded entries independent of k
            @test maximum(abs, d.invLnm.nzval) < 100
            out = JosephsonCircuits.hbnlsolve(wp, (4,), sources, circuit,
                circuitdefs; ftol = 1e-10)
            @test out.solverinfo.converged
            # warm restarting from the nodal solution converges immediately
            # to the same answer (the coupled inductor auxiliary currents
            # are initialized from the constitutive relations)
            out2 = JosephsonCircuits.hbnlsolve(wp, (4,), sources, circuit,
                circuitdefs; ftol = 1e-10, x0 = out.nodeflux)
            @test out2.solverinfo.converged
            @test isapprox(out2.nodeflux[:], out.nodeflux[:], atol = 1e-8)
            # the linearized solve produces finite, passive scattering
            ws = 2*pi*(4.5:0.25:5.0)*1e9
            hb = JosephsonCircuits.hbsolve(ws, wp, sources, (2,), (4,),
                circuit, circuitdefs)
            S = hb.linearized.S((0,),1,(0,),1,:)
            @test all(isfinite, S)
            @test all(abs.(S) .<= 1 + 1e-6)
        end

        # perfect coupling (k = 1 exactly): the nodal inverse inductance
        # does not exist, but the promoted formulation keeps the branch
        # inductance matrix un-inverted and the system remains well posed.
        # the flux biased rf-SQUID scalar equation still holds, and the
        # solution is continuous with the k -> 1 limit.
        phi0 = JosephsonCircuits.phi0
        L1 = 300.0e-12
        Lj1v = 500.0e-12
        ck1 = Tuple{String,String,String,Union{Complex{Float64},Symbol,Int64}}[]
        push!(ck1,("P1","1","0",1))
        push!(ck1,("R1","1","0",:Rl))
        push!(ck1,("L1","1","0",:Ll))
        push!(ck1,("Lj1","1","0",:Lj))
        push!(ck1,("C1","1","0",:Cs))
        push!(ck1,("L2","2","0",:Ll))
        push!(ck1,("P2","2","0",2))
        push!(ck1,("R2","2","0",:Rl))
        push!(ck1,("K1","L1","L2",:K1))
        dk1 = Dict(:Rl=>50.0, :Ll=>L1, :Lj=>Lj1v, :Cs=>100.0e-15, :K1=>1.0)
        Iflux = 2.0e-6
        outk1 = JosephsonCircuits.hbnlsolve((2*pi*5.0e9,), (2,),
            [(mode=(0,),port=2,current=Iflux)], ck1, dk1;
            dc = true, odd = true, even = false, ftol = 1e-12)
        @test outk1.solverinfo.converged
        phik1 = real(outk1.nodeflux[1])
        @test isapprox(phi0*phik1 + (L1*phi0/Lj1v)*sin(phik1), L1*Iflux,
            rtol = 1e-10)
        # linearized solve at k = 1
        ws = 2*pi*(4.5:0.25:5.0)*1e9
        hbk1 = JosephsonCircuits.hbsolve(ws, (2*pi*4.75001*1e9,),
            [(mode=(1,),port=1,current=1.0e-6),
             (mode=(1,),port=2,current=20.0e-6)], (2,), (4,), ck1, dk1)
        @test all(isfinite, hbk1.linearized.S)

        # mutual coupling between inductors sharing a single branch is
        # rejected with an informative error, because the parallel
        # inductors are combined into one branch inductance before the
        # coupling is applied, which silently misrepresents the pair
        csame = Tuple{String,String,String,Union{Complex{Float64},Symbol,Int64}}[]
        push!(csame,("P1","1","0",1))
        push!(csame,("R1","1","0",:Rl))
        push!(csame,("L1","1","0",:Ll))
        push!(csame,("L2","1","0",:Ll))
        push!(csame,("C1","1","0",:Cs))
        push!(csame,("K1","L1","L2",:K1))
        errsb = try
            JosephsonCircuits.hbnlsolve((2*pi*5.0e9,), (1,),
                [(mode=(1,),port=1,current=1.0e-6)], csame,
                Dict(:Rl=>50.0, :Ll=>L1, :Cs=>100.0e-15, :K1=>0.5))
            nothing
        catch e
            e
        end
        @test errsb isa ArgumentError
        @test occursin("same branch", sprint(showerror, errsb))

        # an uncoupled inductor sharing a branch with a coupled inductor
        # is also rejected: the merged branch inductance misrepresents the
        # effective mutual coupling (by the current division factor)
        cmixed = Tuple{String,String,String,Union{Complex{Float64},Symbol,Int64}}[]
        push!(cmixed,("P1","1","0",1))
        push!(cmixed,("R1","1","0",:Rl))
        push!(cmixed,("L1","1","0",:Ll))
        push!(cmixed,("L2","1","0",:Ll))
        push!(cmixed,("L3","2","0",:Ll))
        push!(cmixed,("P2","2","0",2))
        push!(cmixed,("R2","2","0",:Rl))
        push!(cmixed,("C1","1","0",:Cs))
        push!(cmixed,("K1","L1","L3",:K1))
        errmx = try
            JosephsonCircuits.hbnlsolve((2*pi*5.0e9,), (1,),
                [(mode=(1,),port=1,current=1.0e-6)], cmixed,
                Dict(:Rl=>50.0, :Ll=>L1, :Cs=>100.0e-15, :K1=>0.5))
            nothing
        catch e
            e
        end
        @test errmx isa ArgumentError
        @test occursin("L1", sprint(showerror, errmx))
        @test occursin("L2", sprint(showerror, errmx))
        @test occursin("mutual", sprint(showerror, errmx))

        # uncoupled parallel inductors on a shared branch remain supported
        # and combine by reciprocal sum: two 2L inductors in parallel are
        # equivalent to a single inductor L
        function parcircuit(split::Bool)
            c = Tuple{String,String,String,Union{Complex{Float64},Symbol,Int64}}[]
            push!(c,("P1","1","0",1))
            push!(c,("R1","1","0",:Rl))
            if split
                push!(c,("L1","1","0",:L2x))
                push!(c,("L2","1","0",:L2x))
            else
                push!(c,("L1","1","0",:Lx))
            end
            push!(c,("Lj1","1","0",:Lj))
            push!(c,("C1","1","0",:Cs))
            return c
        end
        dp = Dict(:Rl=>50.0, :Lx=>L1, :L2x=>2*L1, :Lj=>Lj1v, :Cs=>100.0e-15)
        o1 = JosephsonCircuits.hbnlsolve((2*pi*5.0e9,), (2,),
            [(mode=(1,),port=1,current=2.0e-6)], parcircuit(true), dp;
            ftol = 1e-12)
        o2 = JosephsonCircuits.hbnlsolve((2*pi*5.0e9,), (2,),
            [(mode=(1,),port=1,current=2.0e-6)], parcircuit(false), dp;
            ftol = 1e-12)
        @test o1.solverinfo.converged && o2.solverinfo.converged
        @test isapprox(o1.nodeflux[:], o2.nodeflux[:], atol = 1e-10)

        # the flux-pumping regression: a floating loop biased through a
        # mutual inductor at dc still solves with a coupled branch promoted
        circuit = Tuple{String,String,String,Union{Complex{Float64},Symbol,Int64}}[]
        push!(circuit,("P1","1","0",1))
        push!(circuit,("R1","1","0",:Rl))
        push!(circuit,("L1","1","0",:Lm))
        push!(circuit,("K1","L1","L2",:K1))
        push!(circuit,("C1","1","2",:Cc))
        push!(circuit,("L2","2","3",:Lm))
        push!(circuit,("Lj3","3","0",:Lj))
        push!(circuit,("Lj4","2","0",:Lj))
        push!(circuit,("C2","2","0",:Cj))
        circuitdefs = Dict(:Lj=>2000.0e-12, :Lm=>10.0e-12,
            :Cc=>200.0e-15, :Cj=>900.0e-15, :Rl=>50.0, :K1=>0.999)
        out = JosephsonCircuits.hbnlsolve((2*pi*5.0e9,), (4,),
            [(mode=(1,),port=1,current=1.0e-6)], circuit, circuitdefs;
            dc = true, odd = true, even = true, ftol = 1e-10)
        @test out.solverinfo.converged
    end

    @testset "junction and inductor on a shared branch" begin

        # a Josephson junction and a linear inductor between the same two
        # nodes share a branch flux, and their currents add in the
        # Kirchhoff current law equations: the physical parallel
        # combination (an rf-SQUID). These tests lock the correctness of
        # this configuration, including with dc and microwave currents
        # simultaneously and with the shared inductor promoted as a
        # mutually coupled branch.
        phi0 = JosephsonCircuits.phi0
        L = 300.0e-12
        Lj = 500.0e-12
        C = 100.0e-15
        Rl = 50.0
        circuit = Tuple{String,String,String,Union{Complex{Float64},Symbol,Int64}}[]
        push!(circuit,("P1","1","0",1))
        push!(circuit,("R1","1","0",:Rl))
        push!(circuit,("L1","1","0",:Ll))
        push!(circuit,("Lj1","1","0",:Lj))
        push!(circuit,("C1","1","0",:Cs))
        circuitdefs = Dict(:Rl=>Rl, :Ll=>L, :Lj=>Lj, :Cs=>C)

        # dc bias: the node flux satisfies the scalar rf-SQUID equation
        # Idc = (phi0/L)*phi + (phi0/Lj)*sin(phi)
        Idc = 1.0e-6
        out = JosephsonCircuits.hbnlsolve((2*pi*5.0e9,), (1,),
            [(mode=(0,),port=1,current=Idc)], circuit, circuitdefs;
            dc = true, odd = true, even = false, ftol = 1e-12)
        @test out.solverinfo.converged
        phidc = real(out.nodeflux[1])
        @test isapprox(phi0*phidc/L + phi0*sin(phidc)/Lj, Idc,
            rtol = 1e-10)

        # linearized about the dc bias: S11 equals the lossless resonator
        # with the parallel effective inductance
        # Leff = 1/(1/L + cos(phidc)/Lj), with the port resistor as the
        # scattering reference
        ws = 2*pi*(3.03:0.5:7.9)*1e9
        hb = JosephsonCircuits.hbsolve(ws, (2*pi*5.0e9,),
            [(mode=(0,),port=1,current=Idc)], (1,), (1,), circuit,
            circuitdefs; dc = true, threewavemixing = true,
            fourwavemixing = false)
        # linearize the analytic model about the operating point hbsolve
        # itself found, so the comparison is self consistent with its
        # convergence tolerance
        phidclin = real(hb.nonlinear.nodeflux[1])
        Leff = 1/(1/L + cos(phidclin)/Lj)
        for (i, w) in enumerate(ws)
            ZLC = 1/(1/(im*w*Leff) + im*w*C)
            S11a = (ZLC - Rl)/(ZLC + Rl)
            @test isapprox(hb.linearized.S[1,1,1,1,i], S11a, atol = 1e-12)
        end

        # simultaneous dc and weak microwave drive: the fundamental
        # response matches the linear response at the effective inductance
        # to first order
        Irf = 5.0e-9
        out2 = JosephsonCircuits.hbnlsolve((2*pi*5.0e9,), (4,),
            [(mode=(0,),port=1,current=Idc),
             (mode=(1,),port=1,current=Irf)], circuit, circuitdefs;
            dc = true, odd = true, even = true, ftol = 1e-12)
        @test out2.solverinfo.converged
        w = 2*pi*5.0e9
        Z = 1/(1/Rl + 1/(im*w*Leff) + im*w*C)
        @test isapprox(im*w*out2.nodeflux[2]*phi0, Irf*Z, rtol = 1e-4)

        # the shared inductor promoted as a mutually coupled branch: a
        # flux-biased rf-SQUID satisfies
        # phi0*phi + (L1*phi0/Lj)*sin(phi) - M*Iflux = 0, since the dc
        # flux line current is the full drive (resistors carry no direct
        # current in the periodic steady state)
        k = 0.6
        M = k*L
        circuit2 = Tuple{String,String,String,Union{Complex{Float64},Symbol,Int64}}[]
        push!(circuit2,("P1","1","0",1))
        push!(circuit2,("R1","1","0",:Rl))
        push!(circuit2,("L1","1","0",:Ll))
        push!(circuit2,("Lj1","1","0",:Lj))
        push!(circuit2,("C1","1","0",:Cs))
        push!(circuit2,("L2","2","0",:Ll))
        push!(circuit2,("P2","2","0",2))
        push!(circuit2,("R2","2","0",:Rl))
        push!(circuit2,("K1","L1","L2",:K1))
        circuitdefs2 = Dict(:Rl=>Rl, :Ll=>L, :Lj=>Lj, :Cs=>C, :K1=>k)
        Iflux = 2.0e-6
        out3 = JosephsonCircuits.hbnlsolve((2*pi*5.0e9,), (2,),
            [(mode=(0,),port=2,current=Iflux)], circuit2, circuitdefs2;
            dc = true, odd = true, even = false, ftol = 1e-12)
        @test out3.solverinfo.converged
        phif = real(out3.nodeflux[1])
        @test isapprox(phi0*phif + (L*phi0/Lj)*sin(phif), M*Iflux,
            rtol = 1e-10)

        # two junctions on the same branch are rejected (an existing
        # conservative guard in the branch vector construction)
        circuit3 = Tuple{String,String,String,Union{Complex{Float64},Symbol,Int64}}[]
        push!(circuit3,("P1","1","0",1))
        push!(circuit3,("R1","1","0",:Rl))
        push!(circuit3,("Lj1","1","0",:Lj))
        push!(circuit3,("Lj2","1","0",:Lj))
        push!(circuit3,("C1","1","0",:Cs))
        @test_throws Exception JosephsonCircuits.hbnlsolve((2*pi*5.0e9,),
            (1,), [(mode=(1,),port=1,current=1.0e-6)], circuit3,
            circuitdefs)
    end

end
