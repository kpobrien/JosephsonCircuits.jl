using JosephsonCircuits
using LinearAlgebra
using SparseArrays
using Test

@testset verbose=true "hbsolve" begin


    @testset "hbsolve-hbnlsolve comparison" begin

        atol = 5e-17

        JosephsonCircuits.@params R Cc Lj Cj
        circuit = Circuit([
            ("P1", "1", "0", Port(1; Z0 = R)),
            ("C1", "1", "2", Capacitor(Cc)),
            ("Lj1", "2", "0", JosephsonJunction(Lj)),
            ("C2", "2", "0", Capacitor(Cj))])

        for tandelta in [0,1e-3]

            circuitdefs = Dict(
                Lj =>1000.0e-12,
                Cc => 100.0e-15,
                Cj => 1000.0e-15/(1+im*tandelta),
                R => 50.0)

            ws = 2*pi*4.74*1e9
            wp = 2*pi*4.75001*1e9
            Ip = 0.00565e-6
            Is = 1e-14

            # linearized simulation
            sources = [(mode=(1,),port=1,current=Ip)]
            Npumpharmonics = (6,)
            Nmodulationharmonics = (6,)
            sol1 = hbsolve(ws, (wp,), sources, Nmodulationharmonics,
                Npumpharmonics, circuit, circuitdefs, atol = atol)
            S1ss = sol1.linearized.S((0,),1,(0,),1,1)
            S1is = sol1.linearized.S((-2,),1,(0,),1,1)

            # nonlinear simulation with (pump,signal) order
            w = (wp,ws)
            Nharmonics = (6,6)
            sources = [(mode=(1,0),port=1,current=Ip),(mode=(0,1),port=1,current=Is)]
            sol2 = hbnlsolve(w, Nharmonics, sources, circuit, circuitdefs, atol = atol)
            S2ss = sol2.S((0,1),1,(0,1),1)
            S2is = sol2.S((2,-1),1,(0,1),1)

            # nonlinear simulation with (signal,pump) order
            w = (ws,wp)
            Nharmonics = (6,6)
            sources = [(mode=(0,1),port=1,current=Ip),(mode=(1,0),port=1,current=Is)]
            sol3 = hbnlsolve(w, Nharmonics, sources, circuit, circuitdefs, atol = atol)
            S3ss = sol3.S((1,0),1,(1,0),1)
            S3is = sol3.S((1,-2),1,(1,0),1)

            @test(isapprox(S1ss,S2ss))
            # conjugate the idler from this simulation since it has a positive
            # frequency  wi = 2wp - ws where wp > ws. it has a negative
            # frequency in the linearized simulation wi = ws - 2wp.
            @test(isapprox(S1is,conj(S2is)))
            @test(isapprox(S1ss,S3ss))
            # no need to conjugate the idler here since both have negative
            # frequencies.
            @test(isapprox(S1is,S3is))
        end

    end

    @testset "uncommon options give the same numbers" begin

        JosephsonCircuits.@params Rleft Cc Lj Cj w
        circuit = Any[]
        push!(circuit,("P1", "1", "0", Port(1; Z0 = Rleft)))
        push!(circuit,("C1", "1", "2", Capacitor(Cc)))
        push!(circuit,("Lj1", "2", "0", JosephsonJunction(Lj)))
        push!(circuit,("C2", "2", "0", Capacitor(Cj)))
        circuit = Circuit(circuit)
        circuitdefs = Dict(Lj => 1000.0e-12, Cc => 100.0e-15,
            Cj => 1000.0e-15, Rleft => 50.0)
        ws = 2*pi*(4.5:0.05:5.0)*1e9
        wp = (2*pi*4.75001*1e9,)
        sources = [(mode=(1,),port=1,current=0.00565e-6)]

        base = hbsolve(ws, wp, sources, (8,), (16,), circuit, circuitdefs;
            atol = 1e-12, returnnodeflux = true, returnvoltage = true,
            keyedarrays = false)
        # four batches and every output flag the other way round must give
        # the same numbers where both computed them, and nothing where they
        # were not asked for
        other = hbsolve(ws, wp, sources, (8,), (16,), circuit, circuitdefs;
            atol = 1e-12, returnS = false,
            returnSnoise = true, returnQE = false, returnnodeflux = true,
            returnnodefluxadjoint = true, returnCM = false,
            returnvoltage = true, returnvoltageadjoint = true,
            returnSsensitivity = true, sensitivitynames = ["C1"],
            nbatches = 4, keyedarrays = false)
        @test isempty(other.linearized.QE)
        @test isempty(other.linearized.S)
        @test isempty(other.linearized.CM)
        @test isempty(other.linearized.Snoise)      # a lossless circuit
        @test isapprox(other.nonlinear.nodeflux, base.nonlinear.nodeflux;
            rtol = 1e-10)
        @test isapprox(other.linearized.nodeflux, base.linearized.nodeflux;
            rtol = 1e-10)
        @test isapprox(other.linearized.voltage, base.linearized.voltage;
            rtol = 1e-10)
        @test size(other.linearized.Ssensitivity, 3) == 1
        @test all(isfinite, other.linearized.nodefluxadjoint)
    end

    @testset "hbsolve initial nodeflux" begin

        circuit = Any[]
        push!(circuit,("P1", "1", "0", Port(1; Z0 = :Rleft)))
        push!(circuit,("I1", "1", "0", CurrentSource(:Ipump)))
        push!(circuit,("L1", "1", "0", Inductor(:Lm)))
        push!(circuit,("K1", "L1", "L2", MutualInductor(:K1)))
        push!(circuit,("C1", "1", "2", Capacitor(:Cc)))
        push!(circuit,("L2", "2", "3", Inductor(:Lm)))
        push!(circuit,("Lj3", "3", "0", JosephsonJunction(:Lj)))
        push!(circuit,("Lj4", "2", "0", JosephsonJunction(:Lj)))
        push!(circuit,("C2", "2", "0", Capacitor(:Cj)))
        circuit = Circuit(circuit)
        circuitdefs = Dict{Symbol,Complex{Float64}}(
            :Lj =>2000e-12,
            :Lm =>10e-12,
            :Cc => 200.0e-15,
            :Cj => 900e-15,
            :Rleft => 50.0,
            :Rright => 50.0,
            :Ipump => 1.0e-8,
            :K1 => 0.9,
        )

        Idc = 50e-5
        Ip=0.0001e-6
        wp=2*pi*5e9
        Npumpmodes = 2
        out1=hbnlsolve(
            (wp,),
            (Npumpmodes,),
            [
                (mode=(0,),port=1,current=Idc),
                (mode=(1,),port=1,current=Ip),
            ],
            circuit,circuitdefs;dc=true,odd=true,even=false)
        out2=hbnlsolve(
            (wp,),
            (Npumpmodes,),
            [
                (mode=(0,),port=1,current=Idc),
                (mode=(1,),port=1,current=Ip),
            ],
            circuit,circuitdefs;dc=true,odd=true,even=false,
            x0 = out1.nodeflux[:]);
        @test isapprox(out1.nodeflux[:],out2.nodeflux[:])
    end


    @testset verbose=true "hbsolve return flags" begin

        JosephsonCircuits.@params R Cc Lj Cj
        circuit = Circuit([
            ("P1", "1", "0", Port(1; Z0 = R)),
            ("C1", "1", "2", Capacitor(Cc)),
            ("Lj1", "2", "0", JosephsonJunction(Lj)),
            ("C2", "2", "0", Capacitor(Cj))])

        circuitdefs = Dict(
            Lj =>1000.0e-12,
            Cc => 100.0e-15,
            Cj => 1000.0e-15/(1+1e-3im),
            R => 50.0)

        ws = 2*pi*(4.5:0.5:5.0)*1e9
        wp = (2*pi*4.75001*1e9,)
        Ip = 0.00565e-6
        sources = [(mode=(1,),port=1,current=Ip)]
        Npumpharmonics = (16,)
        Nmodulationharmonics = (8,)

        # these are all of the returns we will examine
        flags = ["S","Snoise","QE","CM","nodeflux","voltage","nodefluxadjoint","voltageadjoint",
            "Ssensitivity"]

        # set all of the flags to be true
        returnflags = NamedTuple([(Symbol("return"*flags[i])=>true) for i in 1:length(flags)])
        solalltrue = hbsolve(ws, wp, sources, Nmodulationharmonics,
            Npumpharmonics, circuit, circuitdefs;returnflags...);

        # every flag false leaves every output empty
        returnflags = NamedTuple([(Symbol("return"*flags[i])=>false) for i in 1:length(flags)])
        solallfalse = hbsolve(ws, wp, sources, Nmodulationharmonics,
            Npumpharmonics, circuit, circuitdefs;returnflags...);
        for k in 1:length(flags)
            @test isempty(getfield(solallfalse.linearized, Symbol(flags[k])))
        end

        # one flag true at a time, for the two outputs which need more than
        # the forward solve (the noise channels and the adjoint solve, and
        # the sensitivity contraction): the same numbers as with every flag
        # true, and nothing else
        for j in (findfirst(==("Snoise"), flags), findfirst(==("Ssensitivity"), flags))
            # set one of the flags to be true and the rest false
            returnflags = NamedTuple([(Symbol("return"*flags[i])=>ifelse(i==j,true,false)) for i in 1:length(flags)])
            sol = hbsolve(ws, wp, sources, Nmodulationharmonics,
                Npumpharmonics, circuit, circuitdefs;returnflags...);

            # loop over all of the flags and check if the returned value when all of the flags
            # are true is the same as when only the selected flag is true. check that the returned
            # values for the false flags are empty.
            for k in 1:length(flags)
                # compare whether the all flags true return value is the same as when only
                # one flag is true
                if k == j
                    result = @test(isapprox(
                        getfield(solalltrue.linearized, Symbol(flags[k])),
                        getfield(sol.linearized, Symbol(flags[k])),
                        )
                    )
                    if result isa Test.Fail
                        println("",flags[k]," is not correct when ","return"*flags[j]," = true")
                    end
                # check that the rest of the return values are empty.
                else
                    result = @test(isempty(getfield(sol.linearized, Symbol(flags[k]))))
                    if result isa Test.Fail
                        println("",flags[k]," is not empty when ","return"*flags[j]," = true")
                    end
                end
            end
        end
    end

    @testset verbose=true "hbnlsolve lossless error" begin

        JosephsonCircuits.@params Rleft Cc Lj Cj w L1
        circuit = Any[]
        push!(circuit,("P1", "1", "0", Port(1; Z0 = Rleft)))
        push!(circuit,("C1", "1", "2", Capacitor(Cc)))
        push!(circuit,("Lj1", "2", "0", JosephsonJunction(Lj)))
        push!(circuit,("C2", "2", "0", Capacitor(Cj)))
        circuit = Circuit(circuit)
        circuitdefs = Dict(
            Lj =>1000.0e-12,
            Cc => 100.0e-15,
            Cj => 1000.0e-15,
            Rleft => 50.0,
        )
        ws = 2*pi*(4.5:0.01:5.0)*1e9
        wp = 2*pi*4.75001*1e9
        Ip = 0.00565e-6
        Nsignalmodes = 8
        Npumpmodes = 8

        w = (wp,)
        Nharmonics = (2*Npumpmodes,)
        sources = ((mode=(1,),port=1,current=Ip),)

        r = @test_logs (:warn,) match_mode=:any hbnlsolve(w, Nharmonics,
            sources, circuit, circuitdefs, iterations = 1)
        @test !r.solverinfo.converged
        @test r.solverinfo.stages[end].reason == :iterations
    end

    @testset "hbnlsolve simple testcase" begin

        circuit = Circuit([("P1", "1", "0", Port(1; Z0 = 50.0))])
        circuitdefs = Dict()
        Idc = 50e-5
        Ip = 1.0e-6
        wp=2*pi*5.0*1e9
        Npumpmodes = 1
        out=hbnlsolve(
            (wp,),
            (Npumpmodes,),
            [
        #         (mode=(0,),port=1,current=Idc),
                (mode=(1,),port=1,current=Ip),
            ],
            circuit,circuitdefs;dc=false,odd=true,even=false)
        @test isapprox(im*out.nodeflux[1]*wp*JosephsonCircuits.phi0/(50),Ip)

    end

    @testset "hbnlsolve simple testcase dc" begin

        # a port and resistor with no inductive path to ground: a nodal
        # system matrix with a DC mode would be structurally singular, but
        # with the modified nodal analysis formulation the DC node flux is
        # gauge fixed and the circuit solves exactly.
        circuit = Circuit([("P1", "1", "0", Port(1; Z0 = 50.0))])
        circuitdefs = Dict()
        Idc = 50e-5
        Ip = 1.0e-6
        wp=2*pi*5.0*1e9
        Npumpmodes = 1

        out = JosephsonCircuits.hbnlsolve(
            (wp,),
            (Npumpmodes,),
            [
                (mode=(1,),port=1,current=Ip),
            ],
            circuit,circuitdefs;dc=true,odd=true,even=false)
        @test out.solverinfo.converged
        @test isapprox(out.nodeflux[1], 0.0, atol = 1e-15)
        @test isapprox(im*out.nodeflux[2]*wp*JosephsonCircuits.phi0/(50),Ip)
    end

    @testset "undefined symbolic component values" begin

        # forgetting to assign a value to a symbolic variable in
        # circuitdefs fails immediately with an ArgumentError naming the
        # component and the undefined variable, instead of a downstream
        # error about the symbolic frequency variable.
        JosephsonCircuits.@params Rv Ccv Ljv Cjv
        circuit = Circuit([("P1", "1", "0", Port(1; Z0 = Rv)),("C1", "1", "2", Capacitor(Ccv)),
            ("Lj1", "2", "0", JosephsonJunction(Ljv)),("C2", "2", "0", Capacitor(Cjv))])
        circuitdefs = Dict(Ljv=>1000.0e-12, Cjv=>1000.0e-15, Rv=>50.0)
        wp = (2*pi*4.75001*1e9,)
        sources = [(mode=(1,),port=1,current=0.00565e-6)]
        err = try
            JosephsonCircuits.hbsolve(2*pi*(4.5:0.1:5.0)*1e9, wp, sources,
                (2,), (2,), circuit, circuitdefs)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("C1", sprint(showerror, err))
        @test occursin("Ccv", sprint(showerror, err))
        @test occursin("circuitdefs", sprint(showerror, err))

        # the same check protects hblinsolve directly
        err2 = try
            JosephsonCircuits.hblinsolve(2*pi*(4.5:0.1:5.0)*1e9, circuit,
                circuitdefs)
            nothing
        catch e
            e
        end
        @test err2 isa ArgumentError
        @test occursin("Ccv", sprint(showerror, err2))

        # a frequency dependent value is accepted, and one which also
        # carries an undefined parameter is rejected naming both the
        # component and the parameter
        JosephsonCircuits.@params Rundef
        wfd = FrequencyDependent(identity)
        c2 = Circuit([("P1", "1", "0", Port(1; Z0 = 50.0 + 0.0*wfd)),("C1", "1", "0", Capacitor(100.0e-15)),("L1", "1", "0", Inductor(1.0e-9))])
        out = JosephsonCircuits.hbnlsolve(wp, (1,), sources, c2, Dict())
        @test out.solverinfo.converged
        c3 = Circuit([("P1", "1", "0", Port(1; Z0 = Rundef/(1 + wfd*1e-12))),("C1", "1", "0", Capacitor(100.0e-15)),("L1", "1", "0", Inductor(1.0e-9))])
        err3 = try
            JosephsonCircuits.hbnlsolve(wp, (1,), sources, c3, Dict())
            nothing
        catch e
            e
        end
        @test err3 isa ArgumentError
        @test occursin("Rundef", sprint(showerror, err3))
        @test occursin("P1", sprint(showerror, err3))
    end

    @testset "a value written as a closure of the frequency" begin
        # the same frequency law spelled two ways: a closure passed whole,
        # and the identity closure entering an expression. The scattering
        # parameters must agree to roundoff.
        wp = (2*pi*4.75001*1e9,)
        ws = 2*pi*(4.5:0.1:5.0)*1e9
        sources = [(mode=(1,),port=1,current=0.00565e-6)]
        law(w) = 50.0*(1 + (w/1e11)^2)
        wfd = FrequencyDependent(identity)
        cexp = Circuit([("P1", "1", "0", Port(1; Z0 = law(wfd))),("C1", "1", "2", Capacitor(100.0e-15)),("Lj1", "2", "0", JosephsonJunction(1000.0e-12)),
            ("C2", "2", "0", Capacitor(1000.0e-15))])
        cfun = Circuit([("P1", "1", "0", Port(1; Z0 = FrequencyDependent(law))),
            ("C1", "1", "2", Capacitor(100.0e-15)),("Lj1", "2", "0", JosephsonJunction(1000.0e-12)),
            ("C2", "2", "0", Capacitor(1000.0e-15))])
        Sexp = hbsolve(ws, wp, sources, (2,), (8,), cexp, Dict()).linearized.S
        Sfun = hbsolve(ws, wp, sources, (2,), (8,), cfun).linearized.S
        @test isapprox(Array(Sexp), Array(Sfun), rtol = 1e-12)

        # a lossy capacitor is a noise channel however its value is
        # written: as a plain complex number, or as a frequency dependent
        # closure. The imaginary part is only visible once the value is
        # resolved at a frequency, so both must resolve before the
        # channels are chosen.
        lossy(v) = Circuit([("P1", "1", "0", Port(1; Z0 = 50.0)),
            ("C1", "1", "2", Capacitor(v)),
            ("Lj1", "2", "0", JosephsonJunction(1000.0e-12)),
            ("C2", "2", "0", Capacitor(500.0e-15))])
        C0, tand = 100.0e-15, 1e-3
        outs = (hbsolve(ws, wp, sources, (2,), (8,), lossy(C0*(1 - im*tand)),
                    Dict(); keyedarrays = false),
            hbsolve(ws, wp, sources, (2,), (8,),
                lossy(C0*(1 - im*tand)*(1 + 0*wfd)), Dict();
                keyedarrays = false),
            hbsolve(ws, wp, sources, (2,), (8,),
                lossy(FrequencyDependent(w -> C0*(1 - im*tand))), Dict();
                keyedarrays = false))
        for o in outs
            @test [o.linearized.componentnames[i]
                for i in o.linearized.noiseportimpedanceindices] == ["C1"]
            @test isapprox(o.linearized.QE, first(outs).linearized.QE,
                rtol = 1e-12)
            @test o.linearized.QE[1,1,1] < o.linearized.QEideal[1,1,1]
        end
    end

    @testset "calcsources errors" begin

        modes = [(0,), (1,)]
        portindices = [1]
        portnumbers = [1]
        nodeindices = [2 2 2 2 0 2 3 4 3 3; 1 1 1 1 0 3 4 1 1 1]
        edge2indexdict = Dict((1, 2) => 1, (3, 1) => 2, (1, 3) => 2, (4, 1) => 3, (2, 1) => 1, (1, 4) => 3, (3, 4) => 4, (4, 3) => 4)
        Lscale = 1.005e-9 + 0.0im
        Nnodes = 4
        Nbranches = 4
        Nmodes = 2

        # current source for non-existent port
        sources = [(mode = (0,), port = 1, current = 0.0005), (mode = (1,), port = 2, current = 1.0e-10)]
        @test_throws(
            ArgumentError("Source port 2 not found."),
            JosephsonCircuits.calcsources(modes, sources, portindices, portnumbers,
                nodeindices, edge2indexdict, Lscale, Nnodes, Nbranches, Nmodes))

        # current source for non-existent mode
        sources = [(mode = (0,), port = 1, current = 0.0005), (mode = (2,), port = 1, current = 1.0e-10)]
        @test_throws(
            ArgumentError("Source mode (2,) is not among the retained modes; the truncation (`Nharmonics`, `maxintermodorder`, `dc`, `odd`, `even`, `frequencywindow`) removed it."),
            JosephsonCircuits.calcsources(modes, sources, portindices, portnumbers,
                nodeindices, edge2indexdict, Lscale, Nnodes, Nbranches, Nmodes))

        # a direct current with a phase is refused rather than stamped
        # without its imaginary part
        sources = [(mode = (0,), port = 1, current = 0.0005im), (mode = (1,), port = 1, current = 1.0e-10)]
        @test_throws ArgumentError JosephsonCircuits.calcsources(modes,
            sources, portindices, portnumbers, nodeindices, edge2indexdict,
            Lscale, Nnodes, Nbranches, Nmodes)
        sources = [(mode = (0,), port = 1, current = 0.0005 + 0.0im), (mode = (1,), port = 1, current = 1.0e-10im)]
        @test JosephsonCircuits.calcsources(modes, sources, portindices,
            portnumbers, nodeindices, edge2indexdict, Lscale, Nnodes,
            Nbranches, Nmodes) isa AbstractVector

    end

    # The adjoint solutions the noise, quantum efficiency and commutation
    # relation calculations require are the solutions of the linearized
    # system with the complex conjugate of the pump modulation matrix. That
    # system is a diagonal similarity transformation of the transposed
    # forward system, so hblinsolve obtains the adjoint solutions with a
    # transposed solve on the factorization of the forward system instead of
    # assembling and factorizing the conjugated pump matrix at every signal
    # frequency. These tests check the similarity relation and the
    # equivalence of the two solve strategies.
    @testset verbose=true "hblinsolve adjoint solve" begin

        JosephsonCircuits.@params Rleft Rright Cc Lj Cj Lla Llb Kab

        # a JPA: one port, so one promoted port resistor
        circuitjpa = Any[]
        push!(circuitjpa,("P1", "1", "0", Port(1; Z0 = Rleft)))
        push!(circuitjpa,("C1", "1", "2", Capacitor(Cc)))
        push!(circuitjpa,("Lj1", "2", "0", JosephsonJunction(Lj)))
        push!(circuitjpa,("C2", "2", "0", Capacitor(Cj)))
        circuitjpa = Circuit(circuitjpa)
        circuitdefsjpa = Dict(Lj=>1000.0e-12, Cc=>100.0e-15, Cj=>1000.0e-15,
            Rleft=>50.0)

        # a lossy JPA: the complex capacitance adds a noise port, which is the
        # consumer of the adjoint solution we most care about here
        circuitdefsjpalossy = Dict(Lj=>1000.0e-12, Cc=>100.0e-15,
            Cj=>1000.0e-15/(1+1e-3im), Rleft=>50.0)

        # two ports and a mutually coupled inductor pair, which is promoted to
        # auxiliary branch currents as well. exercises both auxiliary blocks at
        # once, with the coupling coefficient close to one.
        circuitmutual = Any[]
        push!(circuitmutual,("P1", "1", "0", Port(1; Z0 = Rleft)))
        push!(circuitmutual,("C1", "1", "2", Capacitor(Cc)))
        push!(circuitmutual,("Lj1", "2", "0", JosephsonJunction(Lj)))
        push!(circuitmutual,("C2", "2", "0", Capacitor(Cj)))
        push!(circuitmutual,("L1", "2", "0", Inductor(Lla)))
        push!(circuitmutual,("L2", "3", "0", Inductor(Llb)))
        push!(circuitmutual,("P2", "3", "0", Port(2; Z0 = Rright)))
        push!(circuitmutual,("K1", "L1", "L2", MutualInductor(Kab)))
        circuitmutual = Circuit(circuitmutual)
        circuitdefsmutual = Dict(Lj=>500.0e-12, Cc=>100.0e-15, Cj=>1000.0e-15,
            Lla=>300.0e-12, Llb=>300.0e-12, Rleft=>50.0, Rright=>50.0,
            Kab=>0.99)

        testcases = (
            ("single-tone JPA", (2*pi*4.75001e9,),
                [(mode=(1,),port=1,current=0.00565e-6)], (4,), (2,),
                circuitjpa, circuitdefsjpa),
            ("single-tone lossy JPA", (2*pi*4.75001e9,),
                [(mode=(1,),port=1,current=0.00565e-6)], (4,), (2,),
                circuitjpa, circuitdefsjpalossy),
            ("two-tone JPA", (2*pi*4.65001e9, 2*pi*4.85001e9),
                [(mode=(1,0),port=1,current=0.00565e-6*1.7),
                 (mode=(0,1),port=1,current=0.00565e-6*1.7)], (4,4), (2,2),
                circuitjpa, circuitdefsjpa),
            ("mutual inductor", (2*pi*4.75001e9,),
                [(mode=(1,),port=1,current=1.0e-6)], (4,), (2,),
                circuitmutual, circuitdefsmutual),
            )

        for (name, wp, sources, Npumpharmonics, Nmodulationharmonics, circuit,
                circuitdefs) in testcases
            @testset "$name" begin

                ws = 2*pi*[4.5e9, 4.75e9]

                # the pump operating point and the linearized system it defines
                nonlinear = hbnlsolve(wp, Npumpharmonics, sources, circuit,
                    circuitdefs; keyedarrays=false)
                psc = JosephsonCircuits.compile(circuit)
                signalfreq = JosephsonCircuits.truncfreqs(
                    JosephsonCircuits.calcfreqsdft(Nmodulationharmonics);
                    dc=true, odd=false, even=true, maxintermodorder=Inf)
                d = JosephsonCircuits.hblinsolve(ws, psc, circuitdefs,
                    signalfreq; nonlinear=nonlinear, debuglsys=true)
                lsys = d.lsys

                # the diagonal of the similarity transformation: one on the node
                # flux rows, the assembled conductance entry of the constitutive
                # equation on the auxiliary rows of the promoted resistors, and
                # one on the auxiliary rows of the promoted coupled inductors.
                # the promoted resistances are constant and real, so the negative
                # frequency conjugation of sparseaddconjsubst!, which acts on the
                # stored conductance and not on the frequency factor, is trivial.
                function similaritydiagonal(d, wmodes)
                    D = ones(Complex{Float64}, d.Nnodalmna + d.Nauxmna)
                    return Diagonal(D)
                end

                for wsi in ws
                    wmodes = wsi .+ d.wpumpmodes
                    A = copy(lsys.Asparse)
                    Aconj = copy(lsys.Asparse)
                    JosephsonCircuits.assemblesystemmatrix!(A, lsys, wsi)
                    JosephsonCircuits.assemblesystemmatrix!(Aconj, lsys, wsi;
                        conjugatepump = true)
                    D = similaritydiagonal(d, wmodes)

                    # the documented similarity relation
                    @test isapprox(Matrix(Aconj), D*Matrix(transpose(A))*inv(D),
                        rtol = 1e-10, norm = v->maximum(abs,v))

                    # the solutions agree exactly in the node flux rows, which
                    # are the only rows the noise, quantum efficiency and
                    # adjoint output calculations read, and differ by the
                    # similarity diagonal in the auxiliary rows.
                    xconj = Matrix(Aconj)\Matrix(d.bnm)
                    xtrans = Matrix(transpose(A))\Matrix(d.bnm)
                    nodal = 1:d.Nnodalmna
                    @test isapprox(xconj[nodal,:], xtrans[nodal,:], rtol = 1e-8,
                        norm = v->maximum(abs,v))
                    @test isapprox(xconj, D*xtrans, rtol = 1e-8,
                        norm = v->maximum(abs,v))
                end

                # end to end: the transposed solve used by hblinsolve gives the
                # same adjoint node fluxes as an independent solve of the
                # conjugated pump system.
                sol = JosephsonCircuits.hblinsolve(ws, psc, circuitdefs,
                    signalfreq; nonlinear=nonlinear, keyedarrays=false,
                    returnnodefluxadjoint=true, returnSnoise=true, returnQE=true)
                for (i, wsi) in enumerate(ws)
                    Aconj = copy(lsys.Asparse)
                    JosephsonCircuits.assemblesystemmatrix!(Aconj, lsys, wsi;
                        conjugatepump = true)
                    xconj = Matrix(Aconj)\Matrix(d.bnm)
                    @test isapprox(sol.nodefluxadjoint[:,:,i],
                        xconj[1:size(sol.nodefluxadjoint,1),:], rtol = 1e-8,
                        norm = v->maximum(abs,v))
                end
            end
        end
    end


    @testset verbose=true "scattering parameter sensitivities" begin
        @testset "integer component values" begin
            # a value written as an integer assembles in the floating point
            # storage its group takes, in the sensitivity stamps as in the
            # matrices, so a reciprocal never asks an integer to hold a
            # fraction
            mk(L, R) = Circuit([(:p, 1, 0, Port(1)), (:l, 1, 0, Inductor(L)),
                (:r, 1, 0, Resistor(R)), (:c, 1, 0, Capacitor(1e-12))])
            wsi = [2pi*5e9]
            si = hblinsolve(wsi, mk(2, 75); keyedarrays = false,
                sensitivitynames = ["l", "r"], returnSsensitivity = true)
            sf = hblinsolve(wsi, mk(2.0, 75.0); keyedarrays = false,
                sensitivitynames = ["l", "r"], returnSsensitivity = true)
            @test si.S == sf.S
            @test si.Ssensitivity == sf.Ssensitivity
            # and through the operating point of a pumped circuit
            jj(L) = Circuit([(:p, 1, 0, Port(1)), (:cc, 1, 2, Capacitor(100e-15)),
                (:jj, 2, 0, JosephsonJunction(1e-9)), (:l, 2, 0, Inductor(L)),
                (:c, 2, 0, Capacitor(1e-12))])
            wpi = (2pi*4.75e9,); srci = [(mode = (1,), port = 1, current = 1e-8)]
            hi = hbsolve(wsi, wpi, srci, (1,), (2,), jj(1); keyedarrays = false,
                sensitivitynames = ["l"], returnSsensitivity = true)
            hf = hbsolve(wsi, wpi, srci, (1,), (2,), jj(1.0); keyedarrays = false,
                sensitivitynames = ["l"], returnSsensitivity = true)
            @test hi.linearized.Ssensitivity == hf.linearized.Ssensitivity
        end

        @testset "components are named as the frontend names them" begin
            # the typed frontend names a component with a symbol, so the
            # sensitivity names take one as readily as a string
            c = Circuit([(:p, 1, 0, Port(1; termination = nothing)),
                (:l, 1, 0, Inductor(2.0)), (:r, 1, 0, Resistor(50.0)),
                (:c, 1, 0, Capacitor(1e-12))])
            ws1 = [2pi*5e9]
            bystring = hblinsolve(ws1, c; keyedarrays = false,
                sensitivitynames = ["l", "r"], returnSsensitivity = true)
            bysymbol = hblinsolve(ws1, c; keyedarrays = false,
                sensitivitynames = [:l, :r], returnSsensitivity = true)
            @test bystring.Ssensitivity == bysymbol.Ssensitivity
            @test bystring.sensitivitynames == bysymbol.sensitivitynames
            # a name the circuit does not have is refused the same way by
            # every entry point, naming it
            wp1 = (2pi*4.75e9,); src1 = [(mode = (1,), port = 1, current = 1e-8)]
            for f in (() -> hblinsolve(ws1, c; sensitivitynames = [:nope],
                          returnSsensitivity = true),
                      () -> hbsolve(ws1, wp1, src1, (1,), (2,), c;
                          sensitivitynames = [:nope], returnSsensitivity = true),
                      () -> JosephsonCircuits.hbnlsolve(wp1, (2,), src1, c;
                          sensitivitynames = ["nope"]))
                e = try f(); nothing catch e; e end
                @test e isa ArgumentError
                @test occursin("nope", sprint(showerror, e))
            end
        end


        # dS/dr, the derivative of the scattering matrix with respect to a
        # relative perturbation of a component value at a fixed pump
        # operating point, against central finite differences. The pump
        # operating point is held fixed in the finite differences to match
        # the definition, by reusing one nonlinear solution while perturbing
        # the component values of the linearized solve.

        @testset "linear network" begin
            JosephsonCircuits.@params R1v R2v R3v C1v L1v C2v
            circuit = Any[]
            push!(circuit,("P1", "1", "0", Port(1; Z0 = R1v)))
            push!(circuit,("C1", "1", "2", Capacitor(C1v))); push!(circuit,("L1", "2", "0", Inductor(L1v)))
            push!(circuit,("C2", "2", "0", Capacitor(C2v))); push!(circuit,("P2", "2", "0", Port(2; Z0 = R2v)))
            push!(circuit,("R3", "1", "2", Resistor(R3v)))
            circuit = Circuit(circuit)
            defs = Dict(R1v=>50.0, R2v=>50.0, R3v=>300.0, C1v=>100e-15,
                L1v=>1e-9, C2v=>200e-15)
            ws = 2*pi*[5.0e9, 7.0e9]
            names = ["C1","L1","C2","R3","P1/termination","P2/termination"]
            syms = Dict("C1"=>C1v,"L1"=>L1v,"C2"=>C2v,"R3"=>R3v,
                "P1/termination"=>R1v,"P2/termination"=>R2v)
            sol = hblinsolve(ws, circuit, defs; keyedarrays=false,
                sensitivitynames=names, returnSsensitivity=true)
            @test size(sol.Ssensitivity) ==
                (size(sol.S,1), size(sol.S,2), length(names), length(ws))
            h = 1e-6
            for (k, name) in enumerate(names)
                dp = copy(defs); dp[syms[name]] *= (1+h)
                dm = copy(defs); dm[syms[name]] *= (1-h)
                Sp = hblinsolve(ws, circuit, dp; keyedarrays=false).S
                Sm = hblinsolve(ws, circuit, dm; keyedarrays=false).S
                fd = (Sp .- Sm)./(2*h)
                for wi in eachindex(ws)
                    @test isapprox(sol.Ssensitivity[:,:,k,wi], fd[:,:,wi],
                        rtol = 1e-6, norm = v->maximum(abs,v))
                end
            end
        end

        @testset "pumped junction" begin
            JosephsonCircuits.@params Rl Cc Lj Cj
            circuit = Any[]
            push!(circuit,("P1", "1", "0", Port(1; Z0 = Rl)))
            push!(circuit,("C1", "1", "2", Capacitor(Cc))); push!(circuit,("Lj1", "2", "0", JosephsonJunction(Lj)))
            push!(circuit,("C2", "2", "0", Capacitor(Cj)))
            circuit = Circuit(circuit)
            defs = Dict(Rl=>50.0, Cc=>100e-15, Lj=>1000e-12, Cj=>1000e-15)
            wp = (2*pi*4.75001e9,)
            sources = [(mode=(1,),port=1,current=0.00565e-6)]
            ws = 2*pi*[4.5e9, 4.75e9]
            names = ["C1","C2","Lj1","P1/termination"]
            syms = Dict("C1"=>Cc,"C2"=>Cj,"Lj1"=>Lj,"P1/termination"=>Rl)

            sol = hbsolve(ws, wp, sources, (4,), (8,), circuit, defs;
                keyedarrays=false, sensitivitynames=names,
                returnSsensitivity=true,sensitivityoperatingpoint=false)

            # one nonlinear solution, reused so the operating point is fixed
            nonlinear = hbnlsolve(wp, (8,), sources, circuit, defs;
                keyedarrays=false)
            psc = JosephsonCircuits.compile(circuit)
            signalfreq = JosephsonCircuits.truncfreqs(
                JosephsonCircuits.calcfreqsdft((4,)); dc=true, odd=false,
                even=true, maxintermodorder=Inf)
            frozen(d) = JosephsonCircuits.hblinsolve(ws, psc, d,
                signalfreq; nonlinear=nonlinear, keyedarrays=false).S

            h = 1e-6
            for (k, name) in enumerate(names)
                dp = copy(defs); dp[syms[name]] *= (1+h)
                dm = copy(defs); dm[syms[name]] *= (1-h)
                fd = (frozen(dp) .- frozen(dm))./(2*h)
                for wi in eachindex(ws)
                    @test isapprox(sol.linearized.Ssensitivity[:,:,k,wi],
                        fd[:,:,wi], rtol = 1e-5, norm = v->maximum(abs,v))
                end
            end

            # keyed array output round trip
            solk = hbsolve(ws, wp, sources, (4,), (8,), circuit, defs;
                sensitivitynames=names, returnSsensitivity=true,
                sensitivityoperatingpoint=false)
            modes = collect(sol.linearized.modes)
            s0 = findfirst(==((0,)), modes)
            @test isapprox(
                solk.linearized.Ssensitivity(outputmode=(0,), outputport=1,
                    inputmode=(0,), inputport=1, component="Lj1",
                    freqindex=2),
                sol.linearized.Ssensitivity[s0,s0,3,2])
        end


        @testset "operating point shift" begin
            # the total derivative, including the shift of the pump operating
            # point, against central finite differences of the full solve, in
            # which the pump is re-solved. At this operating point the
            # operating point contribution is comparable to or larger than the
            # frozen pump term, so the two must differ substantially.
            JosephsonCircuits.@params Rl Cc Lj Cj
            circuit = Any[]
            push!(circuit,("P1", "1", "0", Port(1; Z0 = Rl)))
            push!(circuit,("C1", "1", "2", Capacitor(Cc))); push!(circuit,("Lj1", "2", "0", JosephsonJunction(Lj)))
            push!(circuit,("C2", "2", "0", Capacitor(Cj)))
            circuit = Circuit(circuit)
            defs = Dict(Rl=>50.0, Cc=>100e-15, Lj=>1000e-12, Cj=>1000e-15)
            wp = (2*pi*4.75001e9,)
            sources = [(mode=(1,),port=1,current=0.00565e-6)]
            ws = 2*pi*[4.5e9, 4.75e9]
            names = ["C1","C2","Lj1","P1/termination"]
            syms = Dict("C1"=>Cc,"C2"=>Cj,"Lj1"=>Lj,"P1/termination"=>Rl)
            solve(d; op=false) = hbsolve(ws, wp, sources, (4,), (8,),
                circuit, d; keyedarrays=false, atol=1e-13,
                sensitivitynames=names, returnSsensitivity=true,
                sensitivityoperatingpoint=op)

            total = solve(defs; op=true)
            frozen = solve(defs; op=false)
            h = 1e-6
            for (k, name) in enumerate(names)
                dp = copy(defs); dp[syms[name]] *= (1+h)
                dm = copy(defs); dm[syms[name]] *= (1-h)
                fd = (solve(dp).linearized.S .- solve(dm).linearized.S)./(2*h)
                for wi in eachindex(ws)
                    @test isapprox(total.linearized.Ssensitivity[:,:,k,wi],
                        fd[:,:,wi], rtol = 1e-4, norm = v->maximum(abs,v))
                    # the frozen pump derivative is not the total derivative
                    @test !isapprox(frozen.linearized.Ssensitivity[:,:,k,wi],
                        fd[:,:,wi], rtol = 0.05, norm = v->maximum(abs,v))
                end
            end

            # the operating point is only retained when it is requested
            @test isnothing(hbnlsolve(wp, (16,), sources, circuit, defs;
                keyedarrays=false).operatingpoint)
            op = hbnlsolve(wp, (16,), sources, circuit, defs;
                keyedarrays=false, returnoperatingpoint=true).operatingpoint
            @test !isnothing(op.jacobian)
            @test size(op.jacobian,1) == size(op.jacobian,2)
            # the Jacobian is the exact real Jacobian, so it matches the
            # matrix free Jacobian-vector product at the converged solution
            JosephsonCircuits.setpoint!(op.sys, op.x)
            vr = randn(size(op.jacobian,1))
            Jvr = zeros(size(op.jacobian,1))
            JosephsonCircuits.jacobianvectorproduct!(Jvr, op.sys, vr)
            @test isapprox(Jvr, op.jacobian*vr, atol = 1e-8,
                norm = v->maximum(abs,v))

            # requesting the operating point shift without an operating point
            @test_throws ArgumentError JosephsonCircuits.hblinsolve(ws,
                circuit, defs; keyedarrays=false, sensitivitynames=["C1"],
                returnSsensitivity=true,
                sensitivitynodeflux=zeros(ComplexF64,1,1))
        end


        @testset "reverse mode contraction" begin
            # The contribution of the operating point shift can be contracted
            # in either order. The reverse order pushes the output functional
            # through the transposed pump Jacobian once per output port mode
            # pair instead of contracting each component against the full
            # sparsity structure of the linearized system, so its cost is
            # independent of the number of components. The two must agree.
            JosephsonCircuits.@params Rl Ljx Cg Cjx Cg1v
            Ncells = 3
            circuit = Any[]
            push!(circuit,("P1", "1", "0", Port(1; Z0 = Rl)))
            for i in 1:Ncells
                push!(circuit,("Lj$(i)", "$(i)", "$(i+1)", JosephsonJunction(Ljx)))
                push!(circuit,("Cj$(i)", "$(i)", "$(i+1)", Capacitor(Cjx)))
                push!(circuit,("Cg$(i)", "$(i+1)", "0", Capacitor(i == 1 ? Cg1v : Cg)))
            end
            push!(circuit,("P2", "$(Ncells+1)", "0", Port(2; Z0 = Rl)))
            circuit = Circuit(circuit)
            defs = Dict(Rl=>50.0, Ljx=>JosephsonCircuits.IctoLj(3.4e-6),
                Cg=>45e-15, Cg1v=>45e-15, Cjx=>55e-15)
            ws = 2*pi*[5.5e9, 6.4e9]
            names = ["Cg1","Cg2","Lj2","P1/termination"]

            # The transform of the pump harmonic grid is applied one
            # dimension at a time, so both orders work for any number of
            # pump tones. The two tone grid exercises a second dimension,
            # which holds the full range of harmonics rather than only the
            # non negative ones of the first.
            pumps = (
                ("one tone", (2*pi*6.0e9,),
                    [(mode=(1,),port=1,current=1.0e-6)], (4,), (8,)),
                ("two tone", (2*pi*5.9e9, 2*pi*6.1e9),
                    [(mode=(1,0),port=1,current=0.7e-6),
                     (mode=(0,1),port=1,current=0.7e-6)], (2,2), (3,3)),
                )

            for (label, wp, sources, Nsig, Npump) in pumps
                @testset "$label" begin
                    solve(d, m) = hbsolve(ws, wp, sources, Nsig, Npump,
                        circuit, d; keyedarrays=false, atol=1e-13,
                        sensitivitynames=names, returnSsensitivity=true,
                        sensitivityoperatingpoint=true, sensitivitymode=m)

                    fwd = solve(defs, :forward)
                    rev = solve(defs, :reverse)
                    @test isapprox(fwd.linearized.Ssensitivity,
                        rev.linearized.Ssensitivity, rtol = 1e-10,
                        norm = v->maximum(abs,v))

                    # both orders against central finite differences of the
                    # full solve, in which the pump is re-solved
                    h = 1e-6
                    dp = copy(defs); dp[Cg1v] *= (1+h)
                    dm = copy(defs); dm[Cg1v] *= (1-h)
                    Sp = hbsolve(ws, wp, sources, Nsig, Npump, circuit, dp;
                        keyedarrays=false, atol=1e-13).linearized.S
                    Sm = hbsolve(ws, wp, sources, Nsig, Npump, circuit, dm;
                        keyedarrays=false, atol=1e-13).linearized.S
                    fd = (Sp .- Sm)./(2*h)
                    for wi in eachindex(ws)
                        for sol in (fwd, rev)
                            @test isapprox(
                                sol.linearized.Ssensitivity[:,:,1,wi],
                                fd[:,:,wi], rtol = 1e-4,
                                norm = v->maximum(abs,v))
                        end
                    end

                    # :auto agrees with whichever order it selects
                    @test isapprox(solve(defs, :auto).linearized.Ssensitivity,
                        fwd.linearized.Ssensitivity, rtol = 1e-10,
                        norm = v->maximum(abs,v))

                    # an unknown mode is rejected
                    @test_throws ArgumentError solve(defs, :sideways)
                end
            end
        end

        @testset "operating point input validation" begin
            # malformed low level inputs must be rejected at the boundary,
            # not discovered as out of bounds indexing inside the
            # contractions (the contraction loops are @inbounds).
            JosephsonCircuits.@params Rl Cc Lj Cj
            circuit = Any[]
            push!(circuit,("P1", "1", "0", Port(1; Z0 = Rl)))
            push!(circuit,("C1", "1", "2", Capacitor(Cc))); push!(circuit,("Lj1", "2", "0", JosephsonJunction(Lj)))
            push!(circuit,("C2", "2", "0", Capacitor(Cj)))
            circuit = Circuit(circuit)
            defs = Dict(Rl=>50.0, Cc=>100e-15, Lj=>1000e-12, Cj=>200e-15)
            wp = (2*pi*5e9,)
            sources = [(mode=(1,),port=1,current=1e-7)]
            ws = 2*pi*[4.5e9]
            psc = compile(circuit)
            nl = hbnlsolve(wp, (4,), sources, circuit, defs;
                returnoperatingpoint=true)
            op = nl.operatingpoint
            signalfreq = JosephsonCircuits.truncfreqs(
                JosephsonCircuits.calcfreqsdft((2,)); dc=true, odd=false,
                even=true)
            names = ["C2"]
            good = JosephsonCircuits.calcresidualsensitivity(op, psc, JosephsonCircuits.numericmatrices(psc, defs,
                    Nmodes=length(nl.modes)),
                [psc.componentnamedict[n] for n in names])
            lin(; kwargs...) = hblinsolve(ws, psc, defs, signalfreq;
                nonlinear=nl, keyedarrays=false, sensitivitynames=names,
                returnSsensitivity=true, kwargs...)
            # the residual derivatives are sparse from construction and
            # scale with the touched entries, not with Nstate*Ncomponents
            @test good isa JosephsonCircuits.SparseArrays.SparseMatrixCSC
            @test JosephsonCircuits.SparseArrays.nnz(good) < length(good)
            # the correctly sized input works in both orders
            for mode in (:forward, :reverse)
                @test all(isfinite, lin(sensitivityresidual=good,
                    sensitivitymode=mode).Ssensitivity)
            end
            # wrong column counts (too many and too few), wrong row count
            @test_throws DimensionMismatch lin(
                sensitivityresidual=hcat(good, good))
            @test_throws DimensionMismatch lin(
                sensitivityresidual=good[:, 1:0])
            @test_throws DimensionMismatch lin(
                sensitivityresidual=vcat(good, good))
            # wrong sizes for the node flux derivatives
            @test_throws DimensionMismatch lin(
                sensitivitynodeflux=zeros(Complex{Float64}, 3, 1))
            # both inputs at once is ambiguous between the orders
            @test_throws ArgumentError lin(sensitivityresidual=good,
                sensitivitynodeflux=zeros(Complex{Float64},
                    length(op.x), 1))
        end

        @testset "no junction operating point sensitivity" begin
            # a purely linear circuit: the linearized matrix does not depend
            # on the operating point, so the total sensitivity equals the
            # fixed operating point sensitivity, in every contraction mode,
            # and the operating point machinery must not be constructed at
            # all.
            JosephsonCircuits.@params Rl Ll Cs
            circuit = Any[]
            push!(circuit,("P1", "1", "0", Port(1; Z0 = Rl)))
            push!(circuit,("L1", "1", "2", Inductor(Ll))); push!(circuit,("C1", "2", "0", Capacitor(Cs)))
            push!(circuit,("P2", "2", "0", Port(2; Z0 = Rl)))
            circuit = Circuit(circuit)
            defs = Dict(Rl=>50.0, Ll=>300e-12, Cs=>300e-15)
            wp = (2*pi*6e9,)
            sources = [(mode=(1,),port=1,current=1e-7)]
            ws = 2*pi*[5e9]
            names = ["C1","L1"]
            syms = Dict("C1"=>Cs,"L1"=>Ll)
            solve(d; kwargs...) = hbsolve(ws, wp, sources, (2,), (2,),
                circuit, d; keyedarrays=false, sensitivitynames=names,
                returnSsensitivity=true, kwargs...)
            frozen = solve(defs; sensitivityoperatingpoint=false)
            for mode in (:forward, :reverse, :auto)
                sol = solve(defs; sensitivityoperatingpoint=true,
                    sensitivitymode=mode)
                @test sol.linearized.Ssensitivity ==
                    frozen.linearized.Ssensitivity
            end
            h = 1e-6
            for (k, name) in enumerate(names)
                dp = copy(defs); dp[syms[name]] *= (1+h)
                dm = copy(defs); dm[syms[name]] *= (1-h)
                Sp = solve(dp; sensitivityoperatingpoint=true).linearized.S
                Sm = solve(dm; sensitivityoperatingpoint=true).linearized.S
                fd = (Sp .- Sm)./(2*h)
                @test isapprox(frozen.linearized.Ssensitivity[:,:,k,1],
                    fd[:,:,1], rtol = 1e-6, norm = v->maximum(abs,v))
            end
        end

        @testset "dc pumped reverse contraction" begin
            # a dc bias current plus a pump, so the pump grid contains the
            # self-conjugate (0,) mode, which exercises the real
            # representation branch of the reverse contraction and the
            # transposed transform on a grid with a dc harmonic. the bias is
            # well below the junction critical current so the finite
            # difference re-solves stay on the same solution branch.
            JosephsonCircuits.@params Rl Ll Lj Cj
            circuit = Any[]
            push!(circuit,("P1", "1", "0", Port(1; Z0 = Rl)))
            push!(circuit,("L1", "1", "2", Inductor(Ll))); push!(circuit,("Lj1", "2", "0", JosephsonJunction(Lj)))
            push!(circuit,("C2", "2", "0", Capacitor(Cj)))
            circuit = Circuit(circuit)
            defs = Dict(Rl=>50.0, Ll=>300e-12, Lj=>800e-12, Cj=>1200e-15)
            wp = (2*pi*5.2e9,)
            sources = [(mode=(0,),port=1,current=0.1e-6),
                (mode=(1,),port=1,current=0.5e-6)]
            ws = 2*pi*[4.9e9]
            names = ["C2","Lj1"]
            syms = Dict("C2"=>Cj,"Lj1"=>Lj)
            solve(d, m) = hbsolve(ws, wp, sources, (4,), (8,), circuit, d;
                dc=true, keyedarrays=false, atol=1e-13,
                sensitivitynames=names, returnSsensitivity=true,
                sensitivityoperatingpoint=true, sensitivitymode=m)
            fwd = solve(defs, :forward)
            # the pump grid really does contain the self-conjugate dc mode
            @test (0,) in fwd.nonlinear.modes
            rev = solve(defs, :reverse)
            @test isapprox(fwd.linearized.Ssensitivity,
                rev.linearized.Ssensitivity, rtol = 1e-10,
                norm = v->maximum(abs,v))
            h = 1e-6
            for (k, name) in enumerate(names)
                dp = copy(defs); dp[syms[name]] *= (1+h)
                dm = copy(defs); dm[syms[name]] *= (1-h)
                Sp = hbsolve(ws, wp, sources, (4,), (8,), circuit, dp;
                    dc=true, keyedarrays=false, atol=1e-13).linearized.S
                Sm = hbsolve(ws, wp, sources, (4,), (8,), circuit, dm;
                    dc=true, keyedarrays=false, atol=1e-13).linearized.S
                fd = (Sp .- Sm)./(2*h)
                for sol in (fwd, rev)
                    @test isapprox(sol.linearized.Ssensitivity[:,:,k,1],
                        fd[:,:,1], rtol = 1e-6, norm = v->maximum(abs,v))
                end
            end
        end

        @testset "output combination invariance" begin
            # The sensitivity scaling reads the input waves of the scattering
            # parameter calculation, which later output calculations may
            # refill or overwrite, so the sensitivities must be computed
            # directly after those waves are formed and be independent of
            # which other outputs are requested, including when S itself is
            # not returned.
            JosephsonCircuits.@params Rl Cc Lj Cj
            circuit = Any[]
            push!(circuit,("P1", "1", "0", Port(1; Z0 = Rl)))
            push!(circuit,("C1", "1", "2", Capacitor(Cc))); push!(circuit,("Lj1", "2", "0", JosephsonJunction(Lj)))
            push!(circuit,("C2", "2", "0", Capacitor(Cj)))
            circuit = Circuit(circuit)
            defs = Dict(Rl=>50.0, Cc=>100e-15, Lj=>1000e-12, Cj=>1000e-15)
            wp = (2*pi*4.75001e9,)
            sources = [(mode=(1,),port=1,current=0.00565e-6)]
            ws = 2*pi*[4.5e9, 4.75e9]
            names = ["C1","Lj1"]
            solve(; kwargs...) = hbsolve(ws, wp, sources, (4,), (8,),
                circuit, defs; keyedarrays=false, sensitivitynames=names,
                returnSsensitivity=true, kwargs...).linearized.Ssensitivity
            base = solve()
            @test base == solve(returnSnoise=true)
            @test base == solve(returnnodefluxadjoint=true,
                returnvoltageadjoint=true)
            @test base == solve(returnQE=false, returnCM=false)
            @test base == solve(returnS=false, returnQE=false,
                returnCM=false, returnSnoise=false)
        end

        @testset "sensitivity mode validation" begin
            # an unknown contraction order is rejected even when no operating
            # point derivatives are in play
            JosephsonCircuits.@params Rl Cc Lj Cj
            circuit = Any[]
            push!(circuit,("P1", "1", "0", Port(1; Z0 = Rl)))
            push!(circuit,("C1", "1", "2", Capacitor(Cc))); push!(circuit,("Lj1", "2", "0", JosephsonJunction(Lj)))
            push!(circuit,("C2", "2", "0", Capacitor(Cj)))
            circuit = Circuit(circuit)
            defs = Dict(Rl=>50.0, Cc=>100e-15, Lj=>1000e-12, Cj=>1000e-15)
            @test_throws ArgumentError hblinsolve(2*pi*[4.5e9], circuit,
                defs; keyedarrays=false, sensitivitynames=["C1"],
                returnSsensitivity=true, sensitivitymode=:sideways)
        end

        @testset "unsupported components" begin
            JosephsonCircuits.@params Rl Cc Lj Cj Lla Llb Kab
            circuit = Any[]
            push!(circuit,("P1", "1", "0", Port(1; Z0 = Rl)))
            push!(circuit,("C1", "1", "2", Capacitor(Cc))); push!(circuit,("Lj1", "2", "0", JosephsonJunction(Lj)))
            push!(circuit,("C2", "2", "0", Capacitor(Cj)))
            push!(circuit,("L1", "2", "0", Inductor(Lla))); push!(circuit,("L2", "3", "0", Inductor(Llb)))
            push!(circuit,("P2", "3", "0", Port(2; Z0 = Rl)))
            push!(circuit,("K1", "L1", "L2", MutualInductor(Kab)))
            circuit = Circuit(circuit)
            defs = Dict(Rl=>50.0, Cc=>100e-15, Lj=>500e-12, Cj=>1000e-15,
                Lla=>300e-12, Llb=>300e-12, Kab=>0.5)
            ws = 2*pi*[5.0e9]
            # a mutually coupled inductor is promoted to an auxiliary branch
            # current and is not supported
            @test_throws ArgumentError hblinsolve(ws, circuit, defs;
                keyedarrays=false, sensitivitynames=["L1"],
                returnSsensitivity=true)
            # neither is a port
            @test_throws ArgumentError hblinsolve(ws, circuit, defs;
                keyedarrays=false, sensitivitynames=["P1"],
                returnSsensitivity=true)
        end
    end

    @testset verbose=true "block factorization of the linearized solve" begin
        # the dense node-block direct solve against KLU: every output of
        # the linearized sweep, on a chain with two tones and a direct
        # current mode (the modulation harmonics coupling every mode), with
        # exact factors, with single precision factors refined against the
        # double residual, and across host batches
        circuit, defs = testchaincircuit()
        w1 = 2*pi*5.0e9; w2 = 2*pi*1.19e9
        src = [(mode=(1,0), port=1, current=1.0e-6),
               (mode=(0,1), port=1, current=0.5e-6)]
        ws = 2*pi*collect(range(4.41e9, 5.57e9, length = 4))
        kw = (; dc = true, threewavemixing = true, fourwavemixing = true,
            returnSnoise = true, returnnodeflux = true, keyedarrays = false)
        ra = hbsolve(ws, (w1,w2), src, (2,2), (8,4), circuit, defs; kw...)
        for (f, tol, nb) in ((BlockFactorization(), 1e-10, 1),
                (BlockFactorization(), 1e-10, 3),
                (BlockFactorization(precision = Float32), 1e-8, 1))
            rb = hbsolve(ws, (w1,w2), src, (2,2), (8,4), circuit, defs;
                factorization = f, nbatches = nb, kw...)
            for name in (:S, :Snoise, :QE, :CM, :nodeflux)
                a = getfield(ra.linearized, name); b = getfield(rb.linearized, name)
                @test isapprox(a, b; rtol = tol)
            end
        end
        # batches of a block factorization run with one BLAS thread, which
        # is a setting of the whole process: two sweeps at once leave it
        # as they found it
        let threads = BLAS.get_num_threads()
            jpa, jpadefs = testjpacircuit()
            BLAS.set_num_threads(2)
            try
                both = [Threads.@spawn hbsolve(2*pi*[4.5e9, 4.7e9, 4.9e9],
                    (2*pi*4.75001e9,), [(mode = (1,), port = 1,
                    current = 0.00565e-6)], (2,), (8,), jpa, jpadefs;
                    factorization = BlockFactorization(), nbatches = 3,
                    keyedarrays = false) for _ in 1:2]
                foreach(fetch, both)
                @test BLAS.get_num_threads() == 2
            finally
                BLAS.set_num_threads(threads)
            end
        end
        # a promoted port resistor and a mutual inductor: the modified nodal
        # analysis rows join the node blocks, and the sensitivity path takes
        # its own factorization of the pump Jacobian
        JosephsonCircuits.@params Rl Rr Cc Lj Cj Lla Llb Kab
        circuitm = Any[]
        push!(circuitm,("P1", "1", "0", Port(1; Z0 = Rl)))
        push!(circuitm,("C1", "1", "2", Capacitor(Cc))); push!(circuitm,("Lj1", "2", "0", JosephsonJunction(Lj)))
        push!(circuitm,("C2", "2", "0", Capacitor(Cj)))
        push!(circuitm,("L1", "2", "0", Inductor(Lla))); push!(circuitm,("L2", "3", "0", Inductor(Llb)))
        push!(circuitm,("P2", "3", "0", Port(2; Z0 = Rr)))
        push!(circuitm,("K1", "L1", "L2", MutualInductor(Kab)))
        circuitm = Circuit(circuitm)
        defsm = Dict(Rl=>50.0, Rr=>50.0, Cc=>100e-15, Lj=>500e-12,
            Cj=>1000e-15, Lla=>300e-12, Llb=>300e-12, Kab=>0.5)
        wp = (2*pi*4.75001e9,)
        sources = [(mode=(1,),port=1,current=1.0e-6)]
        wsm = 2*pi*[4.5e9, 4.7e9]
        sa = hbsolve(wsm, wp, sources, (4,), (8,), circuitm, defsm;
            keyedarrays = false, sensitivitynames = ["P1/termination"],
            returnSsensitivity = true, returnSnoise = true)
        sb = hbsolve(wsm, wp, sources, (4,), (8,), circuitm, defsm;
            keyedarrays = false, sensitivitynames = ["P1/termination"],
            returnSsensitivity = true, returnSnoise = true,
            factorization = BlockFactorization())
        for name in (:S, :Snoise, :QE, :CM, :Ssensitivity)
            @test isapprox(getfield(sa.linearized, name),
                getfield(sb.linearized, name); rtol = 1e-9)
        end
        # single precision solutions: the whole solve in single, no
        # refinement, single precision agreement with the double solve
        rs = hbsolve(ws, (w1,w2), src, (2,2), (8,4), circuit, defs;
            factorization = BlockFactorization(precision = Float32,
                refine = false), kw...)
        # the accuracy of a single precision block elimination on this
        # chain: 3e-4 of the norm of S, 1e-4 of the quantum efficiency
        for (name, tol) in ((:S, 2e-3), (:Snoise, 1e-3), (:QE, 1e-3), (:CM, 1e-4))
            @test isapprox(getfield(ra.linearized, name),
                getfield(rs.linearized, name); rtol = tol)
        end
        @test !isapprox(ra.linearized.S, rs.linearized.S; rtol = 1e-12)
        # the refined single precision factorization is double accurate
        rr = hbsolve(ws, (w1,w2), src, (2,2), (8,4), circuit, defs;
            factorization = BlockFactorization(precision = Float32), kw...)
        @test isapprox(ra.linearized.S, rr.linearized.S; rtol = 1e-8)
        @test_throws MethodError hbsolve(ws, (w1,w2), src, (2,2), (8,4),
            circuit, defs; precision = Float32, kw...)
        # the automatic choice: the sparse factorization for one tone, the
        # block factorization for two or more when it fits, and the sweep
        # through it agrees with the explicit choices
        lf = JosephsonCircuits.linearizedfactorization
        d2 = hbsolve(ws[1:1], (w1,w2), src, (2,2), (8,4), circuit, defs;
            kw...)
        @test lf(sprand(ComplexF64, 20, 20, 0.3) + I, 5, 1, JosephsonCircuits.CPU()) isa KLUfactorization
        Asp = sparse(ComplexF64[1 1 0 0; 1 1 1 0; 0 1 1 1; 0 0 1 1])
        @test lf(Asp, 2, 2, JosephsonCircuits.CPU()) isa BlockFactorization
        @test lf(Asp, 2, 2, JosephsonCircuits.CPU(); budget = 0) isa KLUfactorization
        @test lf(Asp, 2, 2, JosephsonCircuits.CPU(); nbatches = 4,
            budget = 4*JosephsonCircuits.blocksystembytes(ComplexF64,
                JosephsonCircuits.clustersymbolic(
                    JosephsonCircuits.blocknodegraph(Asp, 2)...,
                    JosephsonCircuits.klunodeorder(
                        JosephsonCircuits.blocknodegraph(Asp, 2)[2])))) isa BlockFactorization
        auto = hbsolve(ws, (w1,w2), src, (2,2), (8,4), circuit, defs; kw...)
        for name in (:S, :Snoise, :QE, :CM)
            @test isapprox(getfield(ra.linearized, name),
                getfield(auto.linearized, name); rtol = 1e-9)
        end
        # the choice hblinsolve itself makes, through the sweep: the block
        # factorization for two tones and the sparse one for one tone
        function resolved(wp, Npump, Nmod, sources)
            freq = JosephsonCircuits.removeconjfreqs(JosephsonCircuits.truncfreqs(
                JosephsonCircuits.calcfreqsrdft(map(i -> 2i, Npump)); dc = true,
                odd = true, even = true, maxharmonics = Npump, w = wp))
            indices = JosephsonCircuits.fourierindices(freq)
            psc = JosephsonCircuits.compile(circuit; sorting = :number)
            nm = JosephsonCircuits.numericmatrices(psc, defs;
                Nmodes = length(freq.modes))
            nl = JosephsonCircuits.hbnlsolve(wp, sources, freq, indices, psc, nm;
                keyedarrays = false)
            sf = JosephsonCircuits.truncfreqs(JosephsonCircuits.calcfreqsdft(Nmod);
                dc = true, odd = true, even = true, maxharmonics = Nmod)
            return JosephsonCircuits.hblinsolve(ws[1:1], psc, nm.vvn, sf;
                nonlinear = nl, debuglsys = true).factorization
        end
        @test resolved((w1,w2), (8,4), (2,2), src) isa BlockFactorization
        @test resolved((w1,), (8,), (4,), [(mode=(1,), port=1, current=1.0e-6)]) isa KLUfactorization
    end

    @testset "outputs do not depend on whether S is retained" begin
        # every output consumes the per frequency view of S, so whether the
        # scattering cube is retained changes none of them
        JosephsonCircuits.@params R1v R2v C1v L1v C2v Ljv
        circuit = Any[]
        push!(circuit,("P1", "1", "0", Port(1; Z0 = R1v)))
        push!(circuit,("C1", "1", "2", Capacitor(C1v))); push!(circuit,("Lj1", "2", "0", JosephsonJunction(Ljv)))
        push!(circuit,("C2", "2", "0", Capacitor(C2v))); push!(circuit,("P2", "2", "0", Port(2; Z0 = R2v)))
        circuit = Circuit(circuit)
        defs = Dict(R1v=>50.0, R2v=>50.0, C1v=>100e-15, Ljv=>1000e-12,
            C2v=>200e-15)
        ws = 2*pi*[4.5e9, 5.0e9]
        wp = (2*pi*4.75001*1e9,)
        sources = [(mode=(1,),port=1,current=0.00565e-6)]

        withS = hbsolve(ws, wp, sources, (4,), (4,), circuit, defs;
            keyedarrays=false, returnS=true, returnQE=true, returnCM=true)
        withoutS = hbsolve(ws, wp, sources, (4,), (4,), circuit, defs;
            keyedarrays=false, returnS=false, returnQE=true, returnCM=true)

        @test isempty(withoutS.linearized.S)
        @test !isempty(withoutS.linearized.QE)
        @test withoutS.linearized.QE == withS.linearized.QE
        @test withoutS.linearized.QEideal == withS.linearized.QEideal
        @test withoutS.linearized.CM == withS.linearized.CM

        # the sensitivity scaling reads the per frequency scattering matrix,
        # so it must still be computed when only the sensitivities are asked
        # for and S itself is not returned
        names = ["C1","C2","P1/termination"]
        sensS = hblinsolve(ws, circuit, defs; keyedarrays=false,
            sensitivitynames=names, returnSsensitivity=true, returnS=true)
        sensnoS = hblinsolve(ws, circuit, defs; keyedarrays=false,
            sensitivitynames=names, returnSsensitivity=true, returnS=false)
        @test isempty(sensnoS.S)
        @test !isempty(sensnoS.Ssensitivity)
        @test sensnoS.Ssensitivity == sensS.Ssensitivity
    end

end

@testset "inputs refused by name" begin
    circuit, circuitdefs = testjpacircuit()
    ws = 2*pi*[4.5e9, 5.0e9]
    wp = (2*pi*4.75001e9,)
    sources = [(mode = (1,), port = 1, current = 0.00565e-6)]
    # a modulation harmonic count per pump tone, checked before the pump
    # is solved
    @test_throws ArgumentError hbsolve(ws, wp, sources, (2, 2), (4,), circuit,
        circuitdefs)
    # a temperature below zero
    @test_throws ArgumentError hblinsolve(ws, circuit, circuitdefs;
        temperature = -0.05)
    # an operating point of another circuit, with fewer junctions
    nl = hbnlsolve(wp, (4,), sources, circuit, circuitdefs)
    other = Circuit([(:P1, 1, 0, Port(1; Z0 = 50.0)), (:C1, 1, 2, Capacitor(100e-15)),
        (:Lj1, 2, 0, JosephsonJunction(1000e-12)), (:Lj2, 2, 3, JosephsonJunction(1000e-12)),
        (:C2, 3, 0, Capacitor(1000e-15))])
    @test_throws ArgumentError hblinsolve(ws, other; nonlinear = nl,
        Nmodulationharmonics = (2,))
    # a modulation harmonic count per pump tone, against an operating point
    @test_throws ArgumentError hblinsolve(ws, circuit, circuitdefs;
        nonlinear = nl, Nmodulationharmonics = (2, 2))
    # a complex junction inductance, under every method
    lossyjunction = Circuit([(:P1, 1, 0, Port(1; Z0 = 50.0)), (:C1, 1, 2, Capacitor(100e-15)),
        (:Lj1, 2, 0, JosephsonJunction(1000e-12 + 1e-12im)), (:C2, 2, 0, Capacitor(1000e-15))])
    for method in (NewtonKrylov(), Newton(), QuasiNewton())
        @test_throws ArgumentError hbnlsolve(wp, (4,), sources, lossyjunction;
            method = method)
    end
end

@testset "the frequency window of the pump modes" begin
    JC = JosephsonCircuits
    circuit = Any[]
    push!(circuit, ("P1", "1", "0", Port(1; Z0 = :R)))
    for i in 1:6
        push!(circuit, ("Lj$(i)", "$(i)", "$(i+1)", JosephsonJunction(:Lj)))
        push!(circuit, ("C$(i)", "$(i)", "0", Capacitor(:Cg)))
    end
    push!(circuit, ("C7", "7", "0", Capacitor(:Cg))); push!(circuit, ("R2", "7", "0", Resistor(:R)))
    circuit = Circuit(circuit)
    defs = Dict{Symbol,Complex{Float64}}(:Lj => 100e-12, :Cg => 40e-15, :R => 50.0)
    # the same circuit in the typed form, for the cache
    typed = Circuit(Any[(:P1, 1, 0, Port(1; Z0 = :R)),
        [(Symbol(:Lj, i), i, i + 1, JosephsonJunction(:Lj)) for i in 1:6]...,
        [(Symbol(:C, i), i, 0, Capacitor(:Cg)) for i in 1:6]...,
        (:C7, 7, 0, Capacitor(:Cg)), (:R2, 7, 0, Resistor(:R))])
    w = (2*pi*5.0e9, 2*pi*1.19e9)
    src = [(mode=(1,0),port=1,current=0.6e-6), (mode=(0,1),port=1,current=0.6e-6)]
    full = JC.hbnlsolve(w, (8,4), src, circuit, defs; dc = true, odd = true, even = true,
        method = Newton(), keyedarrays = false)
    @test full.solverinfo.converged
    # a floor below every retained frequency changes nothing
    same = JC.hbnlsolve(w, (8,4), src, circuit, defs; dc = true, odd = true, even = true,
        method = Newton(), keyedarrays = false, frequencywindow = (2*pi*0.1e9, Inf))
    @test same.modes == full.modes
    @test same.nodeflux == full.nodeflux
    # a box removes modes and leaves the strong ones within the size of what it dropped
    box = JC.hbnlsolve(w, (8,4), src, circuit, defs; dc = true, odd = true, even = true,
        method = Newton(), keyedarrays = false, frequencywindow = (2*pi*0.5e9, 2*pi*30e9))
    @test box.solverinfo.converged
    @test length(box.modes) < length(full.modes)
    @test all(m -> all(==(0), m) || 2*pi*0.5e9 <= abs(sum(w .* m)) <= 2*pi*30e9, box.modes)
    # a window that drops a source mode is refused with a message naming it
    @test_throws ArgumentError JC.hbnlsolve(w, (8,4), src, circuit, defs; dc = true, odd = true,
        even = true, method = Newton(), keyedarrays = false, frequencywindow = (2*pi*2e9, Inf))
    fullflux = reshape(full.nodeflux, length(full.modes), :)
    boxflux = reshape(box.nodeflux, length(box.modes), :)
    dropped = maximum(norm(fullflux[k, :]) for k in eachindex(full.modes) if !(full.modes[k] in box.modes))
    strong = maximum(norm(fullflux[k, :]) for k in eachindex(full.modes))
    for (kb, m) in enumerate(box.modes)
        kf = findfirst(==(m), full.modes)
        @test norm(boxflux[kb, :] - fullflux[kf, :]) <= 10*dropped + 1e-9*strong
    end
    # the window travels through hbsolve, hbcache and the staged solver
    hs = JC.hbsolve(2*pi*5.1e9, w, src, (1,1), (8,4), circuit, defs; dc = true,
        threewavemixing = true, fourwavemixing = true, keyedarrays = false,
        frequencywindow = (2*pi*0.5e9, 2*pi*30e9))
    @test hs.nonlinear.modes == box.modes
    cache = JC.hbcache(w, (8,4), src, typed, defs; dc = true, odd = true, even = true,
        frequencywindow = (2*pi*0.5e9, 2*pi*30e9), method = Newton())
    @test length(cache.frequencies.modes) == length(box.modes)
    st = JC.hbnlsolve(w, (8,4), src, circuit, defs; dc = true, odd = true, even = true,
        method = Staged(), keyedarrays = false, frequencywindow = (2*pi*0.5e9, 2*pi*30e9))
    @test st.solverinfo.converged
    @test st.modes == box.modes
end

@testset "the evaluation grid of the pump modes" begin
    JC = JosephsonCircuits
    circuit = Any[]
    push!(circuit, ("P1", "1", "0", Port(1; Z0 = :R)))
    for i in 1:6
        push!(circuit, ("Lj$(i)", "$(i)", "$(i+1)", JosephsonJunction(:Lj)))
        push!(circuit, ("C$(i)", "$(i)", "0", Capacitor(:Cg)))
    end
    push!(circuit, ("C7", "7", "0", Capacitor(:Cg))); push!(circuit, ("R2", "7", "0", Resistor(:R)))
    circuit = Circuit(circuit)
    defs = Dict{Symbol,Complex{Float64}}(:Lj => 100e-12, :Cg => 40e-15, :R => 50.0)
    # the same circuit in the typed form, for the cache
    typed = Circuit(Any[(:P1, 1, 0, Port(1; Z0 = :R)),
        [(Symbol(:Lj, i), i, i + 1, JosephsonJunction(:Lj)) for i in 1:6]...,
        [(Symbol(:C, i), i, 0, Capacitor(:Cg)) for i in 1:6]...,
        (:C7, 7, 0, Capacitor(:Cg)), (:R2, 7, 0, Resistor(:R))])
    w = (2*pi*5.0e9, 2*pi*1.19e9)
    src = [(mode=(1,0),port=1,current=0.6e-6), (mode=(0,1),port=1,current=0.6e-6)]
    kw = (; dc = true, odd = true, even = true, method = Newton(), keyedarrays = false)
    # the default grid is twice the retained set, which is Nharmonics
    sol = JC.hbnlsolve(w, (8,4), src, circuit, defs; kw...)
    @test sol.solverinfo.converged
    @test sol.frequencies.Nharmonics == (16,8)
    @test all(m -> all(abs.(m) .<= (8,4)), sol.modes)
    # the retained set does not depend on the grid
    native = JC.hbnlsolve(w, (8,4), src, circuit, defs; kw..., Nevaluationharmonics = (8,4))
    wide = JC.hbnlsolve(w, (8,4), src, circuit, defs; kw..., Nevaluationharmonics = (24,12))
    @test native.frequencies.Nharmonics == (8,4)
    @test wide.frequencies.Nharmonics == (24,12)
    @test native.modes == sol.modes == wide.modes
    # the dealiased default is closer to the wide grid than the native grid is
    @test norm(sol.nodeflux - wide.nodeflux) < norm(native.nodeflux - wide.nodeflux)
    @test norm(sol.nodeflux - wide.nodeflux) < 1e-3*norm(wide.nodeflux)
    # a grid smaller than the retained set is refused everywhere
    @test_throws ArgumentError JC.hbnlsolve(w, (8,4), src, circuit, defs; kw..., Nevaluationharmonics = (8,3))
    @test_throws ArgumentError JC.hbnlsolve(w, (8,4), src, circuit, defs; kw..., method = Staged(), Nevaluationharmonics = (8,3))
    @test_throws ArgumentError JC.hbsolve(2*pi*5.1e9, w, src, (1,1), (8,4), circuit, defs; dc = true,
        threewavemixing = true, fourwavemixing = true, keyedarrays = false, Nevaluationharmonics = (7,4))
    @test_throws ArgumentError JC.hbcache(w, (8,4), src, typed, defs;
        dc = true, odd = true, even = true, Nevaluationharmonics = (8,3))
    # the grid travels through hbsolve, hbcache and the staged solver
    hs = JC.hbsolve(2*pi*5.1e9, w, src, (1,1), (8,4), circuit, defs; dc = true,
        threewavemixing = true, fourwavemixing = true, keyedarrays = false, Nevaluationharmonics = (24,12))
    @test hs.nonlinear.frequencies.Nharmonics == (24,12)
    @test hs.nonlinear.modes == wide.modes
    @test isapprox(hs.nonlinear.nodeflux, wide.nodeflux; rtol = 1e-6)
    hsd = JC.hbsolve(2*pi*5.1e9, w, src, (1,1), (8,4), circuit, defs; dc = true,
        threewavemixing = true, fourwavemixing = true, keyedarrays = false)
    @test hsd.nonlinear.frequencies.Nharmonics == (16,8)
    st = JC.hbnlsolve(w, (8,4), src, circuit, defs; kw..., method = Staged(), Nevaluationharmonics = (24,12))
    @test st.frequencies.Nharmonics == (24,12)
    @test isapprox(st.nodeflux, wide.nodeflux; rtol = 1e-6)
    cache = JC.hbcache(w, (8,4), src, typed, defs;
        dc = true, odd = true, even = true, Nevaluationharmonics = (24,12), keyedarrays = false)
    @test cache.frequencies.Nharmonics == (24,12)
    cached = JC.hbsolve!(cache, (Lj = 100e-12, Cg = 40e-15, R = 50.0))
    @test isapprox(cached.nodeflux, wide.nodeflux; rtol = 1e-6)
end
