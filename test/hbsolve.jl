using JosephsonCircuits
using LinearAlgebra
using SparseArrays
using Test

@testset verbose=true "hbsolve" begin


    @testset "hbsolve-hbnlsolve comparison" begin

        ftol = 5e-17

        JosephsonCircuits.@params R Cc Lj Cj
        circuit = [
            ("P1","1","0",1),
            ("R1","1","0",R),
            ("C1","1","2",Cc),
            ("Lj1","2","0",Lj),
            ("C2","2","0",Cj)]

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
            Npumpharmonics = (10,)
            Nmodulationharmonics = (10,)
            sol1 = hbsolve(ws, (wp,), sources, Nmodulationharmonics,
                Npumpharmonics, circuit, circuitdefs, ftol = ftol)
            S1ss = sol1.linearized.S((0,),1,(0,),1,1)
            S1is = sol1.linearized.S((-2,),1,(0,),1,1)

            # nonlinear simulation with (pump,signal) order
            w = (wp,ws)
            Nharmonics = (10,10)
            sources = [(mode=(1,0),port=1,current=Ip),(mode=(0,1),port=1,current=Is)]
            sol2 = hbnlsolve(w, Nharmonics, sources, circuit, circuitdefs, ftol = ftol)
            S2ss = sol2.S((0,1),1,(0,1),1)
            S2is = sol2.S((2,-1),1,(0,1),1)

            # nonlinear simulation with (signal,pump) order
            w = (ws,wp)
            Nharmonics = (10,10)
            sources = [(mode=(0,1),port=1,current=Ip),(mode=(1,0),port=1,current=Is)]
            sol3 = hbnlsolve(w, Nharmonics, sources, circuit, circuitdefs, ftol = ftol)
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

    @testset "hbsolve method comparison" begin

        JosephsonCircuits.@params Rleft Cc Lj Cj w L1
        circuit = Tuple{String,String,String,Any}[]
        push!(circuit,("P1","1","0",1))
        push!(circuit,("R1","1","0",Rleft))
        push!(circuit,("C1","1","2",Cc)) 
        push!(circuit,("Lj1","2","0",Lj)) 
        push!(circuit,("C2","2","0",Cj))
        circuitdefs = Dict(
            Lj =>1000.0e-12,
            Cc => 100.0e-15,
            Cj => 1000.0e-15,
            Rleft => 50.0,
        )
        ws = 2*pi*(4.5:0.01:5.0)*1e9
        wp = (2*pi*4.75001*1e9,)
        sources = [(mode=(1,),port=1,current=0.00565e-6)]
        Nmodulationharmonics = (8,)
        Npumpharmonics = (16,)

        quasinewton = hbsolve(ws, wp, sources, Nmodulationharmonics,
            Npumpharmonics, circuit, circuitdefs, ftol=1e-12, method = :quasinewton)


        @test isapprox(
            quasinewton.nonlinear.nodeflux[:],
            ComplexF64[-0.013189575486618105 - 0.00865077163136891im, 2.6396823809835998e-5 - 5.667772212498911e-6im, -2.3501629748681806e-8 - 2.2306559257777402e-8im, -7.987057959121183e-12 + 3.750064665528795e-11im, 4.337129743740559e-14 - 1.3397330127340191e-14im, -3.812185723603184e-17 - 3.8456119951452504e-17im, -4.278102807441037e-20 + 3.974761087481868e-20im, -8.24749227891514e-25 + 3.0560138516751506e-24im, 0.12157593753208869 + 0.07973582696437984im, 1.3736443951855489e-5 - 6.463164795537436e-5im, -5.3397980865816756e-8 + 9.191484285021433e-9im, 2.7913096743299775e-11 + 4.514682457284182e-11im, 3.3395873694437265e-14 - 4.569085921919452e-14im, -6.154948861291708e-17 - 1.523212283849711e-17im, -2.2291909276375226e-20 + 6.18003971235134e-20im, 5.405242786243393e-25 + 3.424470385254823e-24im],
            atol = 1e-8)

        @test isapprox(
            10*log10.(abs2.(quasinewton.linearized.S((0,),1,(0,),1,:))),
            [0.002717249265375112, 0.0031881165263016476, 0.003764438371177404, 0.004475612824820614, 0.005360972358329889, 0.006473775210238347, 0.007887075336301673, 0.009702494973065136, 0.012063564112722006, 0.0151763872044898, 0.019342315628701326, 0.025010742306207902, 0.03286645328833334, 0.043977897966981386, 0.060055891091184574, 0.08391851842302066, 0.12035318288996377, 0.17776810664422604, 0.27146258097536335, 0.43031121414540996, 0.7108074673012084, 1.2271520892520107, 2.2163506977088527, 4.180249149477177, 8.16945893982768, 13.302337843265873, 8.180459085825893, 4.185698117374451, 2.2190409577860577, 1.2285225612890744, 0.7115311472717254, 0.4307082670317731, 0.2716887696865355, 0.17790149553267967, 0.12043426616790562, 0.08396907049655332, 0.06008804086742857, 0.04399863098219145, 0.0328799192406124, 0.02501947707672059, 0.019347909243977855, 0.015179860051320881, 0.012065584850116115, 0.009703510073587402, 0.007887387964370132, 0.006473594325612282, 0.005360443913078622, 0.0044748399587568825, 0.003763494713109064, 0.0031870550949444913, 0.002716108513561703],
            atol = 1e-6)

        newton = hbsolve(ws, wp, sources, Nmodulationharmonics,
            Npumpharmonics, circuit, circuitdefs, ftol=1e-12, method = :newton)

        @test isapprox(
            quasinewton.nonlinear.nodeflux[:],
            newton.nonlinear.nodeflux[:],
            atol = 1e-8)

        @test isapprox(
            10*log10.(abs2.(quasinewton.linearized.S((0,),1,(0,),1,:))),
            10*log10.(abs2.(newton.linearized.S((0,),1,(0,),1,:))),
            atol = 1e-6)

        # test some uncommon options
        JosephsonCircuits.@params w
        result = hbsolve(ws, wp, sources, Nmodulationharmonics,
            Npumpharmonics, circuit, circuitdefs, ftol=1e-12, symfreqvar = w,
            returnS=false, returnSnoise=true, returnQE=false,
            returnnodeflux=true,
            returnnodefluxadjoint=true, returnCM=false,
            returnvoltage=true, returnvoltageadjoint=true,
            returnSsensitivity = true,
            sensitivitynames=["C1"],
            nbatches=4)

        @test result.linearized.QE == Array{Float64, 3}(undef, 0, 0, 0)
        @test result.linearized.S == Array{Float64, 3}(undef, 0, 0, 0)
        # @test result.linearized.nodeflux == ComplexF64[]
        @test result.linearized.Snoise[:] == ComplexF64[]
        @test result.linearized.CM == Matrix{Float64}(undef, 0, 0)

    end



    @testset "hbsolve initial nodeflux" begin

        circuit = Array{Tuple{String,String,String,Union{Complex{Float64},Symbol,Int64}},1}(undef,0)
        push!(circuit,("P1","1","0",1))
        push!(circuit,("I1","1","0",:Ipump))
        push!(circuit,("R1","1","0",:Rleft))
        push!(circuit,("L1","1","0",:Lm)) 
        push!(circuit,("K1","L1","L2",:K1))
        push!(circuit,("C1","1","2",:Cc)) 
        push!(circuit,("L2","2","3",:Lm)) 
        push!(circuit,("Lj3","3","0",:Lj)) 
        push!(circuit,("Lj4","2","0",:Lj)) 
        push!(circuit,("C2","2","0",:Cj))
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
        circuit = [
            ("P1","1","0",1),
            ("R1","1","0",R),
            ("C1","1","2",Cc),
            ("Lj1","2","0",Lj),
            ("C2","2","0",Cj)]
        
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

        # loop over all of the flags, setting one of them to be true
        for j in 1:length(flags)
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
        circuit = Tuple{String,String,String,Any}[]
        push!(circuit,("P1","1","0",1))
        push!(circuit,("R1","1","0",Rleft))
        push!(circuit,("C1","1","2",Cc)) 
        push!(circuit,("Lj1","2","0",Lj)) 
        push!(circuit,("C2","2","0",Cj))
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

        @test_warn(
            "Solver did not converge.",
            hbnlsolve(w,Nharmonics,sources,circuit,circuitdefs,iterations=1)
        )
    end

    @testset "hbnlsolve simple testcase" begin

        circuit = [("P1","1","0",1),("R1","1","0",50.0)]
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

        # a port and resistor with no inductive path to ground used to make
        # the nodal system matrix structurally singular when a DC mode was
        # included, producing a NaN in the line search under the previous
        # solver and a SingularException under the current one. with the
        # modified nodal analysis formulation the DC node flux is gauge
        # fixed and the circuit solves exactly.
        circuit = [("P1","1","0",1),("R1","1","0",50.0)]
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
        circuit = [("P1","1","0",1),("R1","1","0",Rv),("C1","1","2",Ccv),
            ("Lj1","2","0",Ljv),("C2","2","0",Cjv)]
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

        # frequency dependent values through symfreqvar are accepted, and
        # a value mixing symfreqvar with an undefined variable is rejected
        # naming both the component and the variable
        JosephsonCircuits.@params wsym Rundef
        c2 = [("P1","1","0",1),("R1","1","0",50.0 + 0.0*wsym),
            ("C1","1","0",100.0e-15),("L1","1","0",1.0e-9)]
        out = JosephsonCircuits.hbnlsolve(wp, (1,), sources, c2, Dict();
            symfreqvar = wsym)
        @test out.solverinfo.converged
        c3 = [("P1","1","0",1),("R1","1","0",Rundef/(1 + wsym*1e-12)),
            ("C1","1","0",100.0e-15),("L1","1","0",1.0e-9)]
        err3 = try
            JosephsonCircuits.hbnlsolve(wp, (1,), sources, c3, Dict();
                symfreqvar = wsym)
            nothing
        catch e
            e
        end
        @test err3 isa ArgumentError
        @test occursin("Rundef", sprint(showerror, err3))
        @test occursin("R1", sprint(showerror, err3))
    end

    @testset "FrequencyDependent matches symfreqvar" begin
        # the same frequency law spelled two ways: an expression in the
        # symbolic frequency variable, and a plain Julia closure through
        # FrequencyDependent. The scattering parameters must agree to
        # roundoff.
        wp = (2*pi*4.75001*1e9,)
        ws = 2*pi*(4.5:0.1:5.0)*1e9
        sources = [(mode=(1,),port=1,current=0.00565e-6)]
        JosephsonCircuits.@params wsym
        law(w) = 50.0*(1 + (w/1e11)^2)
        csym = [("P1","1","0",1),("R1","1","0",law(wsym)),
            ("C1","1","2",100.0e-15),("Lj1","2","0",1000.0e-12),
            ("C2","2","0",1000.0e-15)]
        cfun = [("P1","1","0",1),
            ("R1","1","0",FrequencyDependent(law)),
            ("C1","1","2",100.0e-15),("Lj1","2","0",1000.0e-12),
            ("C2","2","0",1000.0e-15)]
        Ssym = hbsolve(ws, wp, sources, (2,), (8,), csym, Dict();
            symfreqvar = wsym).linearized.S
        Sfun = hbsolve(ws, wp, sources, (2,), (8,), cfun).linearized.S
        @test isapprox(Array(Ssym), Array(Sfun), rtol = 1e-12)
    end

    @testset "calcsources errors" begin

        modes = [(0,), (1,)]
        portindices = [1]
        portnumbers = [1]
        nodeindices = [2 2 2 2 0 2 3 4 3 3; 1 1 1 1 0 3 4 1 1 1]
        edge2indexdict = Dict((1, 2) => 1, (3, 1) => 2, (1, 3) => 2, (4, 1) => 3, (2, 1) => 1, (1, 4) => 3, (3, 4) => 4, (4, 3) => 4)
        Lmean = 1.005e-9 + 0.0im
        Nnodes = 4
        Nbranches = 4
        Nmodes = 2

        # current source for non-existent port
        sources = [(mode = (0,), port = 1, current = 0.0005), (mode = (1,), port = 2, current = 1.0e-10)]
        @test_throws(
            ArgumentError("Source port 2 not found."),
            JosephsonCircuits.calcsources(modes, sources, portindices, portnumbers,
                nodeindices, edge2indexdict, Lmean, Nnodes, Nbranches, Nmodes))

        # current source for non-existent mode
        sources = [(mode = (0,), port = 1, current = 0.0005), (mode = (2,), port = 1, current = 1.0e-10)]
        @test_throws(
            ArgumentError("Source mode (2,) is not among the retained modes; the truncation (`Nharmonics`, `maxintermodorder`, `dc`, `odd`, `even`, `frequencywindow`) removed it."),
            JosephsonCircuits.calcsources(modes, sources, portindices, portnumbers,
                nodeindices, edge2indexdict, Lmean, Nnodes, Nbranches, Nmodes))

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
        circuitjpa = Tuple{String,String,String,Any}[]
        push!(circuitjpa,("P1","1","0",1))
        push!(circuitjpa,("R1","1","0",Rleft))
        push!(circuitjpa,("C1","1","2",Cc))
        push!(circuitjpa,("Lj1","2","0",Lj))
        push!(circuitjpa,("C2","2","0",Cj))
        circuitdefsjpa = Dict(Lj=>1000.0e-12, Cc=>100.0e-15, Cj=>1000.0e-15,
            Rleft=>50.0)

        # a lossy JPA: the complex capacitance adds a noise port, which is the
        # consumer of the adjoint solution we most care about here
        circuitdefsjpalossy = Dict(Lj=>1000.0e-12, Cc=>100.0e-15,
            Cj=>1000.0e-15/(1+1e-3im), Rleft=>50.0)

        # two ports and a mutually coupled inductor pair, which is promoted to
        # auxiliary branch currents as well. exercises both auxiliary blocks at
        # once, with the coupling coefficient close to one.
        circuitmutual = Tuple{String,String,String,Any}[]
        push!(circuitmutual,("P1","1","0",1))
        push!(circuitmutual,("R1","1","0",Rleft))
        push!(circuitmutual,("C1","1","2",Cc))
        push!(circuitmutual,("Lj1","2","0",Lj))
        push!(circuitmutual,("C2","2","0",Cj))
        push!(circuitmutual,("L1","2","0",Lla))
        push!(circuitmutual,("L2","3","0",Llb))
        push!(circuitmutual,("P2","3","0",2))
        push!(circuitmutual,("R2","3","0",Rright))
        push!(circuitmutual,("K1","L1","L2",Kab))
        circuitdefsmutual = Dict(Lj=>500.0e-12, Cc=>100.0e-15, Cj=>1000.0e-15,
            Lla=>300.0e-12, Llb=>300.0e-12, Rleft=>50.0, Rright=>50.0,
            Kab=>0.99)

        testcases = (
            ("single-tone JPA", (2*pi*4.75001e9,),
                [(mode=(1,),port=1,current=0.00565e-6)], (8,), (4,),
                circuitjpa, circuitdefsjpa),
            ("single-tone lossy JPA", (2*pi*4.75001e9,),
                [(mode=(1,),port=1,current=0.00565e-6)], (8,), (4,),
                circuitjpa, circuitdefsjpalossy),
            ("two-tone JPA", (2*pi*4.65001e9, 2*pi*4.85001e9),
                [(mode=(1,0),port=1,current=0.00565e-6*1.7),
                 (mode=(0,1),port=1,current=0.00565e-6*1.7)], (4,4), (2,2),
                circuitjpa, circuitdefsjpa),
            ("mutual inductor", (2*pi*4.75001e9,),
                [(mode=(1,),port=1,current=1.0e-6)], (6,), (4,),
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
                cg = JosephsonCircuits.calccircuitgraph(psc)
                signalfreq = JosephsonCircuits.truncfreqs(
                    JosephsonCircuits.calcfreqsdft(Nmodulationharmonics);
                    dc=true, odd=false, even=true, maxintermodorder=Inf)
                d = JosephsonCircuits.hblinsolve(ws, psc, cg, circuitdefs,
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
                    for (r, ci) in enumerate(d.mnaindices)
                        g = 1/d.vvn[ci]
                        for m in 1:d.Nmodes
                            D[d.Nnodalmna + (r-1)*d.Nmodes + m] =
                                im*wmodes[m]*g
                        end
                    end
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
                    # with no promoted components there are no asymmetric
                    # auxiliary rows; if promoted elements (e.g. voltage
                    # sources) are ever added, their rows differ under the
                    # similarity and this guard reactivates
                    if d.Nauxmnar > 0
                        aux = (d.Nnodalmna+1):(d.Nnodalmna+d.Nauxmnar)
                        @test !isapprox(xconj[aux,:], xtrans[aux,:])
                    end
                end

                # end to end: the transposed solve used by hblinsolve gives the
                # same adjoint node fluxes as an independent solve of the
                # conjugated pump system.
                sol = JosephsonCircuits.hblinsolve(ws, psc, cg, circuitdefs,
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

        # dS/dr, the derivative of the scattering matrix with respect to a
        # relative perturbation of a component value at a fixed pump
        # operating point, against central finite differences. The pump
        # operating point is held fixed in the finite differences to match
        # the definition, by reusing one nonlinear solution while perturbing
        # the component values of the linearized solve.

        @testset "linear network" begin
            JosephsonCircuits.@params R1v R2v R3v C1v L1v C2v
            circuit = Tuple{String,String,String,Any}[]
            push!(circuit,("P1","1","0",1)); push!(circuit,("R1","1","0",R1v))
            push!(circuit,("C1","1","2",C1v)); push!(circuit,("L1","2","0",L1v))
            push!(circuit,("C2","2","0",C2v)); push!(circuit,("P2","2","0",2))
            push!(circuit,("R2","2","0",R2v)); push!(circuit,("R3","1","2",R3v))
            defs = Dict(R1v=>50.0, R2v=>50.0, R3v=>300.0, C1v=>100e-15,
                L1v=>1e-9, C2v=>200e-15)
            ws = 2*pi*[5.0e9, 7.0e9]
            names = ["C1","L1","C2","R3","R1","R2"]
            syms = Dict("C1"=>C1v,"L1"=>L1v,"C2"=>C2v,"R3"=>R3v,
                "R1"=>R1v,"R2"=>R2v)
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
            circuit = Tuple{String,String,String,Any}[]
            push!(circuit,("P1","1","0",1)); push!(circuit,("R1","1","0",Rl))
            push!(circuit,("C1","1","2",Cc)); push!(circuit,("Lj1","2","0",Lj))
            push!(circuit,("C2","2","0",Cj))
            defs = Dict(Rl=>50.0, Cc=>100e-15, Lj=>1000e-12, Cj=>1000e-15)
            wp = (2*pi*4.75001e9,)
            sources = [(mode=(1,),port=1,current=0.00565e-6)]
            ws = 2*pi*[4.5e9, 4.75e9]
            names = ["C1","C2","Lj1","R1"]
            syms = Dict("C1"=>Cc,"C2"=>Cj,"Lj1"=>Lj,"R1"=>Rl)

            sol = hbsolve(ws, wp, sources, (8,), (16,), circuit, defs;
                keyedarrays=false, sensitivitynames=names,
                returnSsensitivity=true,sensitivityoperatingpoint=false)

            # one nonlinear solution, reused so the operating point is fixed
            nonlinear = hbnlsolve(wp, (16,), sources, circuit, defs;
                keyedarrays=false)
            psc = JosephsonCircuits.compile(circuit)
            cg = JosephsonCircuits.calccircuitgraph(psc)
            signalfreq = JosephsonCircuits.truncfreqs(
                JosephsonCircuits.calcfreqsdft((8,)); dc=true, odd=false,
                even=true, maxintermodorder=Inf)
            frozen(d) = JosephsonCircuits.hblinsolve(ws, psc, cg, d,
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
            solk = hbsolve(ws, wp, sources, (8,), (16,), circuit, defs;
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
            circuit = Tuple{String,String,String,Any}[]
            push!(circuit,("P1","1","0",1)); push!(circuit,("R1","1","0",Rl))
            push!(circuit,("C1","1","2",Cc)); push!(circuit,("Lj1","2","0",Lj))
            push!(circuit,("C2","2","0",Cj))
            defs = Dict(Rl=>50.0, Cc=>100e-15, Lj=>1000e-12, Cj=>1000e-15)
            wp = (2*pi*4.75001e9,)
            sources = [(mode=(1,),port=1,current=0.00565e-6)]
            ws = 2*pi*[4.5e9, 4.75e9]
            names = ["C1","C2","Lj1","R1"]
            syms = Dict("C1"=>Cc,"C2"=>Cj,"Lj1"=>Lj,"R1"=>Rl)
            solve(d; op=false) = hbsolve(ws, wp, sources, (8,), (16,),
                circuit, d; keyedarrays=false, ftol=1e-13,
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
            circuit = Tuple{String,String,String,Any}[]
            push!(circuit,("P1","1","0",1)); push!(circuit,("R1","1","0",Rl))
            for i in 1:Ncells
                push!(circuit,("Lj$(i)","$(i)","$(i+1)",Ljx))
                push!(circuit,("Cj$(i)","$(i)","$(i+1)",Cjx))
                push!(circuit,("Cg$(i)","$(i+1)","0", i == 1 ? Cg1v : Cg))
            end
            push!(circuit,("P2","$(Ncells+1)","0",2))
            push!(circuit,("R2","$(Ncells+1)","0",Rl))
            defs = Dict(Rl=>50.0, Ljx=>JosephsonCircuits.IctoLj(3.4e-6),
                Cg=>45e-15, Cg1v=>45e-15, Cjx=>55e-15)
            ws = 2*pi*[5.5e9, 6.4e9]
            names = ["Cg1","Cg2","Lj2","R1"]

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
                        circuit, d; keyedarrays=false, ftol=1e-13,
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
                        keyedarrays=false, ftol=1e-13).linearized.S
                    Sm = hbsolve(ws, wp, sources, Nsig, Npump, circuit, dm;
                        keyedarrays=false, ftol=1e-13).linearized.S
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
            circuit = Tuple{String,String,String,Any}[]
            push!(circuit,("P1","1","0",1)); push!(circuit,("R1","1","0",Rl))
            push!(circuit,("C1","1","2",Cc)); push!(circuit,("Lj1","2","0",Lj))
            push!(circuit,("C2","2","0",Cj))
            defs = Dict(Rl=>50.0, Cc=>100e-15, Lj=>1000e-12, Cj=>200e-15)
            wp = (2*pi*5e9,)
            sources = [(mode=(1,),port=1,current=1e-7)]
            ws = 2*pi*[4.5e9]
            psc = compile(circuit)
            cg = calccircuitgraph(psc)
            nl = hbnlsolve(wp, (4,), sources, circuit, defs;
                returnoperatingpoint=true)
            op = nl.operatingpoint
            signalfreq = JosephsonCircuits.truncfreqs(
                JosephsonCircuits.calcfreqsdft((2,)); dc=true, odd=false,
                even=true)
            names = ["C2"]
            good = JosephsonCircuits.calcresidualsensitivity(op, psc, cg,
                JosephsonCircuits.numericmatrices(psc, cg, defs,
                    Nmodes=length(nl.modes)),
                [psc.componentnamedict[n] for n in names])
            lin(; kwargs...) = hblinsolve(ws, psc, cg, defs, signalfreq;
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
            # and the (formerly junction shaped) operating point machinery
            # must not be constructed at all.
            JosephsonCircuits.@params Rl Ll Cs
            circuit = Tuple{String,String,String,Any}[]
            push!(circuit,("P1","1","0",1)); push!(circuit,("R1","1","0",Rl))
            push!(circuit,("L1","1","2",Ll)); push!(circuit,("C1","2","0",Cs))
            push!(circuit,("P2","2","0",2)); push!(circuit,("R2","2","0",Rl))
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
            circuit = Tuple{String,String,String,Any}[]
            push!(circuit,("P1","1","0",1)); push!(circuit,("R1","1","0",Rl))
            push!(circuit,("L1","1","2",Ll)); push!(circuit,("Lj1","2","0",Lj))
            push!(circuit,("C2","2","0",Cj))
            defs = Dict(Rl=>50.0, Ll=>300e-12, Lj=>800e-12, Cj=>1200e-15)
            wp = (2*pi*5.2e9,)
            sources = [(mode=(0,),port=1,current=0.1e-6),
                (mode=(1,),port=1,current=0.5e-6)]
            ws = 2*pi*[4.9e9]
            names = ["C2","Lj1"]
            syms = Dict("C2"=>Cj,"Lj1"=>Lj)
            solve(d, m) = hbsolve(ws, wp, sources, (4,), (8,), circuit, d;
                dc=true, keyedarrays=false, ftol=1e-13,
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
                    dc=true, keyedarrays=false, ftol=1e-13).linearized.S
                Sm = hbsolve(ws, wp, sources, (4,), (8,), circuit, dm;
                    dc=true, keyedarrays=false, ftol=1e-13).linearized.S
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
            circuit = Tuple{String,String,String,Any}[]
            push!(circuit,("P1","1","0",1)); push!(circuit,("R1","1","0",Rl))
            push!(circuit,("C1","1","2",Cc)); push!(circuit,("Lj1","2","0",Lj))
            push!(circuit,("C2","2","0",Cj))
            defs = Dict(Rl=>50.0, Cc=>100e-15, Lj=>1000e-12, Cj=>1000e-15)
            wp = (2*pi*4.75001e9,)
            sources = [(mode=(1,),port=1,current=0.00565e-6)]
            ws = 2*pi*[4.5e9, 4.75e9]
            names = ["C1","Lj1"]
            solve(; kwargs...) = hbsolve(ws, wp, sources, (8,), (16,),
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

        @testset "port node order" begin
            # The adjoint solve uses the source columns, which carry the
            # canonical orientation of the port branch, while the output
            # functional differences the node fluxes in the node order of the
            # port component. Writing a port with its nodes in the opposite
            # order must not flip the sign of the sensitivities.
            JosephsonCircuits.@params R1v R2v C1v L1v C2v
            for reversed in (false, true)
                circuit = Tuple{String,String,String,Any}[]
                push!(circuit,("P1","1","0",1))
                push!(circuit,("R1","1","0",R1v))
                push!(circuit,("C1","1","2",C1v))
                push!(circuit,("L1","2","0",L1v))
                push!(circuit,("C2","2","0",C2v))
                if reversed
                    push!(circuit,("P2","0","2",2))
                    push!(circuit,("R2","0","2",R2v))
                else
                    push!(circuit,("P2","2","0",2))
                    push!(circuit,("R2","2","0",R2v))
                end
                defs = Dict(R1v=>50.0, R2v=>50.0, C1v=>100e-15, L1v=>1e-9,
                    C2v=>200e-15)
                ws = 2*pi*[5.0e9]
                names = ["C1","L1"]
                syms = Dict("C1"=>C1v,"L1"=>L1v)
                sol = hblinsolve(ws, circuit, defs; keyedarrays=false,
                    sensitivitynames=names, returnSsensitivity=true)
                h = 1e-6
                for (k, name) in enumerate(names)
                    dp = copy(defs); dp[syms[name]] *= (1+h)
                    dm = copy(defs); dm[syms[name]] *= (1-h)
                    Sp = hblinsolve(ws, circuit, dp; keyedarrays=false).S
                    Sm = hblinsolve(ws, circuit, dm; keyedarrays=false).S
                    fd = (Sp .- Sm)./(2*h)
                    @test isapprox(sol.Ssensitivity[:,:,k,1], fd[:,:,1],
                        rtol = 1e-6, norm = v->maximum(abs,v))
                end
            end
        end

        @testset "mutual inductor with promoted port resistor" begin
            # A promoted port resistor sensitivity with the operating point
            # shift, in a circuit which also contains a mutually coupled
            # inductor pair, so the operating point augmentation includes the
            # coupled inductor auxiliary variables as well.
            JosephsonCircuits.@params Rl Rr Cc Lj Cj Lla Llb Kab
            circuit = Tuple{String,String,String,Any}[]
            push!(circuit,("P1","1","0",1)); push!(circuit,("R1","1","0",Rl))
            push!(circuit,("C1","1","2",Cc)); push!(circuit,("Lj1","2","0",Lj))
            push!(circuit,("C2","2","0",Cj))
            push!(circuit,("L1","2","0",Lla)); push!(circuit,("L2","3","0",Llb))
            push!(circuit,("P2","3","0",2)); push!(circuit,("R2","3","0",Rr))
            push!(circuit,("K1","L1","L2",Kab))
            defs = Dict(Rl=>50.0, Rr=>50.0, Cc=>100e-15, Lj=>500e-12,
                Cj=>1000e-15, Lla=>300e-12, Llb=>300e-12, Kab=>0.5)
            wp = (2*pi*4.75001e9,)
            sources = [(mode=(1,),port=1,current=1.0e-6)]
            ws = 2*pi*[4.5e9]
            sol = hbsolve(ws, wp, sources, (4,), (8,), circuit, defs;
                keyedarrays=false, ftol=1e-13, sensitivitynames=["R1"],
                returnSsensitivity=true, sensitivityoperatingpoint=true)
            h = 1e-6
            dp = copy(defs); dp[Rl] *= (1+h)
            dm = copy(defs); dm[Rl] *= (1-h)
            Sp = hbsolve(ws, wp, sources, (4,), (8,), circuit, dp;
                keyedarrays=false, ftol=1e-13).linearized.S
            Sm = hbsolve(ws, wp, sources, (4,), (8,), circuit, dm;
                keyedarrays=false, ftol=1e-13).linearized.S
            fd = (Sp .- Sm)./(2*h)
            @test isapprox(sol.linearized.Ssensitivity[:,:,1,1], fd[:,:,1],
                rtol = 1e-4, norm = v->maximum(abs,v))
        end

        @testset "sensitivity mode validation" begin
            # an unknown contraction order is rejected even when no operating
            # point derivatives are in play
            JosephsonCircuits.@params Rl Cc Lj Cj
            circuit = Tuple{String,String,String,Any}[]
            push!(circuit,("P1","1","0",1)); push!(circuit,("R1","1","0",Rl))
            push!(circuit,("C1","1","2",Cc)); push!(circuit,("Lj1","2","0",Lj))
            push!(circuit,("C2","2","0",Cj))
            defs = Dict(Rl=>50.0, Cc=>100e-15, Lj=>1000e-12, Cj=>1000e-15)
            @test_throws ArgumentError hblinsolve(2*pi*[4.5e9], circuit,
                defs; keyedarrays=false, sensitivitynames=["C1"],
                returnSsensitivity=true, sensitivitymode=:sideways)
        end

        @testset "unsupported components" begin
            JosephsonCircuits.@params Rl Cc Lj Cj Lla Llb Kab
            circuit = Tuple{String,String,String,Any}[]
            push!(circuit,("P1","1","0",1)); push!(circuit,("R1","1","0",Rl))
            push!(circuit,("C1","1","2",Cc)); push!(circuit,("Lj1","2","0",Lj))
            push!(circuit,("C2","2","0",Cj))
            push!(circuit,("L1","2","0",Lla)); push!(circuit,("L2","3","0",Llb))
            push!(circuit,("P2","3","0",2)); push!(circuit,("R2","3","0",Rl))
            push!(circuit,("K1","L1","L2",Kab))
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

    @testset "outputs do not depend on whether S is retained" begin
        # requestS used to be forced true by returnQE and returnSsensitivity,
        # so the scattering cube survived the frequency loop even when the
        # caller never asked for it. Everything now consumes the per frequency
        # view instead, and these outputs must be unchanged either way.
        JosephsonCircuits.@params R1v R2v C1v L1v C2v Ljv
        circuit = Tuple{String,String,String,Any}[]
        push!(circuit,("P1","1","0",1)); push!(circuit,("R1","1","0",R1v))
        push!(circuit,("C1","1","2",C1v)); push!(circuit,("Lj1","2","0",Ljv))
        push!(circuit,("C2","2","0",C2v)); push!(circuit,("P2","2","0",2))
        push!(circuit,("R2","2","0",R2v))
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
        names = ["C1","C2","R1"]
        sensS = hblinsolve(ws, circuit, defs; keyedarrays=false,
            sensitivitynames=names, returnSsensitivity=true, returnS=true)
        sensnoS = hblinsolve(ws, circuit, defs; keyedarrays=false,
            sensitivitynames=names, returnSsensitivity=true, returnS=false)
        @test isempty(sensnoS.S)
        @test !isempty(sensnoS.Ssensitivity)
        @test sensnoS.Ssensitivity == sensS.Ssensitivity
    end

end

@testset "the frequency window of the pump modes" begin
    JC = JosephsonCircuits
    circuit = Tuple{String,String,String,Union{Complex{Float64},Symbol,Int64}}[]
    push!(circuit, ("P1","1","0",1)); push!(circuit, ("R1","1","0",:R))
    for i in 1:6
        push!(circuit, ("Lj$(i)","$(i)","$(i+1)",:Lj))
        push!(circuit, ("C$(i)","$(i)","0",:Cg))
    end
    push!(circuit, ("C7","7","0",:Cg)); push!(circuit, ("R2","7","0",:R))
    defs = Dict{Symbol,Complex{Float64}}(:Lj => 100e-12, :Cg => 40e-15, :R => 50.0)
    w = (2*pi*5.0e9, 2*pi*1.19e9)
    src = [(mode=(1,0),port=1,current=0.6e-6), (mode=(0,1),port=1,current=0.6e-6)]
    full = JC.hbnlsolve(w, (8,4), src, circuit, defs; dc = true, odd = true, even = true,
        method = :newton, keyedarrays = false)
    @test full.solverinfo.converged
    # a floor below every retained frequency changes nothing
    same = JC.hbnlsolve(w, (8,4), src, circuit, defs; dc = true, odd = true, even = true,
        method = :newton, keyedarrays = false, frequencywindow = (2*pi*0.1e9, Inf))
    @test same.modes == full.modes
    @test same.nodeflux == full.nodeflux
    # a box removes modes and leaves the strong ones within the size of what it dropped
    box = JC.hbnlsolve(w, (8,4), src, circuit, defs; dc = true, odd = true, even = true,
        method = :newton, keyedarrays = false, frequencywindow = (2*pi*0.5e9, 2*pi*30e9))
    @test box.solverinfo.converged
    @test length(box.modes) < length(full.modes)
    @test all(m -> all(==(0), m) || 2*pi*0.5e9 <= abs(sum(w .* m)) <= 2*pi*30e9, box.modes)
    # a window that drops a source mode is refused with a message naming it
    @test_throws ArgumentError JC.hbnlsolve(w, (8,4), src, circuit, defs; dc = true, odd = true,
        even = true, method = :newton, keyedarrays = false, frequencywindow = (2*pi*2e9, Inf))
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
    builder(; Lj) = [(n, a, b, v isa Symbol ? (v === :Lj ? Lj : defs[v]) : v) for (n, a, b, v) in circuit]
    cache = JC.hbcache(w, (8,4), src, builder, (; Lj = 100e-12); dc = true, odd = true, even = true,
        frequencywindow = (2*pi*0.5e9, 2*pi*30e9), method = :newton)
    @test length(cache.frequencies.modes) == length(box.modes)
    st = JC.hbnlsolve(w, (8,4), src, circuit, defs; dc = true, odd = true, even = true,
        method = :staged, keyedarrays = false, frequencywindow = (2*pi*0.5e9, 2*pi*30e9))
    @test st.solverinfo.converged
    @test st.modes == box.modes
end

@testset "the evaluation grid of the pump modes" begin
    JC = JosephsonCircuits
    circuit = Tuple{String,String,String,Union{Complex{Float64},Symbol,Int64}}[]
    push!(circuit, ("P1","1","0",1)); push!(circuit, ("R1","1","0",:R))
    for i in 1:6
        push!(circuit, ("Lj$(i)","$(i)","$(i+1)",:Lj))
        push!(circuit, ("C$(i)","$(i)","0",:Cg))
    end
    push!(circuit, ("C7","7","0",:Cg)); push!(circuit, ("R2","7","0",:R))
    defs = Dict{Symbol,Complex{Float64}}(:Lj => 100e-12, :Cg => 40e-15, :R => 50.0)
    w = (2*pi*5.0e9, 2*pi*1.19e9)
    src = [(mode=(1,0),port=1,current=0.6e-6), (mode=(0,1),port=1,current=0.6e-6)]
    kw = (; dc = true, odd = true, even = true, method = :newton, keyedarrays = false)
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
    @test_throws ArgumentError JC.hbnlsolve(w, (8,4), src, circuit, defs; kw..., method = :staged, Nevaluationharmonics = (8,3))
    @test_throws ArgumentError JC.hbsolve(2*pi*5.1e9, w, src, (1,1), (8,4), circuit, defs; dc = true,
        threewavemixing = true, fourwavemixing = true, keyedarrays = false, Nevaluationharmonics = (7,4))
    make(; Lj, Cg, R) = [(String(n), a, b, v isa Symbol ? Dict(:Lj => Lj, :Cg => Cg, :R => R)[v] : v)
        for (n, a, b, v) in circuit]
    @test_throws ArgumentError JC.hbcache(w, (8,4), src, make, (Lj = 100e-12, Cg = 40e-15, R = 50.0);
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
    st = JC.hbnlsolve(w, (8,4), src, circuit, defs; kw..., method = :staged, Nevaluationharmonics = (24,12))
    @test st.frequencies.Nharmonics == (24,12)
    @test isapprox(st.nodeflux, wide.nodeflux; rtol = 1e-6)
    cache = JC.hbcache(w, (8,4), src, make, (Lj = 100e-12, Cg = 40e-15, R = 50.0);
        dc = true, odd = true, even = true, Nevaluationharmonics = (24,12), keyedarrays = false)
    @test cache.frequencies.Nharmonics == (24,12)
    cached = JC.hbsolve!(cache, (Lj = 100e-12, Cg = 40e-15, R = 50.0))
    @test isapprox(cached.nodeflux, wide.nodeflux; rtol = 1e-6)
end
