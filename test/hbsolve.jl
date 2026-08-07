using JosephsonCircuits
using Test

@testset verbose=true "hbsolve" begin


    @testset "hbsolve-hbnlsolve comparison" begin

        @variables R Cc Lj Cj
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
                Npumpharmonics, circuit, circuitdefs, ftol = 1e-15)
            S1ss = sol1.linearized.S((0,),1,(0,),1,1)
            S1is = sol1.linearized.S((-2,),1,(0,),1,1)

            # nonlinear simulation with (pump,signal) order
            w = (wp,ws)
            Nharmonics = (10,10)
            sources = [(mode=(1,0),port=1,current=Ip),(mode=(0,1),port=1,current=Is)]
            sol2 = hbnlsolve(w, Nharmonics, sources, circuit, circuitdefs, ftol = 1e-15)
            S2ss = sol2.S((0,1),1,(0,1),1)
            S2is = sol2.S((2,-1),1,(0,1),1)

            # nonlinear simulation with (signal,pump) order
            w = (ws,wp)
            Nharmonics = (10,10)
            sources = [(mode=(0,1),port=1,current=Ip),(mode=(1,0),port=1,current=Is)]
            sol3 = hbnlsolve(w, Nharmonics, sources, circuit, circuitdefs, ftol = 1e-15)
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

    @testset "hbsolve lossless new syntax" begin

        @variables Rleft Cc Lj Cj w L1
        circuit = Tuple{String,String,String,Num}[]
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

        result = hbsolve(ws, wp, sources, Nmodulationharmonics,
            Npumpharmonics, circuit, circuitdefs, ftol=1e-12)

        @test isapprox(
            result.nonlinear.nodeflux[:],
            ComplexF64[-0.013189575486618105 - 0.00865077163136891im, 2.6396823809835998e-5 - 5.667772212498911e-6im, -2.3501629748681806e-8 - 2.2306559257777402e-8im, -7.987057959121183e-12 + 3.750064665528795e-11im, 4.337129743740559e-14 - 1.3397330127340191e-14im, -3.812185723603184e-17 - 3.8456119951452504e-17im, -4.278102807441037e-20 + 3.974761087481868e-20im, -8.24749227891514e-25 + 3.0560138516751506e-24im, 0.12157593753208869 + 0.07973582696437984im, 1.3736443951855489e-5 - 6.463164795537436e-5im, -5.3397980865816756e-8 + 9.191484285021433e-9im, 2.7913096743299775e-11 + 4.514682457284182e-11im, 3.3395873694437265e-14 - 4.569085921919452e-14im, -6.154948861291708e-17 - 1.523212283849711e-17im, -2.2291909276375226e-20 + 6.18003971235134e-20im, 5.405242786243393e-25 + 3.424470385254823e-24im],
            atol = 1e-8)

        @test isapprox(
            10*log10.(abs2.(result.linearized.S((0,),1,(0,),1,:))),
            [0.002717249265375112, 0.0031881165263016476, 0.003764438371177404, 0.004475612824820614, 0.005360972358329889, 0.006473775210238347, 0.007887075336301673, 0.009702494973065136, 0.012063564112722006, 0.0151763872044898, 0.019342315628701326, 0.025010742306207902, 0.03286645328833334, 0.043977897966981386, 0.060055891091184574, 0.08391851842302066, 0.12035318288996377, 0.17776810664422604, 0.27146258097536335, 0.43031121414540996, 0.7108074673012084, 1.2271520892520107, 2.2163506977088527, 4.180249149477177, 8.16945893982768, 13.302337843265873, 8.180459085825893, 4.185698117374451, 2.2190409577860577, 1.2285225612890744, 0.7115311472717254, 0.4307082670317731, 0.2716887696865355, 0.17790149553267967, 0.12043426616790562, 0.08396907049655332, 0.06008804086742857, 0.04399863098219145, 0.0328799192406124, 0.02501947707672059, 0.019347909243977855, 0.015179860051320881, 0.012065584850116115, 0.009703510073587402, 0.007887387964370132, 0.006473594325612282, 0.005360443913078622, 0.0044748399587568825, 0.003763494713109064, 0.0031870550949444913, 0.002716108513561703],
            atol = 1e-6)

        # test some uncommon options
        @variables w
        result = hbsolve(ws, wp, sources, Nmodulationharmonics,
            Npumpharmonics, circuit, circuitdefs, ftol=1e-12, symfreqvar = w,
            returnS=false, returnSnoise=true, returnQE=false,
            returnnodeflux=true,
            returnnodefluxadjoint=true, returnCM=false,
            returnvoltage=true, returnvoltageadjoint=true,
            returnSsensitivity = true,
            returnZsensitivity=true, returnZsensitivityadjoint=true,
            sensitivitynames=["C1"],
            returnZ = true, returnZadjoint = true,
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

        @variables R Cc Lj Cj
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
            "Ssensitivity","Z","Zadjoint","Zsensitivity","Zsensitivityadjoint"]

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

        @variables Rleft Cc Lj Cj w L1
        circuit = Tuple{String,String,String,Num}[]
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
            "Solver did not converge after maximum iterations of 1.",
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

    @testset "hbnlsolve simple testcase error" begin

        circuit = [("P1","1","0",1),("R1","1","0",50.0)]
        circuitdefs = Dict()
        Idc = 50e-5
        Ip = 1.0e-6
        wp=2*pi*5.0*1e9
        Npumpmodes = 1
        @test_throws(
            ErrorException("NaN in nonlinear solver."),
            hbnlsolve((wp,),(Npumpmodes,),[(mode=(1,),port=1,current=Ip),],
                circuit,circuitdefs;dc=true,odd=true,even=false))
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
            ArgumentError("Source mode (2,) not found."),
            JosephsonCircuits.calcsources(modes, sources, portindices, portnumbers,
                nodeindices, edge2indexdict, Lmean, Nnodes, Nbranches, Nmodes))

    end

end
