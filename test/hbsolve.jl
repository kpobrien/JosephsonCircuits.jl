using JosephsonCircuits
using Test

@testset verbose=true "hbsolve" begin

    @testset "hbsolve lossless" begin

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
        result=hbsolve(ws, wp, Ip, Nsignalmodes, Npumpmodes, circuit,
            circuitdefs, pumpports=[1], ftol=1e-12)

        @test isapprox(
            result.pump.nodeflux,
            ComplexF64[-0.013189575486618105 - 0.00865077163136891im, 2.6396823809835998e-5 - 5.667772212498911e-6im, -2.3501629748681806e-8 - 2.2306559257777402e-8im, -7.987057959121183e-12 + 3.750064665528795e-11im, 4.337129743740559e-14 - 1.3397330127340191e-14im, -3.812185723603184e-17 - 3.8456119951452504e-17im, -4.278102807441037e-20 + 3.974761087481868e-20im, -8.24749227891514e-25 + 3.0560138516751506e-24im, 0.12157593753208869 + 0.07973582696437984im, 1.3736443951855489e-5 - 6.463164795537436e-5im, -5.3397980865816756e-8 + 9.191484285021433e-9im, 2.7913096743299775e-11 + 4.514682457284182e-11im, 3.3395873694437265e-14 - 4.569085921919452e-14im, -6.154948861291708e-17 - 1.523212283849711e-17im, -2.2291909276375226e-20 + 6.18003971235134e-20im, 5.405242786243393e-25 + 3.424470385254823e-24im],
            atol = 1e-8)

        @test isapprox(
            # 10*log10.(abs2.(result.signal.S[result.signal.signalindex,result.signal.signalindex,:])),
            10*log10.(abs2.(result.signal.S[1,1,:])),
            [0.002717249265375112, 0.0031881165263016476, 0.003764438371177404, 0.004475612824820614, 0.005360972358329889, 0.006473775210238347, 0.007887075336301673, 0.009702494973065136, 0.012063564112722006, 0.0151763872044898, 0.019342315628701326, 0.025010742306207902, 0.03286645328833334, 0.043977897966981386, 0.060055891091184574, 0.08391851842302066, 0.12035318288996377, 0.17776810664422604, 0.27146258097536335, 0.43031121414540996, 0.7108074673012084, 1.2271520892520107, 2.2163506977088527, 4.180249149477177, 8.16945893982768, 13.302337843265873, 8.180459085825893, 4.185698117374451, 2.2190409577860577, 1.2285225612890744, 0.7115311472717254, 0.4307082670317731, 0.2716887696865355, 0.17790149553267967, 0.12043426616790562, 0.08396907049655332, 0.06008804086742857, 0.04399863098219145, 0.0328799192406124, 0.02501947707672059, 0.019347909243977855, 0.015179860051320881, 0.012065584850116115, 0.009703510073587402, 0.007887387964370132, 0.006473594325612282, 0.005360443913078622, 0.0044748399587568825, 0.003763494713109064, 0.0031870550949444913, 0.002716108513561703],
            atol = 1e-6)
    end

    @testset "hbsolve lossless two ports" begin

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
        result=hbsolve(ws, wp, [Ip/2,Ip/2], Nsignalmodes, Npumpmodes, circuit,
            circuitdefs, pumpports=[1,1], ftol=1e-12)

        @test isapprox(
            result.pump.nodeflux,
            ComplexF64[-0.013189575486618105 - 0.00865077163136891im, 2.6396823809835998e-5 - 5.667772212498911e-6im, -2.3501629748681806e-8 - 2.2306559257777402e-8im, -7.987057959121183e-12 + 3.750064665528795e-11im, 4.337129743740559e-14 - 1.3397330127340191e-14im, -3.812185723603184e-17 - 3.8456119951452504e-17im, -4.278102807441037e-20 + 3.974761087481868e-20im, -8.24749227891514e-25 + 3.0560138516751506e-24im, 0.12157593753208869 + 0.07973582696437984im, 1.3736443951855489e-5 - 6.463164795537436e-5im, -5.3397980865816756e-8 + 9.191484285021433e-9im, 2.7913096743299775e-11 + 4.514682457284182e-11im, 3.3395873694437265e-14 - 4.569085921919452e-14im, -6.154948861291708e-17 - 1.523212283849711e-17im, -2.2291909276375226e-20 + 6.18003971235134e-20im, 5.405242786243393e-25 + 3.424470385254823e-24im],
            atol = 1e-8)

        @test isapprox(
            # 10*log10.(abs2.(result.signal.S[result.signal.signalindex,result.signal.signalindex,:])),
            10*log10.(abs2.(result.signal.S[1,1,:])),
            [0.002717249265375112, 0.0031881165263016476, 0.003764438371177404, 0.004475612824820614, 0.005360972358329889, 0.006473775210238347, 0.007887075336301673, 0.009702494973065136, 0.012063564112722006, 0.0151763872044898, 0.019342315628701326, 0.025010742306207902, 0.03286645328833334, 0.043977897966981386, 0.060055891091184574, 0.08391851842302066, 0.12035318288996377, 0.17776810664422604, 0.27146258097536335, 0.43031121414540996, 0.7108074673012084, 1.2271520892520107, 2.2163506977088527, 4.180249149477177, 8.16945893982768, 13.302337843265873, 8.180459085825893, 4.185698117374451, 2.2190409577860577, 1.2285225612890744, 0.7115311472717254, 0.4307082670317731, 0.2716887696865355, 0.17790149553267967, 0.12043426616790562, 0.08396907049655332, 0.06008804086742857, 0.04399863098219145, 0.0328799192406124, 0.02501947707672059, 0.019347909243977855, 0.015179860051320881, 0.012065584850116115, 0.009703510073587402, 0.007887387964370132, 0.006473594325612282, 0.005360443913078622, 0.0044748399587568825, 0.003763494713109064, 0.0031870550949444913, 0.002716108513561703],
            atol = 1e-6)
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
            result.pump.nodeflux[:],
            ComplexF64[-0.013189575486618105 - 0.00865077163136891im, 2.6396823809835998e-5 - 5.667772212498911e-6im, -2.3501629748681806e-8 - 2.2306559257777402e-8im, -7.987057959121183e-12 + 3.750064665528795e-11im, 4.337129743740559e-14 - 1.3397330127340191e-14im, -3.812185723603184e-17 - 3.8456119951452504e-17im, -4.278102807441037e-20 + 3.974761087481868e-20im, -8.24749227891514e-25 + 3.0560138516751506e-24im, 0.12157593753208869 + 0.07973582696437984im, 1.3736443951855489e-5 - 6.463164795537436e-5im, -5.3397980865816756e-8 + 9.191484285021433e-9im, 2.7913096743299775e-11 + 4.514682457284182e-11im, 3.3395873694437265e-14 - 4.569085921919452e-14im, -6.154948861291708e-17 - 1.523212283849711e-17im, -2.2291909276375226e-20 + 6.18003971235134e-20im, 5.405242786243393e-25 + 3.424470385254823e-24im],
            atol = 1e-8)

        @test isapprox(
            10*log10.(abs2.(result.signal.S((0,),1,(0,),1,:))),
            [0.002717249265375112, 0.0031881165263016476, 0.003764438371177404, 0.004475612824820614, 0.005360972358329889, 0.006473775210238347, 0.007887075336301673, 0.009702494973065136, 0.012063564112722006, 0.0151763872044898, 0.019342315628701326, 0.025010742306207902, 0.03286645328833334, 0.043977897966981386, 0.060055891091184574, 0.08391851842302066, 0.12035318288996377, 0.17776810664422604, 0.27146258097536335, 0.43031121414540996, 0.7108074673012084, 1.2271520892520107, 2.2163506977088527, 4.180249149477177, 8.16945893982768, 13.302337843265873, 8.180459085825893, 4.185698117374451, 2.2190409577860577, 1.2285225612890744, 0.7115311472717254, 0.4307082670317731, 0.2716887696865355, 0.17790149553267967, 0.12043426616790562, 0.08396907049655332, 0.06008804086742857, 0.04399863098219145, 0.0328799192406124, 0.02501947707672059, 0.019347909243977855, 0.015179860051320881, 0.012065584850116115, 0.009703510073587402, 0.007887387964370132, 0.006473594325612282, 0.005360443913078622, 0.0044748399587568825, 0.003763494713109064, 0.0031870550949444913, 0.002716108513561703],
            atol = 1e-6)

        # test some uncommon options
        result = hbsolve(ws, wp, sources, Nmodulationharmonics,
            Npumpharmonics, circuit, circuitdefs, ftol=1e-12,
            returnS=false, returnSnoise=true, returnQE=false,
            returnnodefluxadjoint=true, returnCM=false,
            returnvoltage=true, nbatches=4)

        @test result.signal.QE == Array{Float64, 3}(undef, 0, 0, 0)
        @test result.signal.S == Array{Float64, 3}(undef, 0, 0, 0)
        @test result.signal.nodeflux == ComplexF64[]
        @test result.signal.Snoise[:] == ComplexF64[]
        @test result.signal.CM == Matrix{Float64}(undef, 0, 0)

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

    @testset "hbsolve hbsolveold comparison" begin

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

        # reduce the tolerance in order to get consistent results between different solvers
        result1=hbsolveold(ws, wp, Ip, Nsignalmodes, Npumpmodes, circuit, circuitdefs, pumpports=[1], ftol=1e-14)
        result2=hbsolve(ws, wp, Ip, Nsignalmodes, Npumpmodes, circuit, circuitdefs, pumpports=[1], ftol=1e-14)

        @test isapprox(
            result1.pump.nodeflux,
            result2.pump.nodeflux,
            atol = 1e-8)

        @test isapprox(
            result1.signal.S[result1.signal.signalindex,result1.signal.signalindex,:],
            result2.signal.S[1,1,:],
            # result2.signal.S[result2.signal.signalindex,result2.signal.signalindex,:],
            atol = 1e-8)
    end

    @testset "hbsolve lossy" begin

        @variables Rleft Cc Lj Cj w L1 Rloss
        circuit = Tuple{String,String,String,Num}[]
        push!(circuit,("P1","1","0",1))
        push!(circuit,("R1","1","0",Rleft))
        push!(circuit,("C1","1","2",Cc)) 
        push!(circuit,("Lj1","2","0",Lj)) 
        push!(circuit,("C2","2","0",Cj))
        push!(circuit,("R2","2","0",Rloss))

        circuitdefs = Dict(
            Lj =>1000.0e-12,
            Cc => 100.0e-15,
            Cj => 1000.0e-15,
            Rloss => 1000.0,
            Rleft => 50.0,
        )
        ws = 2*pi*(4.5:0.01:5.0)*1e9
        wp = 2*pi*4.75001*1e9
        Ip = 0.00565e-6
        Nsignalmodes = 8
        Npumpmodes = 8
        result=hbsolve(ws, wp, Ip, Nsignalmodes, Npumpmodes, circuit,
            circuitdefs,pumpports=[1], ftol=1e-12)

        ## There is something wrong with hbsolve2 because the difference with v1 is
        # so large
        @test isapprox(
            result.pump.nodeflux,
            ComplexF64[-0.005404178767784759 - 0.02037005108680524im, 6.815437369745635e-7 - 1.1470158059579708e-6im, -1.7534339277181034e-10 + 1.2486894513346514e-10im, 3.38803774302399e-14 - 5.129245778951382e-15im, -5.137203687009167e-18 - 1.9118072527307488e-18im, -1.205634914643852e-20 - 2.015256691761906e-20im, -3.2555615922027985e-21 - 1.5327643629048948e-20im, 2.0932125762839933e-20 - 7.109266968211845e-21im, 0.05082755201433778 + 0.015844683797761848im, -1.8806015402318049e-6 - 2.669413368248357e-6im, -7.987872622070442e-12 + 3.598728099517483e-10im, 2.897004298551031e-14 - 3.756364090558882e-14im, -6.560702793144546e-18 + 1.913267082586061e-18im, -2.4333378384957925e-20 - 1.2807787883479442e-20im, -1.1156663175354565e-20 - 1.364946507318454e-20im, 1.7756062330690395e-20 - 1.6460689483328482e-20im],
            atol = 1e-8)

        @test isapprox(
            # 10*log10.(abs2.(result.signal.S[result.signal.signalindex,result.signal.signalindex,:])),
            10*log10.(abs2.(result.signal.S[1,1,:])),

             [-0.35881759555167236, -0.3851174237205393, -0.41414939782595106, -0.4462870004715635, -0.48196714994752377, -0.5217030289528287, -0.5660998957855977, -0.615874647719035, -0.6718801165794636, -0.7351353418633039, -0.8068633960271241, -0.8885387338154525, -0.9819464930235025, -1.0892566478443222, -1.2131163061265446, -1.3567635294998903, -1.5241654076166113, -1.7201809145493425, -1.9507437712551328, -2.223049149697846, -2.545704634827017, -2.9287583966275106, -3.3834221104310536, -3.9211171199862456, -4.551109822158283, -5.275362707148965, -6.078361819718249, -6.90969626897045, -7.662867009688683, -8.174297596028158, -8.287679621767847, -7.968447775092126, -7.33403618231908, -6.555662038768009, -5.76428432885789, -5.030262400951592, -4.380877506283952, -3.820225944113212, -3.3419579811879623, -2.935990186717743, -2.5916846723187135, -2.2992025795746533, -2.049972412052101, -1.8367527622127993, -1.6535251011101937, -1.4953289744293266, -1.3580907826248823, -1.2384674502534234, -1.1337120449250095, -1.0415619181962141, -0.9601472565293714],
            atol = 1e-8)
    end


    @testset "hbsolve hbsolveold lossy comparison" begin

        @variables Rleft Cc Lj Cj w L1 Rloss
        circuit = Tuple{String,String,String,Num}[]
        push!(circuit,("P1","1","0",1))
        push!(circuit,("R1","1","0",Rleft))
        push!(circuit,("C1","1","2",Cc)) 
        push!(circuit,("Lj1","2","0",Lj)) 
        push!(circuit,("C2","2","0",Cj))
        push!(circuit,("R2","2","0",Rloss))

        circuitdefs = Dict(
            Lj =>1000.0e-12,
            Cc => 100.0e-15,
            Cj => 1000.0e-15,
            Rloss => 1000.0,
            Rleft => 50.0,
        )
        ws = 2*pi*(4.5:0.01:5.0)*1e9
        wp = 2*pi*4.75001*1e9
        Ip = 0.00565e-6
        Nsignalmodes = 8
        Npumpmodes = 8
        
        # reduce the tolerance in order to get consistent results between different solvers
        result1=hbsolveold(ws, wp, Ip, Nsignalmodes, Npumpmodes, circuit, circuitdefs, pumpports=[1], ftol=1e-14)
        result2=hbsolve(ws, wp, Ip, Nsignalmodes, Npumpmodes, circuit, circuitdefs, pumpports=[1], ftol=1e-14)

        @test isapprox(
            result1.pump.nodeflux,
            result2.pump.nodeflux,
            atol = 1e-8)

        @test isapprox(
            result1.signal.S[result1.signal.signalindex,result1.signal.signalindex,:],
            # result2.signal.S[result2.signal.signalindex,result2.signal.signalindex,:],
            result2.signal.S[1,1,:],
            atol = 1e-8)
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