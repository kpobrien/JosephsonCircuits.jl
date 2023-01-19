using JosephsonCircuits
using Test

@testset verbose=true "hbsolve2" begin

    @testset "hbsolve2 lossless" begin

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
        result=hbsolve2(ws, wp, Ip, Nsignalmodes, Npumpmodes, circuit, circuitdefs, pumpports=[1])

        @test isapprox(
            result.pump.nodeflux,
            ComplexF64[-0.013189715720865495 - 0.008651167511265932im, 2.639595341567038e-5 - 5.666452835660498e-6im, -2.34984668856955e-8 - 2.2306610719944793e-8im, -7.989866762787069e-12 + 3.749622154543886e-11im, 4.336736558999377e-14 - 1.339035559272133e-14im, -3.8074599154274697e-17 - 3.8493874617243354e-17im, -5.4742642304096124e-21 + 5.311896333888763e-20im, 2.7877342893915465e-20 + 1.40385643300831e-20im, 0.12157314440900607 + 0.0797363708287639im, 1.373852071423108e-5 - 6.462838433653931e-5im, -5.3394886974965865e-8 + 9.187193797853318e-9im, 2.7906051689266846e-11 + 4.5145088389588067e-11im, 3.339713496666733e-14 - 4.5680957097985434e-14im, -6.152523083368712e-17 - 1.5298667328128304e-17im, 2.19075261797519e-20 + 5.594083986630819e-20im, 3.41490683767286e-20 + 1.5843676109621666e-21im],
            atol = 1e-8)

        @test isapprox(
            10*log10.(abs2.(result.signal.S[result.signal.signalindex,result.signal.signalindex,:])),
            [0.0027171001162159705, 0.0031879418548192346, 0.00376423255414956, 0.004475368703477247, 0.005360680732607021, 0.006473424136187171, 0.007886649137566306, 0.009701972830963664, 0.012062918026015033, 0.015175578990910845, 0.019341292444897646, 0.025009429835687373, 0.03286474518958955, 0.04397563920753866, 0.06005285108782671, 0.08391434685195018, 0.12034733547692913, 0.17775971803082047, 0.27145024266250045, 0.4302925781788617, 0.7107785121246843, 1.2271056200958508, 2.216272404072828, 4.1801014508670296, 8.169071205571798, 13.301089937850207, 8.180070357250763, 4.185550179545209, 2.2189625576398537, 1.2284760310770502, 0.7115021541352169, 0.4306896066075445, 0.2716764151818288, 0.17789309591540733, 0.12042841108231785, 0.08396489344330892, 0.06008499685718375, 0.04399636923319631, 0.032878208869787345, 0.0250181628508142, 0.019346884684173782, 0.015179050745275594, 0.012064937886309808, 0.009702987220389104, 0.007886961184202037, 0.0064732427726544774, 0.005360151890417837, 0.004474595506671817, 0.0037632886192849233, 0.0031868801909932577, 0.002715959168639594],
            atol = 1e-8)
    end

    @testset "hbsolve2 hbsolve comparison" begin

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
        result1=hbsolve(ws, wp, Ip, Nsignalmodes, Npumpmodes, circuit, circuitdefs, pumpports=[1], ftol=1e-14)
        result2=hbsolve2(ws, wp, Ip, Nsignalmodes, Npumpmodes, circuit, circuitdefs, pumpports=[1], ftol=1e-14)

        @test isapprox(
            result1.pump.nodeflux,
            result2.pump.nodeflux,
            atol = 1e-8)

        @test isapprox(
            result1.signal.S[result1.signal.signalindex,result1.signal.signalindex,:],
            result2.signal.S[result2.signal.signalindex,result2.signal.signalindex,:],
            atol = 1e-8)
    end

    @testset "hbsolve2 lossy" begin

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
        result=hbsolve2(ws, wp, Ip, Nsignalmodes, Npumpmodes, circuit, circuitdefs,pumpports=[1])

        ## There is something wrong with hbsolve2 because the difference with v1 is
        # so large
        @test isapprox(
            result.pump.nodeflux,
            ComplexF64[-0.005404178815221458 - 0.020370051121134776im, 6.815437530731759e-7 - 1.147015781026838e-6im, -1.7534339362832887e-10 + 1.2486893692222313e-10im, 3.3880455752540875e-14 - 5.129258669286749e-15im, -5.098160353117816e-18 - 1.9161091589686673e-18im, -1.3169582712868142e-21 + 2.378490670885394e-21im, -1.74210248492438e-20 - 1.4464020662586226e-21im, 2.2026817534197045e-25 + 4.844501217305843e-25im, 0.05082755173685041 + 0.015844684081317333im, -1.8806014684433024e-6 - 2.669413379277482e-6im, -7.987884483700569e-12 + 3.5987280288845466e-10im, 2.8970108967623174e-14 - 3.756372877550381e-14im, -6.5248625853400994e-18 + 1.8798941751708544e-18im, 1.320283332174409e-22 + 3.1807872353177106e-21im, -1.8166616926680373e-20 + 7.533796739748396e-21im, 4.36696159156451e-25 + 3.8604536353286075e-25im],
            atol = 1e-8)

        @test isapprox(
            10*log10.(abs2.(result.signal.S[result.signal.signalindex,result.signal.signalindex,:])),
            [-0.35881759544600905, -0.38511742360391643, -0.41414939769688863, -0.4462870003283318, -0.4819671497880933, -0.5217030287748109, -0.5660998955861667, -0.6158746474948269, -0.6718801163264685, -0.7351353415766995, -0.8068633957010887, -0.8885387334429317, -0.9819464925958608, -1.089256647350981, -1.213116305554448, -1.3567635288328814, -1.5241654068346333, -1.7201809136275068, -1.9507437701627133, -2.2230491483972665, -2.5457046332732225, -2.9287583947682676, -3.383422108209453, -3.92111711734785, -4.551109819067609, -5.27536270362385, -6.078361815894438, -6.909696265201516, -7.66286700663108, -8.174297594526012, -8.287679622361011, -7.968447777549189, -7.3340361858580945, -6.555662042628274, -5.764284332567529, -5.030262404295716, -4.380877509200552, -3.820225946614721, -3.34195798331703, -2.935990188525385, -2.59168467385448, -2.2992025808827408, -2.0499724131703427, -1.836752763172862, -1.6535251019382815, -1.4953289751470045, -1.3580907832498634, -1.2384674508002553, -1.1337120454056717, -1.0415619186205975, -0.9601472569056663],
            atol = 1e-8)
    end


    @testset "hbsolve2 hbsolve lossy comparison" begin

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
        result1=hbsolve(ws, wp, Ip, Nsignalmodes, Npumpmodes, circuit, circuitdefs, pumpports=[1], ftol=1e-14)
        result2=hbsolve2(ws, wp, Ip, Nsignalmodes, Npumpmodes, circuit, circuitdefs, pumpports=[1], ftol=1e-14)

        @test isapprox(
            result1.pump.nodeflux,
            result2.pump.nodeflux,
            atol = 1e-8)

        @test isapprox(
            result1.signal.S[result1.signal.signalindex,result1.signal.signalindex,:],
            result2.signal.S[result2.signal.signalindex,result2.signal.signalindex,:],
            atol = 1e-8)
    end

end