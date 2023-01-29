using JosephsonCircuits
using Test

@testset verbose=true "hbsolve" begin


    @testset verbose=true "hbsolve lossless" begin

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
        result=hbsolve(ws, wp, Ip, Nsignalmodes, Npumpmodes, circuit, circuitdefs,pumpports=[1])

        @test isapprox(
            result.pump.nodeflux,
            ComplexF64[-0.013189736218636441 - 0.008651273221683266im, 2.6395682539011797e-5 - 5.6661514233313726e-6im, -2.3497663047131403e-8 - 2.2306503448330988e-8im, -7.990410788185173e-12 + 3.749499317166309e-11im, 4.336610521409153e-14 - 1.3388741572794773e-14im, -3.8059160268943944e-17 - 3.852723364098488e-17im, 6.07942296115915e-20 + 6.60227155809359e-21im, 5.977496795928632e-20 + 2.938144480595029e-20im, 0.12157241551965468 + 0.07973640247896553im, 1.3738923117026488e-5 - 6.46274778537762e-5im, -5.3393939365693146e-8 + 9.186223725369537e-9im, 2.790433171591809e-11 + 4.514438082272354e-11im, 3.3397076362502e-14 - 4.567840462363693e-14im, -6.153011440710407e-17 - 1.534143178629367e-17im, 6.419757193364206e-20 - 2.4735969038803813e-20im, 7.29011218720814e-20 + 2.6769912044348367e-21im],
            atol = 1e-8)

        @test isapprox(
            10*log10.(abs2.(result.signal.S[result.signal.signalindex,result.signal.signalindex,:])),
            [0.0027170567684318336, 0.0031878910893791113, 0.003764172736760887, 0.0044752977535568025, 0.0053605959763061576, 0.0064733221021938975, 0.00788652526979293, 0.009701821078803697, 0.012062730251351263, 0.015175344096677876, 0.019340995073025724, 0.025009048387281324, 0.03286424875818319, 0.04397498273540677, 0.06005196755960492, 0.08391313445153363, 0.12034563601945997, 0.17775728001306818, 0.2714466567248869, 0.43028716192043454, 0.7107700967255848, 1.2270921144417903, 2.2162496487671177, 4.180058522554849, 8.168958508626085, 13.30072725995452, 8.179957371303482, 4.185507181701968, 2.2189397713773675, 1.228462507677842, 0.7114937277036677, 0.4306841832409538, 0.27167282453838204, 0.1778906546996149, 0.12042670939499213, 0.08396367944962557, 0.06008411216445837, 0.04399571189222961, 0.032877711778053494, 0.02501778089224741, 0.019346586912408068, 0.015178815533550454, 0.01206474985673973, 0.009702835261564286, 0.007886837147456062, 0.00647314059948161, 0.005360067018761243, 0.0044745244606336065, 0.0037632287214528454, 0.0031868293579977005, 0.0027159157639621524],
            atol = 1e-8)

    end

    @testset verbose=true "hbsolve lossless nodeflux, no QE, CM" begin

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
        result=hbsolve(ws, wp, Ip, Nsignalmodes, Npumpmodes, circuit, circuitdefs,
            pumpports=[1],returnnodeflux=true,returnS = false, returnQE=false,returnCM=false,
            returnvoltage=true,returnSnoise=true)

        @test isapprox(
            result.pump.nodeflux,
            ComplexF64[-0.013189736218636441 - 0.008651273221683266im, 2.6395682539011797e-5 - 5.6661514233313726e-6im, -2.3497663047131403e-8 - 2.2306503448330988e-8im, -7.990410788185173e-12 + 3.749499317166309e-11im, 4.336610521409153e-14 - 1.3388741572794773e-14im, -3.8059160268943944e-17 - 3.852723364098488e-17im, 6.07942296115915e-20 + 6.60227155809359e-21im, 5.977496795928632e-20 + 2.938144480595029e-20im, 0.12157241551965468 + 0.07973640247896553im, 1.3738923117026488e-5 - 6.46274778537762e-5im, -5.3393939365693146e-8 + 9.186223725369537e-9im, 2.790433171591809e-11 + 4.514438082272354e-11im, 3.3397076362502e-14 - 4.567840462363693e-14im, -6.153011440710407e-17 - 1.534143178629367e-17im, 6.419757193364206e-20 - 2.4735969038803813e-20im, 7.29011218720814e-20 + 2.6769912044348367e-21im],
            atol = 1e-8)


    end

    @testset verbose=true "hbsolve lossy" begin

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
            circuitdefs,pumpports=[1], symfreqvar = w)

        @test isapprox(
            result.pump.nodeflux,
            ComplexF64[-0.0054041783971746 - 0.020370057382730806im, 6.815424075183076e-7 - 1.147012900366793e-6im, -1.753427490430056e-10 + 1.2486837493074984e-10im, 3.388035426864977e-14 - 5.129280878553909e-15im, -5.1298428933921494e-18 - 1.739836240646162e-18im, -4.596535957104933e-20 - 7.912243213092834e-20im, 2.0030609892972045e-20 + 5.0036217228674835e-20im, 1.794196518026824e-19 + 8.871654522001758e-20im, 0.05082751019439782 + 0.015844675018286195im, -1.8805963793270084e-6 - 2.6694074929857497e-6im, -7.987993107073292e-12 + 3.598713769919333e-10im, 2.8969986222335313e-14 - 3.756365383211586e-14im, -6.4252953140502435e-18 + 2.0797573731599905e-18im, -9.416708080413478e-20 - 5.112014012913594e-20im, 4.582330476037865e-20 + 3.9710828188645537e-20im, 2.190537485892898e-19 + 8.560855801192207e-21im],
            atol = 1e-8)

        @test isapprox(
            10*log10.(abs2.(result.signal.S[result.signal.signalindex,result.signal.signalindex,:])),
            [-0.35881757065389075, -0.3851173962405351, -0.41414936741470204, -0.44628696672146195, -0.48196711238056517, -0.52170298700604, -0.56609984879292, -0.6158745948887158, -0.6718800569656551, -0.7351352743301984, -0.8068633192034746, -0.8885386460366157, -0.9819463922574563, -1.089256531597212, -1.2131161713226868, -1.3567633723310317, -1.5241652233566767, -1.7201806973353808, -1.9507435138455569, -2.223048843239626, -2.545704268702978, -2.9287579585295886, -3.38342158695043, -3.9211164982952336, -4.551109093896298, -5.275361876519062, -6.078360918705144, -6.909695380887424, -7.662866289218854, -8.174297242072996, -8.287679761537564, -7.968448354057507, -7.334037016227039, -6.555662948371419, -5.764285202966598, -5.030263188934635, -4.380878193531114, -3.820226533549771, -3.341958482864788, -2.9359906126572595, -2.5916850341946667, -2.2992028878025774, -2.049972675545674, -1.8367529884342717, -1.6535252962345233, -1.4953291435376725, -1.3580909298903507, -1.2384675791047277, -1.133712158184686, -1.0415620181945617, -0.9601473451970821],
            atol = 1e-8)

    end

    @testset verbose=true "hbsolve lossy freq dep" begin

        @variables Rleft Cc Lj Cj w L1 Rloss
        circuit = Tuple{String,String,String,Num}[]
        push!(circuit,("P1","1","0",1))
        push!(circuit,("R1","1","0",Rleft))
        push!(circuit,("C1","1","2",Cc)) 
        push!(circuit,("Lj1","2","0",Lj)) 
        push!(circuit,("C2","2","0",Cj))
        push!(circuit,("R2","2","0",w/(w+1)*Rloss))

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
            circuitdefs,pumpports=[1], symfreqvar = w, returnSnoise = true)

        @test isapprox(
            result.pump.nodeflux,
            ComplexF64[-0.005404178397088291 - 0.020370057382871846im, 6.815424074495673e-7 - 1.147012900315888e-6im, -1.7534274887574658e-10 + 1.2486837491319458e-10im, 3.3880300476756544e-14 - 5.12912465901757e-15im, -5.1175038715858325e-18 - 1.6622598304465073e-18im, 6.573030899482611e-20 - 8.485356182489419e-20im, 3.731089679352745e-20 - 2.916242040788456e-20im, -6.028650696017834e-20 - 2.983335821764548e-20im, 0.05082751019353899 + 0.01584467501756678im, -1.8805963792820397e-6 - 2.669407492781297e-6im, -7.987992963342683e-12 + 3.5987137675020927e-10im, 2.8970081982677106e-14 - 3.756344611647316e-14im, -6.355194218991686e-18 + 2.1481463579605694e-18im, 1.4037159244957018e-20 - 1.2489674484397012e-19im, 2.227823739527584e-20 - 4.839546058409615e-20im, -7.361455308883769e-20 - 2.900374681593027e-21im],
            atol = 1e-8)

        @test isapprox(
            10*log10.(abs2.(result.signal.S[result.signal.signalindex,result.signal.signalindex,:])),
            [-0.35881757066465897, -0.3851173962519461, -0.41414936742681174, -0.4462869667343253, -0.4819671123942377, -0.5217029870205777, -0.5660998488084007, -0.6158745949051934, -0.6718800569832, -0.7351352743488744, -0.8068633192233412, -0.888538646057728, -0.9819463922798333, -1.0892565316208547, -1.2131161713475298, -1.3567633723569377, -1.5241652233833656, -1.7201806973623992, -1.9507435138721474, -2.2230488432646025, -2.5457042687245037, -2.9287579585448764, -3.383421586955312, -3.921116498283568, -4.551109093859279, -5.27536187644465, -6.078360918578419, -6.909695380693948, -7.662866288954396, -8.174297241758815, -8.2876797612243, -7.968448353798562, -7.334037016046925, -6.5556629482650886, -5.764285202916524, -5.030263188923042, -4.380878193544146, -3.8202265335777112, -3.3419584829011564, -2.93599061269787, -2.591685034236857, -2.2992028878446895, -2.049972675586726, -1.836752988473692, -1.6535252962720193, -1.4953291435731257, -1.3580909299237387, -1.2384675791361006, -1.133712158214135, -1.0415620182221883, -0.9601473452230064],
            atol = 1e-8)

    end

    @testset verbose=true "hbnlsolve lossless" begin

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
        result=hbnlsolve(wp, Ip, Npumpmodes, circuit, circuitdefs,ports=[1])

        @test isapprox(
            result.nodeflux,
            ComplexF64[-0.013189736218636441 - 0.008651273221683266im, 2.6395682539011797e-5 - 5.6661514233313726e-6im, -2.3497663047131403e-8 - 2.2306503448330988e-8im, -7.990410788185173e-12 + 3.749499317166309e-11im, 4.336610521409153e-14 - 1.3388741572794773e-14im, -3.8059160268943944e-17 - 3.852723364098488e-17im, 6.07942296115915e-20 + 6.60227155809359e-21im, 5.977496795928632e-20 + 2.938144480595029e-20im, 0.12157241551965468 + 0.07973640247896553im, 1.3738923117026488e-5 - 6.46274778537762e-5im, -5.3393939365693146e-8 + 9.186223725369537e-9im, 2.790433171591809e-11 + 4.514438082272354e-11im, 3.3397076362502e-14 - 4.567840462363693e-14im, -6.153011440710407e-17 - 1.534143178629367e-17im, 6.419757193364206e-20 - 2.4735969038803813e-20im, 7.29011218720814e-20 + 2.6769912044348367e-21im],
            atol = 1e-8)

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

        @test_throws(
            DimensionMismatch("Number of currents Ip must be equal to number of pump ports"),
            hbnlsolve(wp, Ip, Npumpmodes, circuit, circuitdefs,ports=[1,2])
        )

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

        @test_warn(
            "Solver did not converge after maximum iterations of 1.",
            hbnlsolve(wp, Ip, Npumpmodes, circuit, circuitdefs,ports=[1],
                iterations=1)
        )

    end

    @testset "calcAoLjbm error" begin
        @variables Lj1 Lj2 A11 A12 A21 A22 A31 A32;
        @test_throws(
            DimensionMismatch("The length of Ljb cannot be larger than the number of branches."),
            JosephsonCircuits.calcAoLjbm([A11 A12;A21 A22;A31 A32],JosephsonCircuits.SparseArrays.sparsevec([1,2],[Lj1,Lj2]),1,2,1))
    end

   @testset "calcAoLjbm error" begin
        @variables Lj1 Lj2 A11 A12 A21 A22 A31 A32;
        @test_throws(
            DimensionMismatch("The second axis of Am must equal the number of nonzero elements in Ljb (the number of JJs)."),
            JosephsonCircuits.calcAoLjbm([A11;A21;A31],JosephsonCircuits.SparseArrays.sparsevec([1,2],[Lj1,Lj2]),1,2,2))
    end

    @testset "updateAoLjbm! error" begin
        @variables Lj1 Lj2 A11 A12 A21 A22 A31 A32;
        AoLjbm = JosephsonCircuits.calcAoLjbm([A11 A12;A21 A22;A31 A32],JosephsonCircuits.SparseArrays.sparsevec([1,2],[Lj1,Lj2]),1,2,2);
        AoLjbmcopy = copy(AoLjbm);
        AoLjbmcopy.nzval .= 0;

        @test_throws(
            DimensionMismatch("The number of nonzero elements in AoLjbm are not consistent with nnz(Ljb) and Nmodes."),
            JosephsonCircuits.updateAoLjbm!(AoLjbmcopy,[A11 A12;A21 A22;A31 A32],JosephsonCircuits.SparseArrays.sparsevec([1,2],[Lj1,Lj2]),1,3,2)
        )
    end

    @testset "updateAoLjbm! error" begin
        @variables Lj1 Lj2 A11 A12 A21 A22 A31 A32;
        AoLjbm = JosephsonCircuits.calcAoLjbm([A11 A12;A21 A22;A31 A32],JosephsonCircuits.SparseArrays.sparsevec([1,2],[Lj1,Lj2]),1,2,2);
        AoLjbmcopy = copy(AoLjbm);
        AoLjbmcopy.nzval .= 0;

        @test_throws(
            DimensionMismatch("The second axis of Am must equal the number of nonzero elements in Ljb (the number of JJs)."),
            JosephsonCircuits.updateAoLjbm!(AoLjbmcopy,[A11 A12 0;A21 A22 0;A31 A32 0],JosephsonCircuits.SparseArrays.sparsevec([1,2],[Lj1,Lj2]),1,2,2)
        )
    end

    @testset "updateAoLjbm! error" begin
        @variables Lj1 Lj2 A11 A12 A21 A22 A31 A32;
        AoLjbm = JosephsonCircuits.calcAoLjbm([A11 A12;A21 A22;A31 A32],JosephsonCircuits.SparseArrays.sparsevec([1,2],[Lj1,Lj2]),1,2,2);
        AoLjbmcopy = copy(AoLjbm);
        AoLjbmcopy.nzval .= 0;

        @test_throws(
            DimensionMismatch("The length of Ljb cannot be larger than the number of branches."),
            JosephsonCircuits.updateAoLjbm!(AoLjbmcopy,[A11 A12;A21 A22;A31 A32],JosephsonCircuits.SparseArrays.sparsevec([1,2],[Lj1,Lj2]),1,2,1)
        )
    end

    @testset "sincosnloddtoboth" begin
        @test_throws(
            DimensionMismatch("Length of node flux vector not consistent with number of modes number of frequencies"),
            JosephsonCircuits.sincosnloddtoboth([0.5+0.0im,0,0,0],1,5)
        )
    end

end