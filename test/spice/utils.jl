using JosephsonCircuits
using Test

@testset verbose=true "spiceutils" begin

    @testset "wrspice_input_paramp" begin

        wp=2*pi*4.75001*1e9
        Ip=0.00565e-6
        wswrspice=2*pi*(4.5:0.5:5.0)*1e9
        # the drives into port node 1, from ground
        sourcenodes = (0,1)

        netlist = "* SPICE Simulation\nR1 1 0 50\nC1 1 2 100.0f\nB1 2 0 3 jjk ics=0.32910597599999997u\nC2 2 0 670.8940240000001f\n.model jjk jj(rtype=0,cct=1,icrit=0.32910597599999997u,cap=329.105976f,force=1,vm=9.9"

        input1 =  JosephsonCircuits.wrspice_input_paramp(netlist, wswrspice, wp, 2*Ip, sourcenodes,sourcenodes)

        input2 = ["* SPICE Simulation\nR1 1 0 50\nC1 1 2 100.0f\nB1 2 0 3 jjk ics=0.32910597599999997u\nC2 2 0 670.8940240000001f\n.model jjk jj(rtype=0,cct=1,icrit=0.32910597599999997u,cap=329.105976f,force=1,vm=9.9\n* Current source\n* 1-hyperbolic secant rise\nisrc1 0 1 0.011300000000000001u*cos(29.845193040956104g*x+0.0)*(1-2/(exp(x/1.0e-8)+exp(-x/1.0e-8)))\n* Set up the transient simulation\n* .tran 5p 10n\n.tran 2.6315734072138794p 200.0n uic\n\n* The control block\n.control\nset maxdata=2.0e9\nset jjaccel=1\nset dphimax=0.01\nrun\nset filetype=binary\nwrite\n.endc\n\n", "* SPICE Simulation\nR1 1 0 50\nC1 1 2 100.0f\nB1 2 0 3 jjk ics=0.32910597599999997u\nC2 2 0 670.8940240000001f\n.model jjk jj(rtype=0,cct=1,icrit=0.32910597599999997u,cap=329.105976f,force=1,vm=9.9\n* Current source\n* 1-hyperbolic secant rise\nisrc1 0 1 0.011300000000000001u*cos(29.845193040956104g*x+0.0)*(1-2/(exp(x/1.0e-8)+exp(-x/1.0e-8)))\nisrc2 0 1 1.0000000000000001e-7u*cos(28.27433388230814g*x+-1.5707963267948966)*(1-2/(exp(x/1.0e-8)+exp(-x/1.0e-8)))\n* Set up the transient simulation\n* .tran 5p 10n\n.tran 2.6315734072138794p 200.0n uic\n\n* The control block\n.control\nset maxdata=2.0e9\nset jjaccel=1\nset dphimax=0.01\nrun\nset filetype=binary\nwrite\n.endc\n\n", "* SPICE Simulation\nR1 1 0 50\nC1 1 2 100.0f\nB1 2 0 3 jjk ics=0.32910597599999997u\nC2 2 0 670.8940240000001f\n.model jjk jj(rtype=0,cct=1,icrit=0.32910597599999997u,cap=329.105976f,force=1,vm=9.9\n* Current source\n* 1-hyperbolic secant rise\nisrc1 0 1 0.011300000000000001u*cos(29.845193040956104g*x+0.0)*(1-2/(exp(x/1.0e-8)+exp(-x/1.0e-8)))\nisrc2 0 1 1.0000000000000001e-7u*cos(28.27433388230814g*x+0.0)*(1-2/(exp(x/1.0e-8)+exp(-x/1.0e-8)))\n* Set up the transient simulation\n* .tran 5p 10n\n.tran 2.6315734072138794p 200.0n uic\n\n* The control block\n.control\nset maxdata=2.0e9\nset jjaccel=1\nset dphimax=0.01\nrun\nset filetype=binary\nwrite\n.endc\n\n", "* SPICE Simulation\nR1 1 0 50\nC1 1 2 100.0f\nB1 2 0 3 jjk ics=0.32910597599999997u\nC2 2 0 670.8940240000001f\n.model jjk jj(rtype=0,cct=1,icrit=0.32910597599999997u,cap=329.105976f,force=1,vm=9.9\n* Current source\n* 1-hyperbolic secant rise\nisrc1 0 1 0.011300000000000001u*cos(29.845193040956104g*x+0.0)*(1-2/(exp(x/1.0e-8)+exp(-x/1.0e-8)))\nisrc2 0 1 1.0000000000000001e-7u*cos(31.41592653589793g*x+-1.5707963267948966)*(1-2/(exp(x/1.0e-8)+exp(-x/1.0e-8)))\n* Set up the transient simulation\n* .tran 5p 10n\n.tran 2.6315734072138794p 200.0n uic\n\n* The control block\n.control\nset maxdata=2.0e9\nset jjaccel=1\nset dphimax=0.01\nrun\nset filetype=binary\nwrite\n.endc\n\n", "* SPICE Simulation\nR1 1 0 50\nC1 1 2 100.0f\nB1 2 0 3 jjk ics=0.32910597599999997u\nC2 2 0 670.8940240000001f\n.model jjk jj(rtype=0,cct=1,icrit=0.32910597599999997u,cap=329.105976f,force=1,vm=9.9\n* Current source\n* 1-hyperbolic secant rise\nisrc1 0 1 0.011300000000000001u*cos(29.845193040956104g*x+0.0)*(1-2/(exp(x/1.0e-8)+exp(-x/1.0e-8)))\nisrc2 0 1 1.0000000000000001e-7u*cos(31.41592653589793g*x+0.0)*(1-2/(exp(x/1.0e-8)+exp(-x/1.0e-8)))\n* Set up the transient simulation\n* .tran 5p 10n\n.tran 2.6315734072138794p 200.0n uic\n\n* The control block\n.control\nset maxdata=2.0e9\nset jjaccel=1\nset dphimax=0.01\nrun\nset filetype=binary\nwrite\n.endc\n\n"]

        @test input1 == input2
    end

    @testset "wrspice_calcS_paramp" begin

        Nnodes = 2
        ws = 2*pi*5e9
        wp = 2*pi*6e9
        stepsperperiod = 80
        t = LinRange(0,2*pi/wp,stepsperperiod)
        Is = 1e-13
        Vpump = zeros(1,length(t))
        Vsignalsin = zeros(1,length(t))
        Vsignalcos = zeros(1,length(t))
        Vpump[1,:] .= sin.(2*pi*wp*t)
        Vsignalsin[1,:] .= sin.(2*pi*wp*t)+Is/50*sin.(2*pi*ws*t)
        Vsignalcos[1,:] .= sin.(2*pi*wp*t)+Is/50*cos.(2*pi*ws*t)

        out = [(values=Dict("S"=>t,"V"=>Vpump),),
                (values=Dict("S"=>t,"V"=>Vsignalsin),),
                (values=Dict("S"=>t,"V"=>Vsignalcos),),
                ];
        out[1].values["V"];
        out[2].values["V"];
        out[3].values["V"];

        @test_throws(
            ErrorException("Number of WRspice simulations not consistent with number of frequencies."),
            JosephsonCircuits.wrspice_calcS_paramp(out, [2*pi,2*pi*2], Nnodes;stepsperperiod = stepsperperiod, Is = Is),
            )

        # the port traces are v(1) and v(Nnodes-1), found by name: a JPA
        # whose nets are named by words has neither, and the rows of its
        # rawfile, sorted numbered nodes first, start with the junction's
        # phase node v(3)
        V = zeros(3, length(t))
        worded = [(variables = Dict("V" => ["v(3)", "v(in)", "v(mid)"]),
            values = Dict("S" => t, "V" => V)) for _ in 1:3]
        @test_throws ArgumentError JosephsonCircuits.wrspice_calcS_paramp(
            worded, [ws], 3; stepsperperiod = stepsperperiod, Is = Is)

    end
end