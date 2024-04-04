using JosephsonCircuits
using Test

@testset verbose=true "spiceutils" begin

    @testset "wrspice_input_paramp" begin

        wp=2*pi*4.75001*1e9
        Ip=0.00565e-6
        wswrspice=2*pi*(4.5:0.5:5.0)*1e9
        sourcenodes = (1,0)

        netlist = "* SPICE Simulation\nR1 1 0 50\nC1 1 2 100.0f\nB1 2 0 3 jjk ics=0.32910597599999997u\nC2 2 0 670.8940240000001f\n.model jjk jj(rtype=0,cct=1,icrit=0.32910597599999997u,cap=329.105976f,force=1,vm=9.9"

        input1 =  JosephsonCircuits.wrspice_input_paramp(netlist, wswrspice, wp, 2*Ip, sourcenodes,sourcenodes)

        input2 = ["* SPICE Simulation\nR1 1 0 50\nC1 1 2 100.0f\nB1 2 0 3 jjk ics=0.32910597599999997u\nC2 2 0 670.8940240000001f\n.model jjk jj(rtype=0,cct=1,icrit=0.32910597599999997u,cap=329.105976f,force=1,vm=9.9\n* Current source\n* 1-hyperbolic secant rise\nisrc1 1 0 0.011300000000000001u*cos(29.84519304095611g*x+0.0)*(1-2/(exp(x/1.0e-8)+exp(-x/1.0e-8)))\n* Set up the transient simulation\n* .tran 5p 10n\n.tran 2.6315734072138794p 200.0n uic\n\n* The control block\n.control\nset maxdata=1.0e7\nset jjaccel=1\nset dphimax=0.01\nrun\nset filetype=binary\nwrite\n.endc\n\n", "* SPICE Simulation\nR1 1 0 50\nC1 1 2 100.0f\nB1 2 0 3 jjk ics=0.32910597599999997u\nC2 2 0 670.8940240000001f\n.model jjk jj(rtype=0,cct=1,icrit=0.32910597599999997u,cap=329.105976f,force=1,vm=9.9\n* Current source\n* 1-hyperbolic secant rise\nisrc1 1 0 0.011300000000000001u*cos(29.84519304095611g*x+0.0)*(1-2/(exp(x/1.0e-8)+exp(-x/1.0e-8)))\nisrc2 1 0 1.0000000000000001e-7u*cos(28.27433388230814g*x+-1.5707963267948966)*(1-2/(exp(x/1.0e-8)+exp(-x/1.0e-8)))\n* Set up the transient simulation\n* .tran 5p 10n\n.tran 2.6315734072138794p 200.0n uic\n\n* The control block\n.control\nset maxdata=1.0e7\nset jjaccel=1\nset dphimax=0.01\nrun\nset filetype=binary\nwrite\n.endc\n\n", "* SPICE Simulation\nR1 1 0 50\nC1 1 2 100.0f\nB1 2 0 3 jjk ics=0.32910597599999997u\nC2 2 0 670.8940240000001f\n.model jjk jj(rtype=0,cct=1,icrit=0.32910597599999997u,cap=329.105976f,force=1,vm=9.9\n* Current source\n* 1-hyperbolic secant rise\nisrc1 1 0 0.011300000000000001u*cos(29.84519304095611g*x+0.0)*(1-2/(exp(x/1.0e-8)+exp(-x/1.0e-8)))\nisrc2 1 0 1.0000000000000001e-7u*cos(28.27433388230814g*x+0.0)*(1-2/(exp(x/1.0e-8)+exp(-x/1.0e-8)))\n* Set up the transient simulation\n* .tran 5p 10n\n.tran 2.6315734072138794p 200.0n uic\n\n* The control block\n.control\nset maxdata=1.0e7\nset jjaccel=1\nset dphimax=0.01\nrun\nset filetype=binary\nwrite\n.endc\n\n", "* SPICE Simulation\nR1 1 0 50\nC1 1 2 100.0f\nB1 2 0 3 jjk ics=0.32910597599999997u\nC2 2 0 670.8940240000001f\n.model jjk jj(rtype=0,cct=1,icrit=0.32910597599999997u,cap=329.105976f,force=1,vm=9.9\n* Current source\n* 1-hyperbolic secant rise\nisrc1 1 0 0.011300000000000001u*cos(29.84519304095611g*x+0.0)*(1-2/(exp(x/1.0e-8)+exp(-x/1.0e-8)))\nisrc2 1 0 1.0000000000000001e-7u*cos(31.41592653589793g*x+-1.5707963267948966)*(1-2/(exp(x/1.0e-8)+exp(-x/1.0e-8)))\n* Set up the transient simulation\n* .tran 5p 10n\n.tran 2.6315734072138794p 200.0n uic\n\n* The control block\n.control\nset maxdata=1.0e7\nset jjaccel=1\nset dphimax=0.01\nrun\nset filetype=binary\nwrite\n.endc\n\n", "* SPICE Simulation\nR1 1 0 50\nC1 1 2 100.0f\nB1 2 0 3 jjk ics=0.32910597599999997u\nC2 2 0 670.8940240000001f\n.model jjk jj(rtype=0,cct=1,icrit=0.32910597599999997u,cap=329.105976f,force=1,vm=9.9\n* Current source\n* 1-hyperbolic secant rise\nisrc1 1 0 0.011300000000000001u*cos(29.84519304095611g*x+0.0)*(1-2/(exp(x/1.0e-8)+exp(-x/1.0e-8)))\nisrc2 1 0 1.0000000000000001e-7u*cos(31.41592653589793g*x+0.0)*(1-2/(exp(x/1.0e-8)+exp(-x/1.0e-8)))\n* Set up the transient simulation\n* .tran 5p 10n\n.tran 2.6315734072138794p 200.0n uic\n\n* The control block\n.control\nset maxdata=1.0e7\nset jjaccel=1\nset dphimax=0.01\nrun\nset filetype=binary\nwrite\n.endc\n\n"]

        @test input1 == input2
    end

end