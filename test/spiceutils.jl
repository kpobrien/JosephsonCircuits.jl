using JosephsonCircuits
using Test

@testset verbose=true "spiceutils" begin

    @testset "wrspice_input_paramp" begin

        wp=2*pi*4.75001*1e9
        Ip=0.00565e-6
        wswrspice=2*pi*(4.5:0.5:5.0)*1e9

        netlist = "* SPICE Simulation\nR1 1 0 50\nC1 1 2 100.0f\nB1 2 0 3 jjk ics=0.32910597599999997u\nC2 2 0 670.8940240000001f\n.model jjk jj(rtype=0,cct=1,icrit=0.32910597599999997u,cap=329.105976f,force=1,vm=9.9"

        input1 =  JosephsonCircuits.wrspice_input_paramp(netlist, wswrspice, wp, 2*Ip)

        input2 = ["* SPICE Simulation\nR1 1 0 50\nC1 1 2 100.0f\nB1 2 0 3 jjk ics=0.32910597599999997u\nC2 2 0 670.8940240000001f\n.model jjk jj(rtype=0,cct=1,icrit=0.32910597599999997u,cap=329.105976f,force=1,vm=9.9\n* Current source\n* 1-hyperbolic secant rise\nisrc 0 1 0.011300000000000001u*sin(29.84519304095611g*x+0.0)*(1-2/(exp(x/1.0e-8)+exp(-x/1.0e-8)))\nisrc2 0 1 0.0u*sin(0.0g*x+0.0)*(1-2/(exp(x/1.0e-8)+exp(-x/1.0e-8)))\n\n* isrc 1 0 sin(0 0.011300000000000001u 4.7500100000000005g 0.0 0.0)\n* isrc2 1 0 sin(0 0.0u 0.0g 0.0 0.0)\n\n* Set up the transient simulation\n* .tran 5p 10n\n.tran 2.6315734072138794p 200.0n uic\n\n* The control block\n.control\nset maxdata=1.0e7\nset jjaccel=1\nset dphimax=0.01\nrun\nset filetype=binary\nwrite\n.endc\n\n", "* SPICE Simulation\nR1 1 0 50\nC1 1 2 100.0f\nB1 2 0 3 jjk ics=0.32910597599999997u\nC2 2 0 670.8940240000001f\n.model jjk jj(rtype=0,cct=1,icrit=0.32910597599999997u,cap=329.105976f,force=1,vm=9.9\n* Current source\n* 1-hyperbolic secant rise\nisrc 0 1 0.011300000000000001u*sin(29.84519304095611g*x+0.0)*(1-2/(exp(x/1.0e-8)+exp(-x/1.0e-8)))\nisrc2 0 1 1.0000000000000001e-7u*sin(28.27433388230814g*x+0.0)*(1-2/(exp(x/1.0e-8)+exp(-x/1.0e-8)))\n\n* isrc 1 0 sin(0 0.011300000000000001u 4.7500100000000005g 0.0 0.0)\n* isrc2 1 0 sin(0 1.0000000000000001e-7u 4.5g 0.0 0.0)\n\n* Set up the transient simulation\n* .tran 5p 10n\n.tran 2.6315734072138794p 200.0n uic\n\n* The control block\n.control\nset maxdata=1.0e7\nset jjaccel=1\nset dphimax=0.01\nrun\nset filetype=binary\nwrite\n.endc\n\n", "* SPICE Simulation\nR1 1 0 50\nC1 1 2 100.0f\nB1 2 0 3 jjk ics=0.32910597599999997u\nC2 2 0 670.8940240000001f\n.model jjk jj(rtype=0,cct=1,icrit=0.32910597599999997u,cap=329.105976f,force=1,vm=9.9\n* Current source\n* 1-hyperbolic secant rise\nisrc 0 1 0.011300000000000001u*sin(29.84519304095611g*x+0.0)*(1-2/(exp(x/1.0e-8)+exp(-x/1.0e-8)))\nisrc2 0 1 1.0000000000000001e-7u*sin(28.27433388230814g*x+1.5707963267948966)*(1-2/(exp(x/1.0e-8)+exp(-x/1.0e-8)))\n\n* isrc 1 0 sin(0 0.011300000000000001u 4.7500100000000005g 0.0 0.0)\n* isrc2 1 0 sin(0 1.0000000000000001e-7u 4.5g 0.0 90.0)\n\n* Set up the transient simulation\n* .tran 5p 10n\n.tran 2.6315734072138794p 200.0n uic\n\n* The control block\n.control\nset maxdata=1.0e7\nset jjaccel=1\nset dphimax=0.01\nrun\nset filetype=binary\nwrite\n.endc\n\n", "* SPICE Simulation\nR1 1 0 50\nC1 1 2 100.0f\nB1 2 0 3 jjk ics=0.32910597599999997u\nC2 2 0 670.8940240000001f\n.model jjk jj(rtype=0,cct=1,icrit=0.32910597599999997u,cap=329.105976f,force=1,vm=9.9\n* Current source\n* 1-hyperbolic secant rise\nisrc 0 1 0.011300000000000001u*sin(29.84519304095611g*x+0.0)*(1-2/(exp(x/1.0e-8)+exp(-x/1.0e-8)))\nisrc2 0 1 1.0000000000000001e-7u*sin(31.41592653589793g*x+0.0)*(1-2/(exp(x/1.0e-8)+exp(-x/1.0e-8)))\n\n* isrc 1 0 sin(0 0.011300000000000001u 4.7500100000000005g 0.0 0.0)\n* isrc2 1 0 sin(0 1.0000000000000001e-7u 5.0g 0.0 0.0)\n\n* Set up the transient simulation\n* .tran 5p 10n\n.tran 2.6315734072138794p 200.0n uic\n\n* The control block\n.control\nset maxdata=1.0e7\nset jjaccel=1\nset dphimax=0.01\nrun\nset filetype=binary\nwrite\n.endc\n\n", "* SPICE Simulation\nR1 1 0 50\nC1 1 2 100.0f\nB1 2 0 3 jjk ics=0.32910597599999997u\nC2 2 0 670.8940240000001f\n.model jjk jj(rtype=0,cct=1,icrit=0.32910597599999997u,cap=329.105976f,force=1,vm=9.9\n* Current source\n* 1-hyperbolic secant rise\nisrc 0 1 0.011300000000000001u*sin(29.84519304095611g*x+0.0)*(1-2/(exp(x/1.0e-8)+exp(-x/1.0e-8)))\nisrc2 0 1 1.0000000000000001e-7u*sin(31.41592653589793g*x+1.5707963267948966)*(1-2/(exp(x/1.0e-8)+exp(-x/1.0e-8)))\n\n* isrc 1 0 sin(0 0.011300000000000001u 4.7500100000000005g 0.0 0.0)\n* isrc2 1 0 sin(0 1.0000000000000001e-7u 5.0g 0.0 90.0)\n\n* Set up the transient simulation\n* .tran 5p 10n\n.tran 2.6315734072138794p 200.0n uic\n\n* The control block\n.control\nset maxdata=1.0e7\nset jjaccel=1\nset dphimax=0.01\nrun\nset filetype=binary\nwrite\n.endc\n\n"]

        @test input1 == input2
    end

    # @testset "calcgwrspice" begin



    #     # filepath = joinpath(dirname(Base.source_path()),"spiceraw","test01.raw")

    #     # out1 = JosephsonCircuits.spice_raw_load(filepath)

    #     # out2 = (variables = Dict("V" => ["v(1)", "v(2)", "v(3)"], "Hz" => ["frequency"]), values = Dict{Any, Any}("V" => ComplexF64[48.87562301047733 - 7.413126995337487im 49.97131616467212 + 1.1949290155299537im 49.02611690128596 - 6.90980805243651im; -10.116167243319213 + 1.534380793728424im 57.578470543293086 + 1.3775359827006193im 12.368446655904192 - 1.743197747303436im; 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im], "Hz" => ComplexF64[4.0e9 + 0.0im 5.0e9 + 0.0im 6.0e9 + 0.0im]))

    #     # @test out1.variables == out2.variables
    #     # @test out1.values == out2.values
    # end


end