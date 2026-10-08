using JosephsonCircuits
using Test
using XicTools_jll

@testset verbose=true "spicewrapper" begin

    @testset "wrspice_input_transient" begin
        @testset "wrspice_input_transient errors" begin

            @test_throws(
                ArgumentError("Source nodes not strings or integers."),
                JosephsonCircuits.wrspice_input_transient("* SPICE Simulation",1e-6,2pi*5e9,3.14,(1.1,0),1e-9,100e-9,10e-9))

            @test_throws(
                ArgumentError("Input vector lengths not equal."),
                JosephsonCircuits.wrspice_input_transient("* SPICE Simulation",1e-6,2pi*5e9,3.14,(1,0,0),1e-9,100e-9,10e-9))

            @test_throws(
                ArgumentError("Input vector lengths not equal."),
                JosephsonCircuits.wrspice_input_transient("* SPICE Simulation",[1e-6,1e-3],2pi*[5e9,6e9],[3.14,6.28],[(1,0)],1e-9,100e-9,10e-9))

            @test_throws(
                ArgumentError("Two nodes are required per source."),
                JosephsonCircuits.wrspice_input_transient("* SPICE Simulation",[1e-6,1e-3],2pi*[5e9,6e9],[3.14,6.28],[(1,0),(1,0,2)],1e-9,100e-9,10e-9))

            @test_throws(
                ArgumentError("Nodes are not an integer or string."),
                JosephsonCircuits.wrspice_input_transient("* SPICE Simulation",[1e-6,1e-3],2pi*[5e9,6e9],[3.14,6.28],[(1,0),(1.1,0)],1e-9,100e-9,10e-9))

            # test various combinations of vctor and scalar inputs
            @test(JosephsonCircuits.wrspice_input_transient("* SPICE Simulation",1e-6,2pi*5e9,3.14,(1,0),1e-9,100e-9,10e-9) == JosephsonCircuits.wrspice_input_transient("* SPICE Simulation",[1e-6],[2pi*5e9],[3.14],(1,0),1e-9,100e-9,10e-9))
            @test(JosephsonCircuits.wrspice_input_transient("* SPICE Simulation",1e-6,2pi*5e9,3.14,(1,0),1e-9,100e-9,10e-9) == JosephsonCircuits.wrspice_input_transient("* SPICE Simulation",[1e-6],[2pi*5e9],[3.14],[(1,0)],1e-9,100e-9,10e-9))

        end

        @testset "the angular frequency" begin
            # a source is written as SPICE's cos(w*t), with the angular
            # frequency it is given, in units of 1e9 radians per second
            w = 2pi*5e9
            input = JosephsonCircuits.wrspice_input_transient("* SPICE Simulation",
                1e-6, w, 0.0, (0, 1), 1e-12, 1e-9, 1e-10)
            @test occursin("cos($(w*1e-9)g*x", input)
        end


    end

    @testset "wrspice_input_ac array" begin
        out1 = JosephsonCircuits.wrspice_input_ac("* SPICE Simulation",collect(2pi*(4:0.01:5)*1e9),[1,2],1e-6)
        out2 = "* SPICE Simulation\n* AC current source into the port\nisrc 0 1 ac 1.0e-6 0.0\n\n* Set up the AC small signal simulation\n.ac lin 99 4.0g 5.0g\n\n* The control block\n.control\n\n* Maximum size of data to export in kilobytes from 1e3 to 2e9 with\n* default 2.56e5. This has to come before the run command\nset maxdata=2.0e9\n\n* Run the simulation\nrun\n\n* Binary files are faster to save and load.\nset filetype=binary\n\n* Leave filename empty so we can add that as a command line argument.\n* Don't specify any variables so it saves everything.\nwrite\n\n.endc\n\n"
        @test out1 == out2
    end

    @testset "wrspice_input_ac float" begin
        out1 = JosephsonCircuits.wrspice_input_ac("* SPICE Simulation",2pi*4.0e9,[1,2],1e-6)
        out2 = "* SPICE Simulation\n* AC current source into the port\nisrc 0 1 ac 1.0e-6 0.0\n\n* Set up the AC small signal simulation\n.ac lin 1 4.0g 4.0g\n\n* The control block\n.control\n\n* Maximum size of data to export in kilobytes from 1e3 to 2e9 with\n* default 2.56e5. This has to come before the run command\nset maxdata=2.0e9\n\n* Run the simulation\nrun\n\n* Binary files are faster to save and load.\nset filetype=binary\n\n* Leave filename empty so we can add that as a command line argument.\n* Don't specify any variables so it saves everything.\nwrite\n\n.endc\n\n"
        @test out1 == out2
    end

    @testset "wrspice_input_ac float array" begin
        out1 = JosephsonCircuits.wrspice_input_ac("* SPICE Simulation",2pi*[4.0e9],[1,2],1e-6)
        out2 = "* SPICE Simulation\n* AC current source into the port\nisrc 0 1 ac 1.0e-6 0.0\n\n* Set up the AC small signal simulation\n.ac lin 1 4.0g 4.0g\n\n* The control block\n.control\n\n* Maximum size of data to export in kilobytes from 1e3 to 2e9 with\n* default 2.56e5. This has to come before the run command\nset maxdata=2.0e9\n\n* Run the simulation\nrun\n\n* Binary files are faster to save and load.\nset filetype=binary\n\n* Leave filename empty so we can add that as a command line argument.\n* Don't specify any variables so it saves everything.\nwrite\n\n.endc\n\n"
        @test out1 == out2
    end

    @testset "wrspice_input_ac points" begin
        # one frequency from a range is written as from a vector, and two
        # are refused, since WRSPICE answers a sweep of no intermediate
        # steps with three points
        ac(ws) = JosephsonCircuits.wrspice_input_ac("* SPICE Simulation", ws, [1, 2], 1e-6)
        @test ac(2pi*(4:1:4)*1e9) == ac(2pi*[4.0e9]) == ac(2pi*4.0e9)
        @test_throws ArgumentError ac(2pi*[4.0e9, 5.0e9])
        @test_throws ArgumentError ac(2pi*(4:1:5)*1e9)
        @test_throws ArgumentError ac(Float64[])
        # WRSPICE answers a linear sweep at equal steps from the first
        # frequency to the last, so frequencies which are not uniformly
        # spaced are refused
        @test_throws ArgumentError ac(2pi*[4.0e9, 4.1e9, 5.0e9])
    end

    @testset "wrspice_cmd" begin
        # The refusal when there is nothing to find: no binary in the JLL
        # for this platform and no installation at WRSPICE's standard
        # path. Those are the two places `wrspice_cmd` looks, so a
        # machine with either cannot test the refusal. The condition is
        # on those rather than on the continuous integration environment,
        # so that a developer on a platform the JLL does not cover runs
        # it too.
        standard = JosephsonCircuits.wrspicestandardpath()
        if !isdefined(XicTools_jll, :wrspice) &&
                !isfile(standard) && !islink(standard)
            @test_throws(
                ErrorException("WRSPICE executable not found. Please install WRSPICE, load XicTools_jll, or supply a path manually if installed elsewhere."),
                JosephsonCircuits.wrspice_cmd())
        end
    end

    # only run this test on Linux where WRspice is compiled by BinaryBuilder.jl
    if isdefined(XicTools_jll,:wrspice)
        @testset "XicTools.wrspice" begin
            input = "* SPICE Simulation\nR1 1 0 50.0\nC1 1 2 100.0f\nB1 2 0 3 jjk ics=0.32910597599999997u\nC2 2 0 674.18508376f\n.model jjk jj(rtype=0,cct=1,icrit=0.32910597599999997u,cap=325.81491624f,force=1,vm=9.9\n* Current source\n* 1-hyperbolic secant rise\nisrc 0 1 0.011300000000000001u*sin(29.84519304095611g*x+0.0)*(1-2/(exp(x/1.0e-8)+exp(-x/1.0e-8)))\nisrc2 0 1 0.0u*sin(29.84519304095611g*x+0.0)*(1-2/(exp(x/1.0e-8)+exp(-x/1.0e-8)))\nisrc3 0 1 0.0u*sin(0.0g*x+0.0)*(1-2/(exp(x/1.0e-8)+exp(-x/1.0e-8)))\n\n* Set up the transient simulation\n* .tran 5p 10n\n.tran 2.6315734072138794p 20.0n uic\n\n* The control block\n.control\nset maxdata=2.0e9\nset jjaccel=1\nset dphimax=0.01\nrun\nset filetype=binary\nwrite\n.endc\n\n"
            savedoutput = [1.9613450192807734e-24 2.0074538412437603e-16 2.864451080543789e-15 1.3120773989412112e-14 3.768124499190959e-14 8.399007993917563e-14 1.5954632153709466e-13 2.71557278630489e-13 4.266497702086937e-13 6.305226098969998e-13; 1.7830403984905592e-25 1.8245472126928553e-17 2.6017079896847483e-16 1.1903318606729214e-15 3.4127007259546297e-15 7.589659388906407e-15 1.4376315642898046e-14 2.4384847008116373e-14 3.815428549472124e-14 5.611610966214862e-14; 9.779370376757019e-24 3.016205159873438e-14 8.644788812965574e-13 6.021131461378511e-12 2.33445773551018e-11 6.582428233515808e-11 1.5170260251549065e-10 3.0433217417470633e-10 5.517474693786351e-10 9.256640408918309e-10]

            # the saved node voltages are of order 1e-13 and below, so they
            # are compared relative to their size, which a drive of the
            # other sign or a wrong value fails
            output1 = JosephsonCircuits.spice_run(input,JosephsonCircuits.wrspice_cmd())
            @test(
                isapprox(
                 output1.values["V"][:,1:10],
                 savedoutput,
                rtol = 1e-6)
                )

            output2 = JosephsonCircuits.spice_run([input],JosephsonCircuits.wrspice_cmd())
            @test(
                isapprox(
                 output2[1].values["V"][:,1:10],
                 savedoutput,
                rtol = 1e-6)
                )

            # the AC source injects its current into the port with its
            # phase: the node voltage of a resistor and a capacitor in
            # parallel is the current times their impedance
            rc = "* RC\nR1 1 0 50\nC1 1 0 1p"
            for I in (1e-6, 1e-6*cis(pi/2), -1e-6)
                ac = JosephsonCircuits.spice_run(
                    JosephsonCircuits.wrspice_input_ac(rc, 2pi*5.0e9, [1, 2], I),
                    JosephsonCircuits.wrspice_cmd())
                @test ac.values["V"][1, 1] ≈ I/(1/50 + im*2pi*5e9*1e-12) rtol = 1e-6
            end
            # a sweep is answered at the frequencies asked for, which
            # the rawfile gives in Hz
            for ws in (2pi*[4e9, 4.5e9, 5e9], 2pi*(4:0.1:5)*1e9)
                ac = JosephsonCircuits.spice_run(
                    JosephsonCircuits.wrspice_input_ac(rc, ws, [1, 2], 1e-6),
                    JosephsonCircuits.wrspice_cmd())
                @test vec(real.(ac.values["Hz"])) ≈ ws/(2pi) rtol = 1e-12
            end

            # an element of a model the input does not define: WRSPICE
            # exits normally and writes the plot of its constants in place
            # of the analysis, which is an error
            bad = "* an element with an unknown model\nB1 1 0 2 nomodel ics=1u\nR1 1 0 50\n.tran 1p 10p uic\n.control\nrun\nset filetype=binary\nwrite\n.endc\n"
            @test_throws ErrorException JosephsonCircuits.spice_run(bad,
                JosephsonCircuits.wrspice_cmd())
        end
    end

end