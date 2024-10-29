using JosephsonCircuits
using Test
using XicTools_jll

@testset verbose=true "spicewrapper" begin

    @testset "spice_hb_load" begin

        filepath = joinpath(dirname(Base.source_path()),"spicewrapper","invert_hb_lapack.cir.HB.FD.prn")

        out1 = JosephsonCircuits.spice_hb_load(filepath)

        out2 = (data = ComplexF64[0.993780083 + 3.99999905im 0.995269201 + 4.0000034im 0.999724201 + 4.00000209im 1.00106072 + 4.00000154im 1.00784458 + 4.0000023im 1.00456778 + 3.9999961im 1.00841795 + 4.00000033im 1.00289474 + 3.99999266im 0.998122333 + 3.99999862im 0.998782188 + 3.99999631im 0.98648007 + 3.99999862im 0.996688342 + 4.00000526im 0.988432616 + 3.99999987im 0.99814709 + 4.00001079im 1.0073127 + 4.00000099im 1.00101954 + 4.00000565im 1.02541765 + 4.0000011im 1.00237865 + 3.99999379im 1.01655617 + 4.00000043im 1.00135982 + 3.99999012im 0.975240926 + 3.99999978im 0.999355567 + 4.00000676im 0.937732855 + 3.99999973im 0.998205852 + 4.00003608im 0.971872324 + 4.00000011im 0.99865033 + 4.00004765im 1.13962865 + 4.00000034im 1.00019626 + 4.00000431im 1.39587828 + 4.00000011im 1.00164295 + 3.99928246im 3.00678901 + 4.0im 1.00164295 + 3.99928246im 1.39587828 + 4.00000011im 1.00019626 + 4.00000431im 1.13962865 + 4.00000034im 0.99865033 + 4.00004765im 0.971872324 + 4.00000011im 0.998205852 + 4.00003608im 0.937732855 + 3.99999973im 0.999355567 + 4.00000676im 0.975240926 + 3.99999978im 1.00135982 + 3.99999012im 1.01655617 + 4.00000043im 1.00237865 + 3.99999379im 1.02541765 + 4.0000011im 1.00101954 + 4.00000565im 1.0073127 + 4.00000099im 0.99814709 + 4.00001079im 0.988432616 + 3.99999987im 0.996688342 + 4.00000526im 0.98648007 + 3.99999862im 0.998782188 + 3.99999631im 0.998122333 + 3.99999862im 1.00289474 + 3.99999266im 1.00841795 + 4.00000033im 1.00456778 + 3.9999961im 1.00784458 + 4.0000023im 1.00106072 + 4.00000154im 0.999724201 + 4.00000209im 0.995269201 + 4.0000034im 0.993780083 + 3.99999905im; 3.99999905 + 4.0im 4.0000034 + 4.0im 4.00000209 + 4.0im 4.00000154 + 4.0im 4.0000023 + 4.0im 3.9999961 + 4.0im 4.00000033 + 4.0im 3.99999266 + 4.0im 3.99999862 + 4.0im 3.99999631 + 4.0im 3.99999862 + 4.0im 4.00000526 + 4.0im 3.99999987 + 4.0im 4.00001079 + 4.0im 4.00000099 + 4.0im 4.00000565 + 4.0im 4.0000011 + 4.0im 3.99999379 + 4.0im 4.00000043 + 4.0im 3.99999012 + 4.0im 3.99999978 + 4.0im 4.00000676 + 4.0im 3.99999973 + 4.0im 4.00003608 + 4.0im 4.00000011 + 4.0im 4.00004765 + 4.0im 4.00000034 + 4.0im 4.00000431 + 4.0im 4.00000011 + 4.0im 3.99928246 + 4.0im 4.0 + 4.0im 3.99928246 + 4.0im 4.00000011 + 4.0im 4.00000431 + 4.0im 4.00000034 + 4.0im 4.00004765 + 4.0im 4.00000011 + 4.0im 4.00003608 + 4.0im 3.99999973 + 4.0im 4.00000676 + 4.0im 3.99999978 + 4.0im 3.99999012 + 4.0im 4.00000043 + 4.0im 3.99999379 + 4.0im 4.0000011 + 4.0im 4.00000565 + 4.0im 4.00000099 + 4.0im 4.00001079 + 4.0im 3.99999987 + 4.0im 4.00000526 + 4.0im 3.99999862 + 4.0im 3.99999631 + 4.0im 3.99999862 + 4.0im 3.99999266 + 4.0im 4.00000033 + 4.0im 3.9999961 + 4.0im 4.0000023 + 4.0im 4.00000154 + 4.0im 4.00000209 + 4.0im 4.0000034 + 4.0im 3.99999905 + 4.0im; 4.0 + 1.0915949e-8im 4.0 - 9.36300661e-9im 4.0 + 6.48401424e-10im 4.0 + 3.60595283e-9im 4.0 - 7.35854719e-9im 4.0 + 1.34748831e-8im 4.0 - 5.8246132e-9im 4.0 + 8.31378765e-9im 4.0 + 1.14850831e-9im 4.0 - 7.61547701e-9im 4.0 + 5.12376887e-9im 4.0 - 1.70053674e-8im 4.0 + 2.89636083e-9im 4.0 - 7.79947704e-9im 4.0 - 1.61692365e-9im 4.0 + 1.22288148e-8im 4.0 - 3.21459597e-9im 4.0 + 2.12593418e-8im 4.0 - 1.10465795e-9im 4.0 + 5.77272169e-9im 4.0 + 1.45115638e-9im 4.0 - 2.25963228e-8im 4.0 + 1.70405046e-9im 4.0 - 3.43800868e-8im 4.0 + 2.2876213e-10im 4.0 - 9.28222926e-9im 4.0 - 7.21089193e-10im 4.0 + 4.22022329e-8im 4.0 - 3.44478005e-10im 4.0 + 4.50946372e-7im 4.0 + 0.0im 4.0 + 4.50946372e-7im 4.0 - 3.44478005e-10im 4.0 + 4.22022329e-8im 4.0 - 7.21089193e-10im 4.0 - 9.28222926e-9im 4.0 + 2.2876213e-10im 4.0 - 3.43800868e-8im 4.0 + 1.70405046e-9im 4.0 - 2.25963228e-8im 4.0 + 1.45115638e-9im 4.0 + 5.77272169e-9im 4.0 - 1.10465795e-9im 4.0 + 2.12593418e-8im 4.0 - 3.21459597e-9im 4.0 + 1.22288148e-8im 4.0 - 1.61692365e-9im 4.0 - 7.79947704e-9im 4.0 + 2.89636083e-9im 4.0 - 1.70053674e-8im 4.0 + 5.12376887e-9im 4.0 - 7.61547701e-9im 4.0 + 1.14850831e-9im 4.0 + 8.31378765e-9im 4.0 - 5.8246132e-9im 4.0 + 1.34748831e-8im 4.0 - 7.35854719e-9im 4.0 + 3.60595283e-9im 4.0 + 6.48401424e-10im 4.0 - 9.36300661e-9im 4.0 + 1.0915949e-8im; 1.0915949e-8 + 2.09869965e-9im -9.36300661e-9 - 8.79081664e-9im 6.48401424e-10 + 1.22823334e-8im 3.60595283e-9 - 8.04379553e-9im -7.35854719e-9 + 8.84849922e-9im 1.34748831e-8 + 4.50893608e-10im -5.8246132e-9 - 5.49538425e-9im 8.31378765e-9 + 6.17951173e-9im 1.14850831e-9 - 1.51048892e-8im -7.61547701e-9 + 4.16714387e-9im 5.12376887e-9 - 8.0641239e-9im -1.70053674e-8 - 1.50622812e-9im 2.89636083e-9 + 9.68283871e-9im -7.79947704e-9 - 4.14506739e-9im -1.61692365e-9 + 1.88277176e-8im 1.22288148e-8 - 1.908652e-9im -3.21459597e-9 + 6.86908443e-9im 2.12593418e-8 + 1.56560337e-9im -1.10465795e-9 - 1.6478723e-8im 5.77272169e-9 + 2.37874127e-9im 1.45115638e-9 - 2.61468149e-8im -2.25963228e-8 + 5.2964802e-10im 1.70405046e-9 - 6.53503815e-9im -3.43800868e-8 - 1.18458049e-9im 2.2876213e-10 + 2.91763346e-8im -9.28222926e-9 - 1.04704215e-9im -7.21089193e-10 + 4.44005409e-8im 4.22022329e-8 - 2.77413178e-11im -3.44478005e-10 - 1.05683898e-8im 4.50946372e-7 + 3.82019855e-10im 0.0 + 0.0im 4.50946372e-7 - 3.82019855e-10im -3.44478005e-10 + 1.05683898e-8im 4.22022329e-8 + 2.77413178e-11im -7.21089193e-10 - 4.44005409e-8im -9.28222926e-9 + 1.04704215e-9im 2.2876213e-10 - 2.91763346e-8im -3.43800868e-8 + 1.18458049e-9im 1.70405046e-9 + 6.53503815e-9im -2.25963228e-8 - 5.2964802e-10im 1.45115638e-9 + 2.61468149e-8im 5.77272169e-9 - 2.37874127e-9im -1.10465795e-9 + 1.6478723e-8im 2.12593418e-8 - 1.56560337e-9im -3.21459597e-9 - 6.86908443e-9im 1.22288148e-8 + 1.908652e-9im -1.61692365e-9 - 1.88277176e-8im -7.79947704e-9 + 4.14506739e-9im 2.89636083e-9 - 9.68283871e-9im -1.70053674e-8 + 1.50622812e-9im 5.12376887e-9 + 8.0641239e-9im -7.61547701e-9 - 4.16714387e-9im 1.14850831e-9 + 1.51048892e-8im 8.31378765e-9 - 6.17951173e-9im -5.8246132e-9 + 5.49538425e-9im 1.34748831e-8 - 4.50893608e-10im -7.35854719e-9 - 8.84849922e-9im 3.60595283e-9 + 8.04379553e-9im 6.48401424e-10 - 1.22823334e-8im -9.36300661e-9 + 8.79081664e-9im 1.0915949e-8 - 2.09869965e-9im; 2.09869965e-9 - 9.96332926e-9im -8.79081664e-9 + 5.96204411e-9im 1.22823334e-8 - 2.74151651e-9im -8.04379553e-9 - 5.14819657e-9im 8.84849922e-9 + 5.0630892e-9im 4.50893608e-10 - 9.57745551e-9im -5.49538425e-9 + 5.49305919e-9im 6.17951173e-9 - 9.76184383e-10im -1.51048892e-8 + 2.30349565e-10im 4.16714387e-9 + 1.13054389e-8im -8.0641239e-9 - 3.74318787e-9im -1.50622812e-9 + 1.17433684e-8im 9.68283871e-9 - 2.76879858e-9im -4.14506739e-9 - 2.99064776e-9im 1.88277176e-8 + 6.24941455e-10im -1.908652e-9 - 1.78779779e-8im 6.86908443e-9 + 2.11209346e-9im 1.56560337e-9 - 1.50452807e-8im -1.6478723e-8 + 6.79019581e-10im 2.37874127e-9 + 4.10948391e-9im -2.61468149e-8 - 1.23579912e-9im 5.2964802e-10 + 1.58398746e-8im -6.53503815e-9 - 1.43458514e-9im -1.18458049e-9 - 1.69746244e-9im 2.91763346e-8 - 3.34667433e-10im -1.04704215e-9 - 3.83675661e-8im 4.44005409e-8 + 3.76882627e-10im -2.77413178e-11 - 4.65154842e-8im -1.05683898e-8 + 2.31691194e-10im 3.82019855e-10 + 2.66589443e-7im 0.0 + 0.0im -3.82019855e-10 + 2.66589443e-7im 1.05683898e-8 + 2.31691194e-10im 2.77413178e-11 - 4.65154842e-8im -4.44005409e-8 + 3.76882627e-10im 1.04704215e-9 - 3.83675661e-8im -2.91763346e-8 - 3.34667433e-10im 1.18458049e-9 - 1.69746244e-9im 6.53503815e-9 - 1.43458514e-9im -5.2964802e-10 + 1.58398746e-8im 2.61468149e-8 - 1.23579912e-9im -2.37874127e-9 + 4.10948391e-9im 1.6478723e-8 + 6.79019581e-10im -1.56560337e-9 - 1.50452807e-8im -6.86908443e-9 + 2.11209346e-9im 1.908652e-9 - 1.78779779e-8im -1.88277176e-8 + 6.24941455e-10im 4.14506739e-9 - 2.99064776e-9im -9.68283871e-9 - 2.76879858e-9im 1.50622812e-9 + 1.17433684e-8im 8.0641239e-9 - 3.74318787e-9im -4.16714387e-9 + 1.13054389e-8im 1.51048892e-8 + 2.30349565e-10im -6.17951173e-9 - 9.76184383e-10im 5.49538425e-9 + 5.49305919e-9im -4.50893608e-10 - 9.57745551e-9im -8.84849922e-9 + 5.0630892e-9im 8.04379553e-9 - 5.14819657e-9im -1.22823334e-8 - 2.74151651e-9im 8.79081664e-9 + 5.96204411e-9im -2.09869965e-9 - 9.96332926e-9im; -9.96332926e-9 - 4.65417434e-9im 5.96204411e-9 + 5.26980111e-9im -2.74151651e-9 - 1.01345618e-8im -5.14819657e-9 + 7.13238073e-9im 5.0630892e-9 - 2.96900275e-9im -9.57745551e-9 + 1.22525931e-9im 5.49305919e-9 + 9.80046467e-9im -9.76184383e-10 - 4.10928741e-9im 2.30349565e-10 + 1.3607011e-8im 1.13054389e-8 - 3.41433382e-9im -3.74318787e-9 + 3.31301988e-9im 1.17433684e-8 + 1.00017305e-9im -2.76879858e-9 - 9.30250229e-9im -2.99064776e-9 + 3.5434984e-9im 6.24941455e-10 - 7.80518173e-9im -1.78779779e-8 + 2.08235278e-9im 2.11209346e-9 + 8.88129323e-9im -1.50452807e-8 - 7.54758338e-10im 6.79019581e-10 + 2.182614e-8im 4.10948391e-9 - 1.70423008e-9im -1.23579912e-9 + 1.14686435e-8im 1.58398746e-8 - 5.45018866e-10im -1.43458514e-9 - 1.73415323e-8im -1.69746244e-9 + 6.25830374e-10im -3.34667433e-10 - 3.36107146e-8im -3.83675661e-8 + 5.30820116e-10im 3.76882627e-10 - 1.1014799e-8im -4.65154842e-8 - 7.20824389e-11im 2.31691194e-10 + 1.67339046e-8im 2.66589443e-7 - 2.33049142e-11im 0.0 + 0.0im 2.66589443e-7 + 2.33049142e-11im 2.31691194e-10 - 1.67339046e-8im -4.65154842e-8 + 7.20824389e-11im 3.76882627e-10 + 1.1014799e-8im -3.83675661e-8 - 5.30820116e-10im -3.34667433e-10 + 3.36107146e-8im -1.69746244e-9 - 6.25830374e-10im -1.43458514e-9 + 1.73415323e-8im 1.58398746e-8 + 5.45018866e-10im -1.23579912e-9 - 1.14686435e-8im 4.10948391e-9 + 1.70423008e-9im 6.79019581e-10 - 2.182614e-8im -1.50452807e-8 + 7.54758338e-10im 2.11209346e-9 - 8.88129323e-9im -1.78779779e-8 - 2.08235278e-9im 6.24941455e-10 + 7.80518173e-9im -2.99064776e-9 - 3.5434984e-9im -2.76879858e-9 + 9.30250229e-9im 1.17433684e-8 - 1.00017305e-9im -3.74318787e-9 - 3.31301988e-9im 1.13054389e-8 + 3.41433382e-9im 2.30349565e-10 - 1.3607011e-8im -9.76184383e-10 + 4.10928741e-9im 5.49305919e-9 - 9.80046467e-9im -9.57745551e-9 - 1.22525931e-9im 5.0630892e-9 + 2.96900275e-9im -5.14819657e-9 - 7.13238073e-9im -2.74151651e-9 + 1.01345618e-8im 5.96204411e-9 - 5.26980111e-9im -9.96332926e-9 + 4.65417434e-9im], f = [-3.0e6, -2.9e6, -2.8e6, -2.7e6, -2.6e6, -2.5e6, -2.4e6, -2.3e6, -2.2e6, -2.1e6, -2.0e6, -1.9e6, -1.8e6, -1.7e6, -1.6e6, -1.5e6, -1.4e6, -1.3e6, -1.2e6, -1.1e6, -1.0e6, -900000.0, -800000.0, -700000.0, -600000.0, -500000.0, -400000.0, -300000.0, -200000.0, -100000.0, 0.0, 100000.0, 200000.0, 300000.0, 400000.0, 500000.0, 600000.0, 700000.0, 800000.0, 900000.0, 1.0e6, 1.1e6, 1.2e6, 1.3e6, 1.4e6, 1.5e6, 1.6e6, 1.7e6, 1.8e6, 1.9e6, 2.0e6, 2.1e6, 2.2e6, 2.3e6, 2.4e6, 2.5e6, 2.6e6, 2.7e6, 2.8e6, 2.9e6, 3.0e6], index = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0, 31.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0, 39.0, 40.0, 41.0, 42.0, 43.0, 44.0, 45.0, 46.0, 47.0, 48.0, 49.0, 50.0, 51.0, 52.0, 53.0, 54.0, 55.0, 56.0, 57.0, 58.0, 59.0, 60.0], header = SubString{String}["Index", "FREQ", "{V(VOUT)+1.0}", "{V(IN)+4.0}", "{V(1)+4.0}", "Re(IG(MP1))", "Im(IG(MP1))", "Re(IG(MN1))", "Im(IG(MN1))", "Re(IS(MP1))", "Im(IS(MP1))", "Re(ID(MN1))", "Im(ID(MN1))", "Re(ID(MP1))", "Im(ID(MP1))"]);

        @test all(out1.data .== out2.data)
        @test all(out1.f .== out2.f)
        @test all(out1.header .== out2.header)
        @test all(out1.index .== out2.index)

    end

    @testset "wrspice_input_transient" begin
        @testset "wrspice_input_transient errors" begin

            @test_throws(
                ArgumentError("Source nodes not strings or integers."),
                JosephsonCircuits.wrspice_input_transient("* SPICE Simulation",1e-6,5e9,3.14,(1.1,0),1e-9,100e-9,10e-9))

            @test_throws(
                ArgumentError("Input vector lengths not equal."),
                JosephsonCircuits.wrspice_input_transient("* SPICE Simulation",1e-6,5e9,3.14,(1,0,0),1e-9,100e-9,10e-9))

            @test_throws(
                ArgumentError("Input vector lengths not equal."),
                JosephsonCircuits.wrspice_input_transient("* SPICE Simulation",[1e-6,1e-3],[5e9,6e9],[3.14,6.28],[(1,0)],1e-9,100e-9,10e-9))

            @test_throws(
                ArgumentError("Two nodes are required per source."),
                JosephsonCircuits.wrspice_input_transient("* SPICE Simulation",[1e-6,1e-3],[5e9,6e9],[3.14,6.28],[(1,0),(1,0,2)],1e-9,100e-9,10e-9))

            @test_throws(
                ArgumentError("Nodes are not an integer or string."),
                JosephsonCircuits.wrspice_input_transient("* SPICE Simulation",[1e-6,1e-3],[5e9,6e9],[3.14,6.28],[(1,0),(1.1,0)],1e-9,100e-9,10e-9))

            # test various combinations of vctor and scalar inputs
            @test(JosephsonCircuits.wrspice_input_transient("* SPICE Simulation",1e-6,5e9,3.14,(1,0),1e-9,100e-9,10e-9) == JosephsonCircuits.wrspice_input_transient("* SPICE Simulation",[1e-6],[5e9],[3.14],(1,0),1e-9,100e-9,10e-9))
            @test(JosephsonCircuits.wrspice_input_transient("* SPICE Simulation",1e-6,5e9,3.14,(1,0),1e-9,100e-9,10e-9) == JosephsonCircuits.wrspice_input_transient("* SPICE Simulation",[1e-6],[5e9],[3.14],[(1,0)],1e-9,100e-9,10e-9))

        end


    end

    @testset "wrspice_input_ac array" begin
        out1 = JosephsonCircuits.wrspice_input_ac("* SPICE Simulation",collect((4:0.01:5)*1e9),[1,2],1e-6)
        out2 = "* SPICE Simulation\n* AC current source with magnitude 1 and phase 0\nisrc 1 0 ac 1.0e-6 0.0\n\n* Set up the AC small signal simulation\n.ac lin 99 4.0g 5.0g\n\n* The control block\n.control\n\n* Maximum size of data to export in kilobytes from 1e3 to 2e9 with\n* default 2.56e5. This has to come before the run command\nset maxdata=2.0e9\n\n* Run the simulation\nrun\n\n* Binary files are faster to save and load.\nset filetype=binary\n\n* Leave filename empty so we can add that as a command line argument.\n* Don't specify any variables so it saves everything.\nwrite\n\n.endc\n\n"
        @test out1 == out2
    end

    @testset "wrspice_input_ac float" begin
        out1 = JosephsonCircuits.wrspice_input_ac("* SPICE Simulation",4.0*1e9,[1,2],1e-6)
        out2 = "* SPICE Simulation\n* AC current source with magnitude 1 and phase 0\nisrc 1 0 ac 1.0e-6 0.0\n\n* Set up the AC small signal simulation\n.ac lin 1 4.0g 4.0g\n\n* The control block\n.control\n\n* Maximum size of data to export in kilobytes from 1e3 to 2e9 with\n* default 2.56e5. This has to come before the run command\nset maxdata=2.0e9\n\n* Run the simulation\nrun\n\n* Binary files are faster to save and load.\nset filetype=binary\n\n* Leave filename empty so we can add that as a command line argument.\n* Don't specify any variables so it saves everything.\nwrite\n\n.endc\n\n"
        @test out1 == out2
    end

    @testset "wrspice_input_ac float array" begin
        out1 = JosephsonCircuits.wrspice_input_ac("* SPICE Simulation",[4.0]*1e9,[1,2],1e-6)
        out2 = "* SPICE Simulation\n* AC current source with magnitude 1 and phase 0\nisrc 1 0 ac 1.0e-6 0.0\n\n* Set up the AC small signal simulation\n.ac lin 1 4.0g 4.0g\n\n* The control block\n.control\n\n* Maximum size of data to export in kilobytes from 1e3 to 2e9 with\n* default 2.56e5. This has to come before the run command\nset maxdata=2.0e9\n\n* Run the simulation\nrun\n\n* Binary files are faster to save and load.\nset filetype=binary\n\n* Leave filename empty so we can add that as a command line argument.\n* Don't specify any variables so it saves everything.\nwrite\n\n.endc\n\n"
        @test out1 == out2
    end

    if haskey(ENV,"CI")
        @testset "wrspice_cmd" begin
            @test_throws(
                ErrorException("WRSPICE executable not found. Please install WRSPICE or supply a path manually if already installed."),
                JosephsonCircuits.wrspice_cmd())
        end
    end

    # only run this test on Linux where WRspice is compiled by BinaryBuilder.jl
    if isdefined(XicTools_jll,:wrspice)
        @testset "XicTools.wrspice" begin
            input = "* SPICE Simulation\nR1 1 0 50.0\nC1 1 2 100.0f\nB1 2 0 3 jjk ics=0.32910597599999997u\nC2 2 0 674.18508376f\n.model jjk jj(rtype=0,cct=1,icrit=0.32910597599999997u,cap=325.81491624f,force=1,vm=9.9\n* Current source\n* 1-hyperbolic secant rise\nisrc 0 1 0.011300000000000001u*sin(29.84519304095611g*x+0.0)*(1-2/(exp(x/1.0e-8)+exp(-x/1.0e-8)))\nisrc2 0 1 0.0u*sin(29.84519304095611g*x+0.0)*(1-2/(exp(x/1.0e-8)+exp(-x/1.0e-8)))\nisrc3 0 1 0.0u*sin(0.0g*x+0.0)*(1-2/(exp(x/1.0e-8)+exp(-x/1.0e-8)))\n\n* Set up the transient simulation\n* .tran 5p 10n\n.tran 2.6315734072138794p 20.0n uic\n\n* The control block\n.control\nset maxdata=2.0e9\nset jjaccel=1\nset dphimax=0.01\nrun\nset filetype=binary\nwrite\n.endc\n\n"
            savedoutput = [1.9613450192807734e-24 2.0074538412437603e-16 2.864451080543789e-15 1.3120773989412112e-14 3.768124499190959e-14 8.399007993917563e-14 1.5954632153709466e-13 2.71557278630489e-13 4.266497702086937e-13 6.305226098969998e-13; 1.7830403984905592e-25 1.8245472126928553e-17 2.6017079896847483e-16 1.1903318606729214e-15 3.4127007259546297e-15 7.589659388906407e-15 1.4376315642898046e-14 2.4384847008116373e-14 3.815428549472124e-14 5.611610966214862e-14; 9.779370376757019e-24 3.016205159873438e-14 8.644788812965574e-13 6.021131461378511e-12 2.33445773551018e-11 6.582428233515808e-11 1.5170260251549065e-10 3.0433217417470633e-10 5.517474693786351e-10 9.256640408918309e-10]

            output1 = JosephsonCircuits.spice_run(input,XicTools_jll.wrspice())
            @test(
                isapprox(
                 output1.values["V"][:,1:10],
                 savedoutput,
                atol = 1e-6)
                )

            output2 = JosephsonCircuits.spice_run([input],XicTools_jll.wrspice())
            @test(
                isapprox(
                 output2[1].values["V"][:,1:10],
                 savedoutput,
                atol = 1e-6)
                )
        end
    end

end