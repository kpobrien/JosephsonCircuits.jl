using JosephsonCircuits
using Test
import UUIDs

@testset verbose=true "touchstone" begin

    @testset "touchstone_parse" begin

        examples = [
            "!Example 1:\n!1-port S-parameter file, single frequency point\n# MHz S MA R 50\n!freq magS11 angS11\n2.000 0.894 -12.136",
            "!Example 1a:\n!1-port S-parameter file, single frequency point\n# MHz S DB R 50\n!freq magS11 angS11\n2.000 -0.97 -12.136",
            "!Example 2:\n!1-port Z-parameter file, multiple frequency points\n# MHz Z MA R 75\n!freq magZ11 angZ11\n100 0.99 -4\n200 0.80 -22\n300 0.707 -45\n400 0.40 -62\n500 0.01 -89",
            "!Example 3:\n!2-port H-parameter file, single frequency point\n# KHz H MA R 1\n! freq magH11 angH11 magH21 angH21 magH12 angH12 magH22 angH22\n2 .95 -26 3.57 157 .04 76 .66 -1",
            "!Example 4:\n!2-port S-parameter file, three frequency points\n# GHZ S RI R 50.0\n!freq RelS11 ImS11 ReS21 ImS21 ReS12 ImS12 ReS22 ImS22\n1.0000 0.3926 -0.1211 -0.0003 -0.0021 -0.0003 -0.0021 0.3926 -0.1211\n2.0000 0.3517 -0.3054 -0.0096 -0.0298 -0.0096 -0.0298 0.3517 -0.3054\n10.000 0.3419 0.3336 -0.0134 0.0379 -0.0134 0.0379 0.3419 0.3336",
            "!Example 5:\n! 4-port S-parameter data, taken at three frequency points\n# GHZ S MA R 50\n5.00000 0.60 161.24 0.40 -42.20 0.42 -66.58 0.53 -79.34 !row 1\n0.40 -42.20 0.60 161.20 0.53 -79.34 0.42 -66.58 !row 2\n0.42 -66.58 0.53 -79.34 0.60 161.24 0.40 -42.20 !row 3\n0.53 -79.34 0.42 -66.58 0.40 -42.20 0.60 161.24 !row 4\n6.00000 0.57 150.37 0.40 -44.34 0.41 -81.24 0.57 -95.77 !row 1\n0.40 -44.34 0.57 150.37 0.57 -95.77 0.41 -81.24 !row 2\n0.41 -81.24 0.57 -95.77 0.57 150.37 0.40 -44.34 !row 3\n0.57 -95.77 0.41 -81.24 0.40 -44.34 0.57 150.37 !row 4\n7.00000 0.50 136.69 0.45 -46.41 0.37 -99.09 0.62 -114.19 !row 1\n0.45 -46.41 0.50 136.69 0.62 -114.19 0.37 -99.09 !row 2\n0.37 -99.09 0.62 -114.19 0.50 136.69 0.45 -46.41 !row 3\n0.62 -114.19 0.37 -99.09 0.45 -46.41 0.50 136.69 !row 4",
            "!Example 8:\n!2-port network, S-parameter and noise data\n# GHZ S MA R 50\n2 .95 -26 3.57 157 .04 76 .66 -14\n22 .60 -144 1.30 40 .14 40 .56 -85\n! NOISE PARAMETERS\n4 .7 .64 69 .38\n18 2.7 .46 -33 .40",
            "!Example 4 (Version 2.0):\n! 4-port S-parameter data\n! Default impedance is overridden by [Reference]\n! Data cannot be represented using 1.0 syntax\n! Note that the [Reference] keyword arguments appear on a separate line\n[Version] 2.0\n# GHz S MA R 50\n[Number of Ports] 4\n[Reference]\n50 75 0.01 0.01\n[Number of Frequencies] 1\n[Network Data]\n5.00000 0.60 161.24 0.40 -42.20 0.42 -66.58 0.53 -79.34 !row 1\n0.40 -42.20 0.60 161.20 0.53 -79.34 0.42 -66.58 !row 2\n0.42 -66.58 0.53 -79.34 0.60 161.24 0.40 -42.20 !row 3\n0.53 -79.34 0.42 -66.58 0.40 -42.20 0.60 161.24 !row 4",
            "!Example 5 (Version 2.0):\n! 4-port S-parameter data\n! Default impedance is overridden by the [Reference] keyword arguments\n! Data cannot be represented using 1.0 syntax\n[Version] 2.0\n# GHz S MA R 50\n[Number of Ports] 4\n[Number of Frequencies] 1\n[Reference] 50 75 0.01 0.01\n[Matrix Format] Full\n[Network Data]\n5.00000 0.60 161.24 0.40 -42.20 0.42 -66.58 0.53 -79.34 !row 1\n0.40 -42.20 0.60 161.20 0.53 -79.34 0.42 -66.58 !row 2\n0.42 -66.58 0.53 -79.34 0.60 161.24 0.40 -42.20 !row 3\n0.53 -79.34 0.42 -66.58 0.40 -42.20 0.60 161.24 !row 4\n[End]",
            # same as above but no matrix format
            "!Example 5 (Version 2.0):\n! 4-port S-parameter data\n! Default impedance is overridden by the [Reference] keyword arguments\n! Data cannot be represented using 1.0 syntax\n[Version] 2.0\n# GHz S MA R 50\n[Number of Ports] 4\n[Number of Frequencies] 1\n[Reference] 50 75 0.01 0.01\n[Network Data]\n5.00000 0.60 161.24 0.40 -42.20 0.42 -66.58 0.53 -79.34 !row 1\n0.40 -42.20 0.60 161.20 0.53 -79.34 0.42 -66.58 !row 2\n0.42 -66.58 0.53 -79.34 0.60 161.24 0.40 -42.20 !row 3\n0.53 -79.34 0.42 -66.58 0.40 -42.20 0.60 161.24 !row 4\n[End]",
            "!Example 6 (Version 2.0):\n! 4-port S-parameter data\n! Default impedance is overridden by the [Reference] keyword arguments\n! Note that [Reference] arguments are split across two lines\n! Data cannot be represented using 1.0 syntax\n[Version] 2.0\n# GHz S MA R 50\n[Number of Ports] 4\n[Number of Frequencies] 1\n[Reference] 50 75\n0.01 0.01\n[Matrix Format] Lower\n[Network Data]\n5.00000 0.60 161.24 !row 1\n0.40 -42.20 0.60 161.20 !row 2\n0.42 -66.58 0.53 -79.34 0.60 161.24 !row 3\n0.53 -79.34 0.42 -66.58 0.40 -42.20 0.60 161.24 !row 4\n[End]",
            "!Example 7 (Version 2.0):\n!1-port Z-parameter file, multiple frequency points\n[Version] 2.0\n# MHz Z MA\n[Number of Ports] 1\n[Number of Frequencies] 5\n[Reference] 20.0\n[Network Data]\n!freq magZ11 angZ11\n100 74.25 -4\n200 60 -22\n300 53.025 -45\n400 30 -62\n500 0.75 -89\n[End]",
            "!Example 8 (Version 1.0):\n!1-port S-parameter file, single frequency point\n# MHz S MA R 50\n!freq magS11 angS11\n2.000 0.894 -12.136",
            "!Example 9 (Version 1.0):\n!1-port Z-parameter file, multiple frequency points\n# MHz Z MA R 75\n!freq magZ11 angZ11\n100 0.99 -4\n200 0.80 -22\n300 0.707 -45\n400 0.40 -62\n500 0.01 -89",
            "!Example 10 (Version 2.0):\n!1-port Z-parameter file, multiple frequency points\n[Version] 2.0\n# MHz Z MA\n[Number of Ports] 1\n[Number of Frequencies] 5\n[Reference] 20.0\n[Network Data]\n!freq magZ11 angZ11\n100 74.25 -4\n200 60 -22\n300 53.025 -45\n400 30 -62\n500 0.75 -89\n[End]",
            "!Example 11 (Version 1.0):\n!2-port H-parameter file, single frequency point\n# kHz H MA R 1\n! freq magH11 angH11 magH21 angH21 magH12 angH12 magH22 angH22\n2 .95 -26 3.57 157 .04 76 .66 -14",
            "!Example 12 (Version 2.0):\n!2-port H-parameter file, single frequency point\n[Version] 2.0\n# kHz H MA R 1\n[Number of Ports] 2\n[Two-Port Data Order] 21_12\n[Number of Frequencies] 1\n[Matrix Format] Full\n[Network Data]\n! freq magH11 angH11 magH21 angH21 magH12 angH12 magH22 angH22\n2 .95 -26 3.57 157 .04 76 .66 -14\n[End]",
            "!Example 13 (Version 1.0):\n!2-port S-parameter file, three frequency points\n# GHz S RI R 50.0\n!freq ReS11 ImS11 ReS21 ImS21 ReS12 ImS12 ReS22 ImS22\n1.0000 0.3926 -0.1211 -0.0003 -0.0021 -0.0003 -0.0021 0.3926 -0.1211\n2.0000 0.3517 -0.3054 -0.0096 -0.0298 -0.0096 -0.0298 0.3517 -0.3054\n10.000 0.3419 0.3336 -0.0134 0.0379 -0.0134 0.0379 0.3419 0.3336",
            "!Example 14 (Version 1.0):\n! 4-port S-parameter data, taken at three frequency points\n! note that data points need not be aligned\n# GHz S MA R 50\n5.00000 0.60 161.24 0.40 -42.20 0.42 -66.58 0.53 -79.34 !row 1\n0.40 -42.20 0.60 161.20 0.53 -79.34 0.42 -66.58 !row 2\n0.42 -66.58 0.53 -79.34 0.60 161.24 0.40 -42.20 !row 3\n0.53 -79.34 0.42 -66.58 0.40 -42.20 0.60 161.24 !row 4\n6.00000 0.57 150.37 0.40 -44.34 0.41 -81.24 0.57 -95.77 !row 1\n0.40 -44.34 0.57 150.37 0.57 -95.77 0.41 -81.24 !row 2\n0.41 -81.24 0.57 -95.77 0.57 150.37 0.40 -44.34 !row 3\n0.57 -95.77 0.41 -81.24 0.40 -44.34 0.57 150.37 !row 4\n7.00000 0.50 136.69 0.45 -46.41 0.37 -99.09 0.62 -114.19 !row 1\n0.45 -46.41 0.50 136.69 0.62 -114.19 0.37 -99.09 !row 2\n0.37 -99.09 0.62 -114.19 0.50 136.69 0.45 -46.41 !row 3\n0.62 -114.19 0.37 -99.09 0.45 -46.41 0.50 136.69 !row 4",
            "!Example 16 (Version 2.0):\n!6-port component shown; note that all five ports are used in some\n!relationship\n[Version] 2.0\n# MHz Y RI R 50\n[Number of Ports] 6\n[Number of Frequencies] 1\n[Reference] 50 75 75 50 0.01 0.01\n[Mixed-Mode Order] D2,3 D6,5 C2,3 C6,5 S4 S1\n[Network Data]\n5.00 8.0 9.0 2.0 -1.0 3.0 -2.0 1.0 3.0 1.0 0.1 0.2 -0.2\n2.0 -1.0 7.0 7.0 1.8 -2.0 -1.0 -1.0 -0.5 0.5 0.2 -0.1\n3.0 -2.0 1.8 -2.0 5.8 6.0 1.2 0.8 0.9 0.7 0.3 -0.5\n1.0 3.0 -1.0 -1.0 1.2 0.8 6.3 8.0 2.0 -0.5 1.5 0.6\n1.0 0.1 -0.5 0.5 0.9 0.7 2.0 -0.5 4.7 -6.0 -1.0 2.0\n0.2 -0.2 0.2 -0.1 0.3 -0.5 1.5 0.6 -1.0 2.0 5.5 -7.0\n[End]",
            "!Example 17 (Version 2.0):\n!2-port network, S-parameter and noise data\n!Default MA format, GHz frequencies, 50 ohm reference, S-parameters\n[Version] 2.0\n#\n[Number of Ports] 2\n[Two-Port Data Order] 21_12\n[Number of Frequencies] 2\n[Number of Noise Frequencies] 2\n[Reference] 50 25.0\n[Network Data]\n2 .95 -26 3.57 157 .04 76 .66 -14\n22 .60 -144 1.30 40 .14 40 .56 -85\n[Noise Data]\n4 .7 .64 69 19\n18 2.7 .46 -33 20\n[End]",
            "!Example 18 (Version 1.0):\n!2-port network, S-parameter and noise data\n!Default MA format, GHz frequencies, 50 ohm reference, S-parameters\n#\n! NETWORK PARAMETERS\n2 .95 -26 3.57 157 .04 76 .66 -14\n22 .60 -144 1.30 40 .14 40 .56 -85\n! NOISE PARAMETERS\n4 .7 .64 69 .38\n18 2.7 .46 -33 .40",
            "!Example 19 (Version 2.0):\n!2-port network, S-parameter and noise data\n!Default MA format, GHz frequencies, 50 ohm reference, S-parameters\n[Version] 2.0\n#\n[Number of Ports] 2\n[Two-Port Data Order] 21_12\n[Number of Frequencies] 2\n[Number of Noise Frequencies] 2\n[Reference] 50 25.0\n[Network Data]\n! NETWORK PARAMETERS\n2 .95 -26 3.57 157 .04 76 .66 -14\n22 .60 -144 1.30 40 .14 40 .56 -85\n[Noise Data]\n! NOISE PARAMETERS\n4 .7 .64 69 19\n18 2.7 .46 -33 20\n[End]"
        ]

        parsedexamples = Any[
            JosephsonCircuits.TouchstoneFile([2.0e6], [0.874020294860635 - 0.18794819544685323im;;;], "mhz", "s", "ma", 50.0, 1.0, 1, "12_21", 1, 0, [50.0], String[], "Full", Tuple{Char, Vector{Int64}}[], ["Example 1:", "1-port S-parameter file, single frequency point", "freq magS11 angS11"], [2.0, 0.894, -12.136], Float64[]), 
            JosephsonCircuits.TouchstoneFile([2.0e6], [0.8743473504516138 - 0.18801852505875893im;;;], "mhz", "s", "db", 50.0, 1.0, 1, "12_21", 1, 0, [50.0], String[], "Full", Tuple{Char, Vector{Int64}}[], ["Example 1a:", "1-port S-parameter file, single frequency point", "freq magS11 angS11"], [2.0, -0.97, -12.136], Float64[]),
            JosephsonCircuits.TouchstoneFile([1.0e8, 2.0e8, 3.0e8, 4.0e8, 5.0e8], [74.06913073179194 - 5.179418175501303im;;; 55.63103127400724 - 22.47639560495472im;;; 37.494337072416684 - 37.49433707241668im;;; 14.084146883576725 - 26.488427785767808im;;; 0.013089304827962698 - 0.7498857713672935im], "mhz", "z", "ma", 75.0, 1.0, 1, "12_21", 5, 0, [75.0], String[], "Full", Tuple{Char, Vector{Int64}}[], ["Example 2:", "1-port Z-parameter file, multiple frequency points", "freq magZ11 angZ11"], [100.0, 0.99, -4.0, 200.0, 0.8, -22.0, 300.0, 0.707, -45.0, 400.0, 0.4, -62.0, 500.0, 0.01, -89.0], Float64[]),
            JosephsonCircuits.TouchstoneFile([2000.0], [0.8538543439842087 - 0.4164525894496235im 0.009676875823986707 + 0.03881182905103986im; -3.286202326825212 + 1.3949101287067074im 0.6598994788032183 - 0.011518588248607119im;;;], "khz", "h", "ma", 1.0, 1.0, 2, "21_12", 1, 0, [1.0, 1.0], String[], "Full", Tuple{Char, Vector{Int64}}[], ["Example 3:", "2-port H-parameter file, single frequency point", " freq magH11 angH11 magH21 angH21 magH12 angH12 magH22 angH22"], [2.0, 0.95, -26.0, 3.57, 157.0, 0.04, 76.0, 0.66, -1.0], Float64[]),
            JosephsonCircuits.TouchstoneFile([1.0e9, 2.0e9, 1.0e10], [0.3926 - 0.1211im -0.0003 - 0.0021im; -0.0003 - 0.0021im 0.3926 - 0.1211im;;; 0.3517 - 0.3054im -0.0096 - 0.0298im; -0.0096 - 0.0298im 0.3517 - 0.3054im;;; 0.3419 + 0.3336im -0.0134 + 0.0379im; -0.0134 + 0.0379im 0.3419 + 0.3336im], "ghz", "s", "ri", 50.0, 1.0, 2, "21_12", 3, 0, [50.0, 50.0], String[], "Full", Tuple{Char, Vector{Int64}}[], ["Example 4:", "2-port S-parameter file, three frequency points", "freq RelS11 ImS11 ReS21 ImS21 ReS12 ImS12 ReS22 ImS22"], [1.0, 0.3926, -0.1211, -0.0003, -0.0021, -0.0003, -0.0021, 0.3926, -0.1211, 2.0, 0.3517, -0.3054, -0.0096, -0.0298, -0.0096, -0.0298, 0.3517, -0.3054, 10.0, 0.3419, 0.3336, -0.0134, 0.0379, -0.0134, 0.0379, 0.3419, 0.3336], Float64[]),
            JosephsonCircuits.TouchstoneFile([5.0e9, 6.0e9, 7.0e9], [-0.5681244079815996 + 0.1929628385351877im 0.29632183851470006 - 0.2686882357291961im 0.16693665375723588 - 0.38539869438327984im 0.09803970583787712 - 0.5208533537179372im; 0.29632183851470006 - 0.2686882357291961im -0.5679895560694177 + 0.1933594171383067im 0.09803970583787712 - 0.5208533537179372im 0.16693665375723588 - 0.38539869438327984im; 0.16693665375723588 - 0.38539869438327984im 0.09803970583787712 - 0.5208533537179372im -0.5681244079815996 + 0.1929628385351877im 0.29632183851470006 - 0.2686882357291961im; 0.09803970583787712 - 0.5208533537179372im 0.16693665375723588 - 0.38539869438327984im 0.29632183851470006 - 0.2686882357291961im -0.5681244079815996 + 0.1929628385351877im;;; -0.495464624294036 + 0.28180632724119187im 0.286081989392916 - 0.2795659051905141im 0.062441313054034775 - 0.40521732739862937im -0.05730515806890161 - 0.5671120866801361im; 0.286081989392916 - 0.2795659051905141im -0.495464624294036 + 0.28180632724119187im -0.05730515806890161 - 0.5671120866801361im 0.062441313054034775 - 0.40521732739862937im; 0.062441313054034775 - 0.40521732739862937im -0.05730515806890161 - 0.5671120866801361im -0.495464624294036 + 0.28180632724119187im 0.286081989392916 - 0.2795659051905141im; -0.05730515806890161 - 0.5671120866801361im 0.062441313054034775 - 0.40521732739862937im 0.286081989392916 - 0.2795659051905141im -0.495464624294036 + 0.28180632724119187im;;; -0.3638265243449566 + 0.3429726813946975im 0.3102719136297667 - 0.32593149527549903im -0.05845471959176759 - 0.3653533163356367im -0.2540535762162701 - 0.565558821354352im; 0.3102719136297667 - 0.32593149527549903im -0.3638265243449566 + 0.3429726813946975im -0.2540535762162701 - 0.565558821354352im -0.05845471959176759 - 0.3653533163356367im; -0.05845471959176759 - 0.3653533163356367im -0.2540535762162701 - 0.565558821354352im -0.3638265243449566 + 0.3429726813946975im 0.3102719136297667 - 0.32593149527549903im; -0.2540535762162701 - 0.565558821354352im -0.05845471959176759 - 0.3653533163356367im 0.3102719136297667 - 0.32593149527549903im -0.3638265243449566 + 0.3429726813946975im], "ghz", "s", "ma", 50.0, 1.0, 4, "12_21", 3, 0, [50.0, 50.0, 50.0, 50.0], String[], "Full", Tuple{Char, Vector{Int64}}[], ["Example 5:", " 4-port S-parameter data, taken at three frequency points", "row 1", "row 2", "row 3", "row 4", "row 1", "row 2", "row 3", "row 4", "row 1", "row 2", "row 3", "row 4"], [5.0, 0.6, 161.24, 0.4, -42.2, 0.42, -66.58, 0.53, -79.34, 0.4, -42.2, 0.6, 161.2, 0.53, -79.34, 0.42, -66.58, 0.42, -66.58, 0.53, -79.34, 0.6, 161.24, 0.4, -42.2, 0.53, -79.34, 0.42, -66.58, 0.4, -42.2, 0.6, 161.24, 6.0, 0.57, 150.37, 0.4, -44.34, 0.41, -81.24, 0.57, -95.77, 0.4, -44.34, 0.57, 150.37, 0.57, -95.77, 0.41, -81.24, 0.41, -81.24, 0.57, -95.77, 0.57, 150.37, 0.4, -44.34, 0.57, -95.77, 0.41, -81.24, 0.4, -44.34, 0.57, 150.37, 7.0, 0.5, 136.69, 0.45, -46.41, 0.37, -99.09, 0.62, -114.19, 0.45, -46.41, 0.5, 136.69, 0.62, -114.19, 0.37, -99.09, 0.37, -99.09, 0.62, -114.19, 0.5, 136.69, 0.45, -46.41, 0.62, -114.19, 0.37, -99.09, 0.45, -46.41, 0.5, 136.69], Float64[]),
            JosephsonCircuits.TouchstoneFile([2.0e9, 2.2e10], [0.8538543439842087 - 0.4164525894496235im 0.009676875823986707 + 0.03881182905103986im; -3.286202326825212 + 1.3949101287067074im 0.6403951793421577 - 0.1596684510957807im;;; -0.48541019662496837 - 0.35267115137548394im 0.10724622203665693 + 0.0899902653561155im; 0.9958577760546714 + 0.835623892592501im 0.048807215938688565 - 0.5578690309313775im], "ghz", "s", "ma", 50.0, 1.0, 2, "21_12", 2, 2, [50.0, 50.0], String[], "Full", Tuple{Char, Vector{Int64}}[], ["Example 8:", "2-port network, S-parameter and noise data", " NOISE PARAMETERS"], [2.0, 0.95, -26.0, 3.57, 157.0, 0.04, 76.0, 0.66, -14.0, 22.0, 0.6, -144.0, 1.3, 40.0, 0.14, 40.0, 0.56, -85.0], [4.0, 0.7, 0.64, 69.0, 0.38, 18.0, 2.7, 0.46, -33.0, 0.4]),
            JosephsonCircuits.TouchstoneFile([5.0e9], [-0.5681244079815996 + 0.1929628385351877im 0.29632183851470006 - 0.2686882357291961im 0.16693665375723588 - 0.38539869438327984im 0.09803970583787712 - 0.5208533537179372im; 0.29632183851470006 - 0.2686882357291961im -0.5679895560694177 + 0.1933594171383067im 0.09803970583787712 - 0.5208533537179372im 0.16693665375723588 - 0.38539869438327984im; 0.16693665375723588 - 0.38539869438327984im 0.09803970583787712 - 0.5208533537179372im -0.5681244079815996 + 0.1929628385351877im 0.29632183851470006 - 0.2686882357291961im; 0.09803970583787712 - 0.5208533537179372im 0.16693665375723588 - 0.38539869438327984im 0.29632183851470006 - 0.2686882357291961im -0.5681244079815996 + 0.1929628385351877im;;;], "ghz", "s", "ma", 50.0, 2.0, 4, "12_21", 1, 0, [50.0, 75.0, 0.01, 0.01], String[], "Full", Tuple{Char, Vector{Int64}}[], ["Example 4 (Version 2.0):", " 4-port S-parameter data", " Default impedance is overridden by [Reference]", " Data cannot be represented using 1.0 syntax", " Note that the [Reference] keyword arguments appear on a separate line", "row 1", "row 2", "row 3", "row 4"], [5.0, 0.6, 161.24, 0.4, -42.2, 0.42, -66.58, 0.53, -79.34, 0.4, -42.2, 0.6, 161.2, 0.53, -79.34, 0.42, -66.58, 0.42, -66.58, 0.53, -79.34, 0.6, 161.24, 0.4, -42.2, 0.53, -79.34, 0.42, -66.58, 0.4, -42.2, 0.6, 161.24], Float64[]),
            JosephsonCircuits.TouchstoneFile([5.0e9], [-0.5681244079815996 + 0.1929628385351877im 0.29632183851470006 - 0.2686882357291961im 0.16693665375723588 - 0.38539869438327984im 0.09803970583787712 - 0.5208533537179372im; 0.29632183851470006 - 0.2686882357291961im -0.5679895560694177 + 0.1933594171383067im 0.09803970583787712 - 0.5208533537179372im 0.16693665375723588 - 0.38539869438327984im; 0.16693665375723588 - 0.38539869438327984im 0.09803970583787712 - 0.5208533537179372im -0.5681244079815996 + 0.1929628385351877im 0.29632183851470006 - 0.2686882357291961im; 0.09803970583787712 - 0.5208533537179372im 0.16693665375723588 - 0.38539869438327984im 0.29632183851470006 - 0.2686882357291961im -0.5681244079815996 + 0.1929628385351877im;;;], "ghz", "s", "ma", 50.0, 2.0, 4, "12_21", 1, 0, [50.0, 75.0, 0.01, 0.01], String[], "Full", Tuple{Char, Vector{Int64}}[], ["Example 5 (Version 2.0):", " 4-port S-parameter data", " Default impedance is overridden by the [Reference] keyword arguments", " Data cannot be represented using 1.0 syntax", "row 1", "row 2", "row 3", "row 4"], [5.0, 0.6, 161.24, 0.4, -42.2, 0.42, -66.58, 0.53, -79.34, 0.4, -42.2, 0.6, 161.2, 0.53, -79.34, 0.42, -66.58, 0.42, -66.58, 0.53, -79.34, 0.6, 161.24, 0.4, -42.2, 0.53, -79.34, 0.42, -66.58, 0.4, -42.2, 0.6, 161.24], Float64[]),
            JosephsonCircuits.TouchstoneFile([5.0e9], [-0.5681244079815996 + 0.1929628385351877im 0.29632183851470006 - 0.2686882357291961im 0.16693665375723588 - 0.38539869438327984im 0.09803970583787712 - 0.5208533537179372im; 0.29632183851470006 - 0.2686882357291961im -0.5679895560694177 + 0.1933594171383067im 0.09803970583787712 - 0.5208533537179372im 0.16693665375723588 - 0.38539869438327984im; 0.16693665375723588 - 0.38539869438327984im 0.09803970583787712 - 0.5208533537179372im -0.5681244079815996 + 0.1929628385351877im 0.29632183851470006 - 0.2686882357291961im; 0.09803970583787712 - 0.5208533537179372im 0.16693665375723588 - 0.38539869438327984im 0.29632183851470006 - 0.2686882357291961im -0.5681244079815996 + 0.1929628385351877im;;;], "ghz", "s", "ma", 50.0, 2.0, 4, "12_21", 1, 0, [50.0, 75.0, 0.01, 0.01], String[], "Full", Tuple{Char, Vector{Int64}}[], ["Example 5 (Version 2.0):", " 4-port S-parameter data", " Default impedance is overridden by the [Reference] keyword arguments", " Data cannot be represented using 1.0 syntax", "row 1", "row 2", "row 3", "row 4"], [5.0, 0.6, 161.24, 0.4, -42.2, 0.42, -66.58, 0.53, -79.34, 0.4, -42.2, 0.6, 161.2, 0.53, -79.34, 0.42, -66.58, 0.42, -66.58, 0.53, -79.34, 0.6, 161.24, 0.4, -42.2, 0.53, -79.34, 0.42, -66.58, 0.4, -42.2, 0.6, 161.24], Float64[]),
            JosephsonCircuits.TouchstoneFile([5.0e9], [-0.5681244079815996 + 0.1929628385351877im 0.29632183851470006 - 0.2686882357291961im 0.16693665375723588 - 0.38539869438327984im 0.09803970583787712 - 0.5208533537179372im; 0.29632183851470006 - 0.2686882357291961im -0.5679895560694177 + 0.1933594171383067im 0.09803970583787712 - 0.5208533537179372im 0.16693665375723588 - 0.38539869438327984im; 0.16693665375723588 - 0.38539869438327984im 0.09803970583787712 - 0.5208533537179372im -0.5681244079815996 + 0.1929628385351877im 0.29632183851470006 - 0.2686882357291961im; 0.09803970583787712 - 0.5208533537179372im 0.16693665375723588 - 0.38539869438327984im 0.29632183851470006 - 0.2686882357291961im -0.5681244079815996 + 0.1929628385351877im;;;], "ghz", "s", "ma", 50.0, 2.0, 4, "12_21", 1, 0, [50.0, 75.0, 0.01, 0.01], String[], "Lower", Tuple{Char, Vector{Int64}}[], ["Example 6 (Version 2.0):", " 4-port S-parameter data", " Default impedance is overridden by the [Reference] keyword arguments", " Note that [Reference] arguments are split across two lines", " Data cannot be represented using 1.0 syntax", "row 1", "row 2", "row 3", "row 4"], [5.0, 0.6, 161.24, 0.4, -42.2, 0.6, 161.2, 0.42, -66.58, 0.53, -79.34, 0.6, 161.24, 0.53, -79.34, 0.42, -66.58, 0.4, -42.2, 0.6, 161.24], Float64[]),
            JosephsonCircuits.TouchstoneFile([1.0e8, 2.0e8, 3.0e8, 4.0e8, 5.0e8], [74.06913073179194 - 5.179418175501303im;;; 55.63103127400724 - 22.47639560495472im;;; 37.494337072416684 - 37.49433707241668im;;; 14.084146883576725 - 26.488427785767808im;;; 0.013089304827962698 - 0.7498857713672935im], "mhz", "z", "ma", 50.0, 2.0, 1, "12_21", 5, 0, [20.0], String[], "Full", Tuple{Char, Vector{Int64}}[], ["Example 7 (Version 2.0):", "1-port Z-parameter file, multiple frequency points", "freq magZ11 angZ11"], [100.0, 74.25, -4.0, 200.0, 60.0, -22.0, 300.0, 53.025, -45.0, 400.0, 30.0, -62.0, 500.0, 0.75, -89.0], Float64[]),
            JosephsonCircuits.TouchstoneFile([2.0e6], [0.874020294860635 - 0.18794819544685323im;;;], "mhz", "s", "ma", 50.0, 1.0, 1, "12_21", 1, 0, [50.0], String[], "Full", Tuple{Char, Vector{Int64}}[], ["Example 8 (Version 1.0):", "1-port S-parameter file, single frequency point", "freq magS11 angS11"], [2.0, 0.894, -12.136], Float64[]),
            JosephsonCircuits.TouchstoneFile([1.0e8, 2.0e8, 3.0e8, 4.0e8, 5.0e8], [74.06913073179194 - 5.179418175501303im;;; 55.63103127400724 - 22.47639560495472im;;; 37.494337072416684 - 37.49433707241668im;;; 14.084146883576725 - 26.488427785767808im;;; 0.013089304827962698 - 0.7498857713672935im], "mhz", "z", "ma", 75.0, 1.0, 1, "12_21", 5, 0, [75.0], String[], "Full", Tuple{Char, Vector{Int64}}[], ["Example 9 (Version 1.0):", "1-port Z-parameter file, multiple frequency points", "freq magZ11 angZ11"], [100.0, 0.99, -4.0, 200.0, 0.8, -22.0, 300.0, 0.707, -45.0, 400.0, 0.4, -62.0, 500.0, 0.01, -89.0], Float64[]),
            JosephsonCircuits.TouchstoneFile([1.0e8, 2.0e8, 3.0e8, 4.0e8, 5.0e8], [74.06913073179194 - 5.179418175501303im;;; 55.63103127400724 - 22.47639560495472im;;; 37.494337072416684 - 37.49433707241668im;;; 14.084146883576725 - 26.488427785767808im;;; 0.013089304827962698 - 0.7498857713672935im], "mhz", "z", "ma", 50.0, 2.0, 1, "12_21", 5, 0, [20.0], String[], "Full", Tuple{Char, Vector{Int64}}[], ["Example 10 (Version 2.0):", "1-port Z-parameter file, multiple frequency points", "freq magZ11 angZ11"], [100.0, 74.25, -4.0, 200.0, 60.0, -22.0, 300.0, 53.025, -45.0, 400.0, 30.0, -62.0, 500.0, 0.75, -89.0], Float64[]),
            JosephsonCircuits.TouchstoneFile([2000.0], [0.8538543439842087 - 0.4164525894496235im 0.009676875823986707 + 0.03881182905103986im; -3.286202326825212 + 1.3949101287067074im 0.6403951793421577 - 0.1596684510957807im;;;], "khz", "h", "ma", 1.0, 1.0, 2, "21_12", 1, 0, [1.0, 1.0], String[], "Full", Tuple{Char, Vector{Int64}}[], ["Example 11 (Version 1.0):", "2-port H-parameter file, single frequency point", " freq magH11 angH11 magH21 angH21 magH12 angH12 magH22 angH22"], [2.0, 0.95, -26.0, 3.57, 157.0, 0.04, 76.0, 0.66, -14.0], Float64[]),
            JosephsonCircuits.TouchstoneFile([2000.0], [0.8538543439842087 - 0.4164525894496235im 0.009676875823986707 + 0.03881182905103986im; -3.286202326825212 + 1.3949101287067074im 0.6403951793421577 - 0.1596684510957807im;;;], "khz", "h", "ma", 1.0, 2.0, 2, "21_12", 1, 0, [1.0, 1.0], String[], "Full", Tuple{Char, Vector{Int64}}[], ["Example 12 (Version 2.0):", "2-port H-parameter file, single frequency point", " freq magH11 angH11 magH21 angH21 magH12 angH12 magH22 angH22"], [2.0, 0.95, -26.0, 3.57, 157.0, 0.04, 76.0, 0.66, -14.0], Float64[]),
            JosephsonCircuits.TouchstoneFile([1.0e9, 2.0e9, 1.0e10], [0.3926 - 0.1211im -0.0003 - 0.0021im; -0.0003 - 0.0021im 0.3926 - 0.1211im;;; 0.3517 - 0.3054im -0.0096 - 0.0298im; -0.0096 - 0.0298im 0.3517 - 0.3054im;;; 0.3419 + 0.3336im -0.0134 + 0.0379im; -0.0134 + 0.0379im 0.3419 + 0.3336im], "ghz", "s", "ri", 50.0, 1.0, 2, "21_12", 3, 0, [50.0, 50.0], String[], "Full", Tuple{Char, Vector{Int64}}[], ["Example 13 (Version 1.0):", "2-port S-parameter file, three frequency points", "freq ReS11 ImS11 ReS21 ImS21 ReS12 ImS12 ReS22 ImS22"], [1.0, 0.3926, -0.1211, -0.0003, -0.0021, -0.0003, -0.0021, 0.3926, -0.1211, 2.0, 0.3517, -0.3054, -0.0096, -0.0298, -0.0096, -0.0298, 0.3517, -0.3054, 10.0, 0.3419, 0.3336, -0.0134, 0.0379, -0.0134, 0.0379, 0.3419, 0.3336], Float64[]),
            JosephsonCircuits.TouchstoneFile([5.0e9, 6.0e9, 7.0e9], [-0.5681244079815996 + 0.1929628385351877im 0.29632183851470006 - 0.2686882357291961im 0.16693665375723588 - 0.38539869438327984im 0.09803970583787712 - 0.5208533537179372im; 0.29632183851470006 - 0.2686882357291961im -0.5679895560694177 + 0.1933594171383067im 0.09803970583787712 - 0.5208533537179372im 0.16693665375723588 - 0.38539869438327984im; 0.16693665375723588 - 0.38539869438327984im 0.09803970583787712 - 0.5208533537179372im -0.5681244079815996 + 0.1929628385351877im 0.29632183851470006 - 0.2686882357291961im; 0.09803970583787712 - 0.5208533537179372im 0.16693665375723588 - 0.38539869438327984im 0.29632183851470006 - 0.2686882357291961im -0.5681244079815996 + 0.1929628385351877im;;; -0.495464624294036 + 0.28180632724119187im 0.286081989392916 - 0.2795659051905141im 0.062441313054034775 - 0.40521732739862937im -0.05730515806890161 - 0.5671120866801361im; 0.286081989392916 - 0.2795659051905141im -0.495464624294036 + 0.28180632724119187im -0.05730515806890161 - 0.5671120866801361im 0.062441313054034775 - 0.40521732739862937im; 0.062441313054034775 - 0.40521732739862937im -0.05730515806890161 - 0.5671120866801361im -0.495464624294036 + 0.28180632724119187im 0.286081989392916 - 0.2795659051905141im; -0.05730515806890161 - 0.5671120866801361im 0.062441313054034775 - 0.40521732739862937im 0.286081989392916 - 0.2795659051905141im -0.495464624294036 + 0.28180632724119187im;;; -0.3638265243449566 + 0.3429726813946975im 0.3102719136297667 - 0.32593149527549903im -0.05845471959176759 - 0.3653533163356367im -0.2540535762162701 - 0.565558821354352im; 0.3102719136297667 - 0.32593149527549903im -0.3638265243449566 + 0.3429726813946975im -0.2540535762162701 - 0.565558821354352im -0.05845471959176759 - 0.3653533163356367im; -0.05845471959176759 - 0.3653533163356367im -0.2540535762162701 - 0.565558821354352im -0.3638265243449566 + 0.3429726813946975im 0.3102719136297667 - 0.32593149527549903im; -0.2540535762162701 - 0.565558821354352im -0.05845471959176759 - 0.3653533163356367im 0.3102719136297667 - 0.32593149527549903im -0.3638265243449566 + 0.3429726813946975im], "ghz", "s", "ma", 50.0, 1.0, 4, "12_21", 3, 0, [50.0, 50.0, 50.0, 50.0], String[], "Full", Tuple{Char, Vector{Int64}}[], ["Example 14 (Version 1.0):", " 4-port S-parameter data, taken at three frequency points", " note that data points need not be aligned", "row 1", "row 2", "row 3", "row 4", "row 1", "row 2", "row 3", "row 4", "row 1", "row 2", "row 3", "row 4"], [5.0, 0.6, 161.24, 0.4, -42.2, 0.42, -66.58, 0.53, -79.34, 0.4, -42.2, 0.6, 161.2, 0.53, -79.34, 0.42, -66.58, 0.42, -66.58, 0.53, -79.34, 0.6, 161.24, 0.4, -42.2, 0.53, -79.34, 0.42, -66.58, 0.4, -42.2, 0.6, 161.24, 6.0, 0.57, 150.37, 0.4, -44.34, 0.41, -81.24, 0.57, -95.77, 0.4, -44.34, 0.57, 150.37, 0.57, -95.77, 0.41, -81.24, 0.41, -81.24, 0.57, -95.77, 0.57, 150.37, 0.4, -44.34, 0.57, -95.77, 0.41, -81.24, 0.4, -44.34, 0.57, 150.37, 7.0, 0.5, 136.69, 0.45, -46.41, 0.37, -99.09, 0.62, -114.19, 0.45, -46.41, 0.5, 136.69, 0.62, -114.19, 0.37, -99.09, 0.37, -99.09, 0.62, -114.19, 0.5, 136.69, 0.45, -46.41, 0.62, -114.19, 0.37, -99.09, 0.45, -46.41, 0.5, 136.69], Float64[]), 
            JosephsonCircuits.TouchstoneFile([5.0e6], [8.0 + 9.0im 2.0 - 1.0im 3.0 - 2.0im 1.0 + 3.0im 1.0 + 0.1im 0.2 - 0.2im; 2.0 - 1.0im 7.0 + 7.0im 1.8 - 2.0im -1.0 - 1.0im -0.5 + 0.5im 0.2 - 0.1im; 3.0 - 2.0im 1.8 - 2.0im 5.8 + 6.0im 1.2 + 0.8im 0.9 + 0.7im 0.3 - 0.5im; 1.0 + 3.0im -1.0 - 1.0im 1.2 + 0.8im 6.3 + 8.0im 2.0 - 0.5im 1.5 + 0.6im; 1.0 + 0.1im -0.5 + 0.5im 0.9 + 0.7im 2.0 - 0.5im 4.7 - 6.0im -1.0 + 2.0im; 0.2 - 0.2im 0.2 - 0.1im 0.3 - 0.5im 1.5 + 0.6im -1.0 + 2.0im 5.5 - 7.0im;;;], "mhz", "y", "ri", 50.0, 2.0, 6, "12_21", 1, 0, [50.0, 75.0, 75.0, 50.0, 0.01, 0.01], String[], "Full", [('D', [2, 3]), ('D', [6, 5]), ('C', [2, 3]), ('C', [6, 5]), ('S', [4]), ('S', [1])], ["Example 16 (Version 2.0):", "6-port component shown; note that all five ports are used in some", "relationship"], [5.0, 8.0, 9.0, 2.0, -1.0, 3.0, -2.0, 1.0, 3.0, 1.0, 0.1, 0.2, -0.2, 2.0, -1.0, 7.0, 7.0, 1.8, -2.0, -1.0, -1.0, -0.5, 0.5, 0.2, -0.1, 3.0, -2.0, 1.8, -2.0, 5.8, 6.0, 1.2, 0.8, 0.9, 0.7, 0.3, -0.5, 1.0, 3.0, -1.0, -1.0, 1.2, 0.8, 6.3, 8.0, 2.0, -0.5, 1.5, 0.6, 1.0, 0.1, -0.5, 0.5, 0.9, 0.7, 2.0, -0.5, 4.7, -6.0, -1.0, 2.0, 0.2, -0.2, 0.2, -0.1, 0.3, -0.5, 1.5, 0.6, -1.0, 2.0, 5.5, -7.0], Float64[]),
            JosephsonCircuits.TouchstoneFile([2.0e9, 2.2e10], [0.8538543439842087 - 0.4164525894496235im 0.009676875823986707 + 0.03881182905103986im; -3.286202326825212 + 1.3949101287067074im 0.6403951793421577 - 0.1596684510957807im;;; -0.48541019662496837 - 0.35267115137548394im 0.10724622203665693 + 0.0899902653561155im; 0.9958577760546714 + 0.835623892592501im 0.048807215938688565 - 0.5578690309313775im], "ghz", "s", "ma", 50.0, 2.0, 2, "21_12", 2, 2, [50.0, 25.0], String[], "Full", Tuple{Char, Vector{Int64}}[], ["Example 17 (Version 2.0):", "2-port network, S-parameter and noise data", "Default MA format, GHz frequencies, 50 ohm reference, S-parameters"], [2.0, 0.95, -26.0, 3.57, 157.0, 0.04, 76.0, 0.66, -14.0, 22.0, 0.6, -144.0, 1.3, 40.0, 0.14, 40.0, 0.56, -85.0], [4.0, 0.7, 0.64, 69.0, 19.0, 18.0, 2.7, 0.46, -33.0, 20.0]),
            JosephsonCircuits.TouchstoneFile([2.0e9, 2.2e10], [0.8538543439842087 - 0.4164525894496235im 0.009676875823986707 + 0.03881182905103986im; -3.286202326825212 + 1.3949101287067074im 0.6403951793421577 - 0.1596684510957807im;;; -0.48541019662496837 - 0.35267115137548394im 0.10724622203665693 + 0.0899902653561155im; 0.9958577760546714 + 0.835623892592501im 0.048807215938688565 - 0.5578690309313775im], "ghz", "s", "ma", 50.0, 1.0, 2, "21_12", 2, 2, [50.0, 50.0], String[], "Full", Tuple{Char, Vector{Int64}}[], ["Example 18 (Version 1.0):", "2-port network, S-parameter and noise data", "Default MA format, GHz frequencies, 50 ohm reference, S-parameters", " NETWORK PARAMETERS", " NOISE PARAMETERS"], [2.0, 0.95, -26.0, 3.57, 157.0, 0.04, 76.0, 0.66, -14.0, 22.0, 0.6, -144.0, 1.3, 40.0, 0.14, 40.0, 0.56, -85.0], [4.0, 0.7, 0.64, 69.0, 0.38, 18.0, 2.7, 0.46, -33.0, 0.4]),
            JosephsonCircuits.TouchstoneFile([2.0e9, 2.2e10], [0.8538543439842087 - 0.4164525894496235im 0.009676875823986707 + 0.03881182905103986im; -3.286202326825212 + 1.3949101287067074im 0.6403951793421577 - 0.1596684510957807im;;; -0.48541019662496837 - 0.35267115137548394im 0.10724622203665693 + 0.0899902653561155im; 0.9958577760546714 + 0.835623892592501im 0.048807215938688565 - 0.5578690309313775im], "ghz", "s", "ma", 50.0, 2.0, 2, "21_12", 2, 2, [50.0, 25.0], String[], "Full", Tuple{Char, Vector{Int64}}[], ["Example 19 (Version 2.0):", "2-port network, S-parameter and noise data", "Default MA format, GHz frequencies, 50 ohm reference, S-parameters", " NETWORK PARAMETERS", " NOISE PARAMETERS"], [2.0, 0.95, -26.0, 3.57, 157.0, 0.04, 76.0, 0.66, -14.0, 22.0, 0.6, -144.0, 1.3, 40.0, 0.14, 40.0, 0.56, -85.0], [4.0, 0.7, 0.64, 69.0, 19.0, 18.0, 2.7, 0.46, -33.0, 20.0]),
        ]

        for i in 1:length(examples)
            parsedexample=JosephsonCircuits.touchstone_parse(IOBuffer(examples[i]));
            for field in fieldnames(JosephsonCircuits.TouchstoneFile)
              @test getfield(parsedexample,field) == getfield(parsedexamples[i],field)
            end
        end

    end


    @testset "touchstone_parse errors" begin

        # empty input
        message = "Input is empty and thus not a valid touchstone file."
        @test_throws(
            ErrorException(message),
            JosephsonCircuits.touchstone_parse(IOBuffer())
        )

        # two version lines
        example = "!Example 4 (Version 2.0):\n! 4-port S-parameter data\n! Default impedance is overridden by [Reference]\n! Data cannot be represented using 1.0 syntax\n! Note that the [Reference] keyword arguments appear on a separate line\n[Version] 2.0\n[Version] 2.0\n# GHz S MA R 50\n[Number of Ports] 4\n[Reference]\n50 75 0.01 0.01\n[Number of Frequencies] 1\n[Network Data]\n5.00000 0.60 161.24 0.40 -42.20 0.42 -66.58 0.53 -79.34 !row 1\n0.40 -42.20 0.60 161.20 0.53 -79.34 0.42 -66.58 !row 2\n0.42 -66.58 0.53 -79.34 0.60 161.24 0.40 -42.20 !row 3\n0.53 -79.34 0.42 -66.58 0.40 -42.20 0.60 161.24 !row 4"
        message = "Only one [Version] keyword allowed: [version] 2.0"
        @test_throws(
            ErrorException(message),
            JosephsonCircuits.touchstone_parse(IOBuffer(example))
        )

        # second option line
        # the spec says a second options line should be ignored, but the 
        # golden parser tsck2 throws an error.
        example = "!Example 4 (Version 2.0):\n! 4-port S-parameter data\n! Default impedance is overridden by [Reference]\n! Data cannot be represented using 1.0 syntax\n! Note that the [Reference] keyword arguments appear on a separate line\n[Version] 2.0\n# GHz S MA R 50\n# GHz S MA R 50\n[Number of Ports] 4\n[Reference]\n50 75 0.01 0.01\n[Number of Frequencies] 1\n[Network Data]\n5.00000 0.60 161.24 0.40 -42.20 0.42 -66.58 0.53 -79.34 !row 1\n0.40 -42.20 0.60 161.20 0.53 -79.34 0.42 -66.58 !row 2\n0.42 -66.58 0.53 -79.34 0.60 161.24 0.40 -42.20 !row 3\n0.53 -79.34 0.42 -66.58 0.40 -42.20 0.60 161.24 !row 4"
        message = "Invalid secondary options line: # ghz s ma r 50"
        @test_throws(
            ErrorException(message),
            JosephsonCircuits.touchstone_parse(IOBuffer(example))
        )

        # two [Number of Ports] keywords
        example = "!Example 4 (Version 2.0):\n! 4-port S-parameter data\n! Default impedance is overridden by [Reference]\n! Data cannot be represented using 1.0 syntax\n! Note that the [Reference] keyword arguments appear on a separate line\n[Version] 2.0\n# GHz S MA R 50\n[Number of Ports] 4\n[Number of Ports] 4\n[Reference]\n50 75 0.01 0.01\n[Number of Frequencies] 1\n[Network Data]\n5.00000 0.60 161.24 0.40 -42.20 0.42 -66.58 0.53 -79.34 !row 1\n0.40 -42.20 0.60 161.20 0.53 -79.34 0.42 -66.58 !row 2\n0.42 -66.58 0.53 -79.34 0.60 161.24 0.40 -42.20 !row 3\n0.53 -79.34 0.42 -66.58 0.40 -42.20 0.60 161.24 !row 4"
        message = "Only one [Number of Ports] keyword allowed: [number of ports] 4"
        @test_throws(
            ErrorException(message),
            JosephsonCircuits.touchstone_parse(IOBuffer(example))
        )

        # two [Two-Port Data Order] keywords
        example = "!Example 17 (Version 2.0):\n!2-port network, S-parameter and noise data\n!Default MA format, GHz frequencies, 50 ohm reference, S-parameters\n[Version] 2.0\n# ghz s ma R 50.0\n[Number of Ports] 2\n[Two-Port Data Order] 21_12\n[Two-Port Data Order] 21_12\n[Number of Frequencies] 2\n[Number of Noise Frequencies] 2\n[Reference] 50.0 25.0\n[Network Data]\n! freq mags11 angs11 mags21 angs21 mags12 angs12 mags22 angs22 \n2.0 0.95 -26.0 3.57 157.0 0.04 76.0 0.66 -14.0\n22.0 0.6 -144.0 1.3 40.0 0.14 40.0 0.56 -85.0\n[Noise Data]\n4.0 0.7 0.64 69.0 19.0\n18.0 2.7 0.46 -33.0 20.0\n[End]"
        message = "Only one [Two-Port Data Order] line allowed: [two-port data order] 21_12"
        @test_throws(
            ErrorException(message),
            JosephsonCircuits.touchstone_parse(IOBuffer(example))
        )

        # [Two-Port Data Order] line for 4-port file
        example = "!Example 4 (Version 2.0):\n! 4-port S-parameter data\n! Default impedance is overridden by [Reference]\n! Data cannot be represented using 1.0 syntax\n! Note that the [Reference] keyword arguments appear on a separate line\n[Version] 2.0\n# GHz S MA R 50\n[Number of Ports] 4\n[Two-Port Data Order] 21_12\n[Reference]\n50 75 0.01 0.01\n[Number of Frequencies] 1\n[Network Data]\n5.00000 0.60 161.24 0.40 -42.20 0.42 -66.58 0.53 -79.34 !row 1\n0.40 -42.20 0.60 161.20 0.53 -79.34 0.42 -66.58 !row 2\n0.42 -66.58 0.53 -79.34 0.60 161.24 0.40 -42.20 !row 3\n0.53 -79.34 0.42 -66.58 0.40 -42.20 0.60 161.24 !row 4"
        message = "[Two-Port Data Order] is only allowed if [Number of Ports] is 2: [two-port data order] 21_12" 
        @test_throws(
            ErrorException(message),
            JosephsonCircuits.touchstone_parse(IOBuffer(example))
        )

        # [Begin Information] before [Number of Ports]
        example = "!Example 17 (Version 2.0):\n!2-port network, S-parameter and noise data\n!Default MA format, GHz frequencies, 50 ohm reference, S-parameters\n[Version] 2.0\n# ghz s ma R 50.0\n[Begin Information]\n[Number of Ports] 2\n[Two-Port Data Order] 21_12\n[Number of Frequencies] 2\n[Number of Noise Frequencies] 2\n[Reference] 50.0 25.0\n[Network Data]\n! freq mags11 angs11 mags21 angs21 mags12 angs12 mags22 angs22 \n2.0 0.95 -26.0 3.57 157.0 0.04 76.0 0.66 -14.0\n22.0 0.6 -144.0 1.3 40.0 0.14 40.0 0.56 -85.0\n[Noise Data]\n4.0 0.7 0.64 69.0 19.0\n18.0 2.7 0.46 -33.0 20.0\n[End]"
        message = "[Number of Ports] must be before [Begin Information]"
        @test_throws(
            ErrorException(message),
            JosephsonCircuits.touchstone_parse(IOBuffer(example))
        )

        # two [End Information] keywords
        example = "!Example 17 (Version 2.0):\n!2-port network, S-parameter and noise data\n!Default MA format, GHz frequencies, 50 ohm reference, S-parameters\n[Version] 2.0\n# ghz s ma R 50.0\n[Number of Ports] 2\n[Two-Port Data Order] 21_12\n[Begin Information]\ninformation\n[End Information]\n[End Information]\n[Number of Frequencies] 2\n[Number of Noise Frequencies] 2\n[Reference] 50.0 25.0\n[Network Data]\n! freq mags11 angs11 mags21 angs21 mags12 angs12 mags22 angs22 \n2.0 0.95 -26.0 3.57 157.0 0.04 76.0 0.66 -14.0\n22.0 0.6 -144.0 1.3 40.0 0.14 40.0 0.56 -85.0\n[Noise Data]\n4.0 0.7 0.64 69.0 19.0\n18.0 2.7 0.46 -33.0 20.0\n[End]"
        message = "The [End Information] line should have been parsed by the parseinformation() function."
        @test_throws(
            ErrorException(message),
            JosephsonCircuits.touchstone_parse(IOBuffer(example))
        )

        # two [Reference] keywords
        example = "!Example 17 (Version 2.0):\n!2-port network, S-parameter and noise data\n!Default MA format, GHz frequencies, 50 ohm reference, S-parameters\n[Version] 2.0\n# ghz s ma R 50.0\n[Number of Ports] 2\n[Two-Port Data Order] 21_12\n[Number of Frequencies] 2\n[Number of Noise Frequencies] 2\n[Reference] 50.0 25.0\n[Reference] 50.0 25.0\n[Network Data]\n! freq mags11 angs11 mags21 angs21 mags12 angs12 mags22 angs22 \n2.0 0.95 -26.0 3.57 157.0 0.04 76.0 0.66 -14.0\n22.0 0.6 -144.0 1.3 40.0 0.14 40.0 0.56 -85.0\n[Noise Data]\n4.0 0.7 0.64 69.0 19.0\n18.0 2.7 0.46 -33.0 20.0\n[End]"
        message = "Only one [Reference] keyword allowed: [reference] 50.0 25.0"
        @test_throws(
            ErrorException(message),
            JosephsonCircuits.touchstone_parse(IOBuffer(example))
        )

        # [Number of Ports] must be before [Reference]
        example = "!Example 17 (Version 2.0):\n!2-port network, S-parameter and noise data\n!Default MA format, GHz frequencies, 50 ohm reference, S-parameters\n[Version] 2.0\n# ghz s ma R 50.0\n[Reference] 50.0 25.0\n[Number of Ports] 2\n[Two-Port Data Order] 21_12\n[Number of Frequencies] 2\n[Number of Noise Frequencies] 2\n[Network Data]\n! freq mags11 angs11 mags21 angs21 mags12 angs12 mags22 angs22 \n2.0 0.95 -26.0 3.57 157.0 0.04 76.0 0.66 -14.0\n22.0 0.6 -144.0 1.3 40.0 0.14 40.0 0.56 -85.0\n[Noise Data]\n4.0 0.7 0.64 69.0 19.0\n18.0 2.7 0.46 -33.0 20.0\n[End]"
        message = "[Number of Ports] must be before [Reference]"
        @test_throws(
            ErrorException(message),
            JosephsonCircuits.touchstone_parse(IOBuffer(example))
        )

        # Unknown keyword
        example = "!Example 17 (Version 2.0):\n!2-port network, S-parameter and noise data\n!Default MA format, GHz frequencies, 50 ohm reference, S-parameters\n[Version] 2.0\n# ghz s ma R 50.0\n[Number of Ports] 2\n[Two-Port Data Order] 21_12\n[Number of Frequencies] 2\n[Number of Noise Frequencies] 2\n[Reference] 50.0 25.0\n[Unknown keyword example] unknown line example\n[Network Data]\n! freq mags11 angs11 mags21 angs21 mags12 angs12 mags22 angs22 \n2.0 0.95 -26.0 3.57 157.0 0.04 76.0 0.66 -14.0\n22.0 0.6 -144.0 1.3 40.0 0.14 40.0 0.56 -85.0\n[Noise Data]\n4.0 0.7 0.64 69.0 19.0\n18.0 2.7 0.46 -33.0 20.0\n[End]"
        message = "Unknown line type: [unknown keyword example] unknown line example"
        @test_throws(
            ErrorException(message),
            JosephsonCircuits.touchstone_parse(IOBuffer(example))
        )

        # Number of ports in data and header not consistent
        example = "!Example 4 (Version 2.0):\n! 4-port S-parameter data\n! Default impedance is overridden by [Reference]\n! Data cannot be represented using 1.0 syntax\n! Note that the [Reference] keyword arguments appear on a separate line\n[Version] 2.0\n# GHz S MA R 50\n[Number of Ports] 2\n[Reference]\n50 75\n[Number of Frequencies] 1\n[Network Data]\n5.00000 0.60 161.24 0.40 -42.20 0.42 -66.58 0.53 -79.34 !row 1\n0.40 -42.20 0.60 161.20 0.53 -79.34 0.42 -66.58 !row 2\n0.42 -66.58 0.53 -79.34 0.60 161.24 0.40 -42.20 !row 3\n0.53 -79.34 0.42 -66.58 0.40 -42.20 0.60 161.24 !row 4"
        message = "Number of ports not consistent between header and network data."
        @test_throws(
            ErrorException(message),
            JosephsonCircuits.touchstone_parse(IOBuffer(example))
        )

        # More impedance values on refernece line than number of ports
        example = "!Example 4 (Version 2.0):\n! 4-port S-parameter data\n! Default impedance is overridden by [Reference]\n! Data cannot be represented using 1.0 syntax\n! Note that the [Reference] keyword arguments appear on a separate line\n[Version] 2.0\n# GHz S MA R 50\n[Number of Ports] 2\n[Reference]\n50 75 0.01 0.01\n[Number of Frequencies] 1\n[Network Data]\n5.00000 0.60 161.24 0.40 -42.20 0.42 -66.58 0.53 -79.34 !row 1\n0.40 -42.20 0.60 161.20 0.53 -79.34 0.42 -66.58 !row 2\n0.42 -66.58 0.53 -79.34 0.60 161.24 0.40 -42.20 !row 3\n0.53 -79.34 0.42 -66.58 0.40 -42.20 0.60 161.24 !row 4"
        message = "Too many values on [Reference] line: 50 75 0.01 0.01"
        @test_throws(
            ErrorException(message),
            JosephsonCircuits.touchstone_parse(IOBuffer(example))
        )

        # Too few impedance values on reference line relative to number of ports. 
        # Can I make a better error for this?
        example = "!Example 4 (Version 2.0):\n! 4-port S-parameter data\n! Default impedance is overridden by [Reference]\n! Data cannot be represented using 1.0 syntax\n! Note that the [Reference] keyword arguments appear on a separate line\n[Version] 2.0\n# GHz S MA R 50\n[Number of Ports] 4\n[Reference]\n50 75\n[Number of Frequencies] 1\n[Network Data]\n5.00000 0.60 161.24 0.40 -42.20 0.42 -66.58 0.53 -79.34 !row 1\n0.40 -42.20 0.60 161.20 0.53 -79.34 0.42 -66.58 !row 2\n0.42 -66.58 0.53 -79.34 0.60 161.24 0.40 -42.20 !row 3\n0.53 -79.34 0.42 -66.58 0.40 -42.20 0.60 161.24 !row 4"
        message = """cannot parse "[number" as Float64"""
        @test_throws(
            ArgumentError(message),
            JosephsonCircuits.touchstone_parse(IOBuffer(example))
        )

        # # Wrong version number
        # example = "!Example 4 (Version 2.0):\n! 4-port S-parameter data\n! Default impedance is overridden by [Reference]\n! Data cannot be represented using 1.0 syntax\n! Note that the [Reference] keyword arguments appear on a separate line\n[Version] 1.0\n# GHz S MA R 50\n[Number of Ports] 4\n[Reference]\n50 75 0.01 0.01\n[Number of Frequencies] 1\n[Network Data]\n5.00000 0.60 161.24 0.40 -42.20 0.42 -66.58 0.53 -79.34 !row 1\n0.40 -42.20 0.60 161.20 0.53 -79.34 0.42 -66.58 !row 2\n0.42 -66.58 0.53 -79.34 0.60 161.24 0.40 -42.20 !row 3\n0.53 -79.34 0.42 -66.58 0.40 -42.20 0.60 161.24 !row 4"
        # message = "Incorrect version number."
        # @test_throws message JosephsonCircuits.touchstone_parse(IOBuffer(example))

        # No network data.
        example = "!Example 13 (Version 1.0):\n!2-port S-parameter file, three frequency points\n!freq ReS11 ImS11 ReS21 ImS21 ReS12 ImS12 ReS22 ImS22\n# ghz s ri R 50.0\n"
        message = "No network data."
        @test_throws(
            ErrorException(message),
            JosephsonCircuits.touchstone_parse(IOBuffer(example))
        )

        # No [Two-Port Data Order] keyword in v2 file with 2 ports
        # example = "!Example 17 (Version 2.0):\n!2-port network, S-parameter and noise data\n!Default MA format, GHz frequencies, 50 ohm reference, S-parameters\n[Version] 2.0\n# ghz s ma R 50.0\n[Number of Ports] 2\n[Two-Port Data Order] 21_12\n[Number of Frequencies] 2\n[Number of Noise Frequencies] 2\n[Reference] 50.0 25.0\n[Network Data]\n! freq mags11 angs11 mags21 angs21 mags12 angs12 mags22 angs22 \n2.0 0.95 -26.0 3.57 157.0 0.04 76.0 0.66 -14.0\n22.0 0.6 -144.0 1.3 40.0 0.14 40.0 0.56 -85.0\n[Noise Data]\n4.0 0.7 0.64 69.0 19.0\n18.0 2.7 0.46 -33.0 20.0\n[End]"
        example = "!Example 17 (Version 2.0):\n!2-port network, S-parameter and noise data\n!Default MA format, GHz frequencies, 50 ohm reference, S-parameters\n[Version] 2.0\n# ghz s ma R 50.0\n[Number of Ports] 2\n[Number of Frequencies] 2\n[Number of Noise Frequencies] 2\n[Reference] 50.0 25.0\n[Network Data]\n! freq mags11 angs11 mags21 angs21 mags12 angs12 mags22 angs22 \n2.0 0.95 -26.0 3.57 157.0 0.04 76.0 0.66 -14.0\n22.0 0.6 -144.0 1.3 40.0 0.14 40.0 0.56 -85.0\n[Noise Data]\n4.0 0.7 0.64 69.0 19.0\n18.0 2.7 0.46 -33.0 20.0\n[End]"
        message = "two-port data order must be defined for a v2 file with two ports."
        @test_throws(
            ErrorException(message),
            JosephsonCircuits.touchstone_parse(IOBuffer(example))
        )

        # [Two-Port Data Order] keyword not allowed in v2 file with more than 2 ports
        example = "!Example 4 (Version 2.0):\n! 4-port S-parameter data\n! Default impedance is overridden by [Reference]\n! Data cannot be represented using 1.0 syntax\n! Note that the [Reference] keyword arguments appear on a separate line\n[Version] 2.0\n# GHz S MA R 50\n[Number of Ports] 4\n[Two-Port Data Order] 21_12\n[Reference]\n50 75 0.1 0.1\n[Number of Frequencies] 1\n[Network Data]\n5.00000 0.60 161.24 0.40 -42.20 0.42 -66.58 0.53 -79.34 !row 1\n0.40 -42.20 0.60 161.20 0.53 -79.34 0.42 -66.58 !row 2\n0.42 -66.58 0.53 -79.34 0.60 161.24 0.40 -42.20 !row 3\n0.53 -79.34 0.42 -66.58 0.40 -42.20 0.60 161.24 !row 4"
        message = "[Two-Port Data Order] is only allowed if [Number of Ports] is 2: [two-port data order] 21_12"
        @test_throws(
            ErrorException(message),
            JosephsonCircuits.touchstone_parse(IOBuffer(example))
        )

        # For non-full matrix format, number of ports from data and number of ports from network data not consistent
        # example = "!Example 6 (Version 2.0):\n! 4-port S-parameter data\n! Default impedance is overridden by the [Reference] keyword arguments\n! Note that [Reference] arguments are split across two lines\n! Data cannot be represented using 1.0 syntax\n[Version] 2.0\n# GHz S MA R 50\n[Number of Ports] 4\n[Number of Frequencies] 1\n[Reference] 50 75\n0.01 0.01\n[Matrix Format] Lower\n[Network Data]\n5.00000 0.60 161.24 !row 1\n0.40 -42.20 0.60 161.20 !row 2\n0.42 -66.58 0.53 -79.34 0.60 161.24 !row 3\n0.53 -79.34 0.42 -66.58 0.40 -42.20 0.60 161.24 !row 4\n[End]"
        example = "!Example 6 (Version 2.0):\n! 4-port S-parameter data\n! Default impedance is overridden by the [Reference] keyword arguments\n! Note that [Reference] arguments are split across two lines\n! Data cannot be represented using 1.0 syntax\n[Version] 2.0\n# GHz S MA R 50\n[Number of Ports] 4\n[Number of Frequencies] 1\n[Reference] 50 75\n0.01 0.01\n[Matrix Format] Lower\n[Network Data]\n5.00000 0.60 161.24 !row 1\n0.40 -42.20 0.60 161.20 !row 2\n0.42 -66.58 0.53 -79.34 0.60 161.24 !row 3\n[End]"
        message = "Number of ports not consistent between header and network data."
        @test_throws(
            ErrorException(message),
            JosephsonCircuits.touchstone_parse(IOBuffer(example))
        )

        # in network data, doesn't have expected number of entries based on previous lines.
        example = "!Example 6 (Version 2.0):\n! 4-port S-parameter data\n! Default impedance is overridden by the [Reference] keyword arguments\n! Note that [Reference] arguments are split across two lines\n! Data cannot be represented using 1.0 syntax\n[Version] 2.0\n# GHz S MA R 50\n[Number of Ports] 4\n[Number of Frequencies] 1\n[Reference] 50 75\n0.01 0.01\n[Matrix Format] Lower\n[Network Data]\n5.00000 0.60 !row 1\n0.40 -42.20 0.60 !row 2\n0.42 -66.58 0.53 -79.34 0.60 !row 3\n[End]"
        message = "Number of ports are not consistent between lines:0.42 -66.58 0.53 -79.34 0.60 "
        @test_throws(
            ErrorException(message),
            JosephsonCircuits.touchstone_parse(IOBuffer(example))
        )

    end

    @testset "touchstone_file errors" begin

        message = "The size of the last axis of network data N must equal the number of frequecies"
        @test_throws(
            ErrorException(message),
            JosephsonCircuits.touchstone_file([2.0e6,3.0e6], [0.874020294860635 - 0.18794819544685323im;;;])
        )

        message = "Network data matrix must be square"
        @test_throws(
            ErrorException(message),
            JosephsonCircuits.touchstone_file([2.0e6], zeros(Complex{Float64},1,2,1))
        )

        message = "Version must be 1.0, 1.1, or 2.0"
        @test_throws(
            ErrorException(message),
            JosephsonCircuits.touchstone_file([2.0e6], zeros(Complex{Float64},1,1,1),version=1.5)
        )

        message = "Unknown frequency unit"
        @test_throws(
            ErrorException(message),
            JosephsonCircuits.touchstone_file([2.0e6], zeros(Complex{Float64},1,1,1),frequencyunit="THz")
        )

        message = "Unknown format"
        @test_throws(
            ErrorException(message),
            JosephsonCircuits.touchstone_file([2.0e6], zeros(Complex{Float64},1,1,1),format="DA")
        )

        message = "Unknown parameter"
        @test_throws(
            ErrorException(message),
            JosephsonCircuits.touchstone_file([2.0e6], zeros(Complex{Float64},1,1,1),parameter="K")
        )

        message = "Number of per port impedances must equal number of ports"
        @test_throws(
            ErrorException(message),
            JosephsonCircuits.touchstone_file([2.0e6], zeros(Complex{Float64},1,1,1),reference=[50.0,50.0])
        )

        message = "The port impedances are not equal, so we cannot generate a Touchstone file with version < 2.0."
        @test_throws(
            ErrorException(message),
            JosephsonCircuits.touchstone_file([2.0e6], zeros(Complex{Float64},2,2,1),reference=[50.0,40.0],version=1.0)
        )

        message = "For Touchstone files with version less than 2.0 only the Full matrix format is required."
        @test_throws(
            ErrorException(message),
            JosephsonCircuits.touchstone_file([2.0e6], zeros(Complex{Float64},2,2,1),version=1.0,matrixformat="Lower")
        )

        message = "For Touchstone files with version 2.0 the valid formats are Full, Lower, or Upper."
        @test_throws(
            ErrorException(message),
            JosephsonCircuits.touchstone_file([2.0e6], zeros(Complex{Float64},2,2,1),version=2.0,matrixformat="Triangular")
        )

        message = "Unknown two-port data order."
        @test_throws(
            ErrorException(message),
            JosephsonCircuits.touchstone_file([2.0e6], zeros(Complex{Float64},2,2,1),version=2.0,twoportdataorder="23_32")
        )

        message = "Invalid number of elements in noise data."
        @test_throws(
            ErrorException(message),
            JosephsonCircuits.touchstone_file([2.0e6], zeros(Complex{Float64},2,2,1),version=2.0,noisedata=[1.0,2.0,3.0,4.0,5.0,6.0])
        )

    end

    @testset "arraytonetworkdata errors" begin

        begin
            frequencies = 4.0e9:5.0e8:6.0e9
            N = [0.9546262517670427 - 0.296397700700921im;;; 0.8915960960938982 - 0.44358732281729774im;;; 0.9857309246425359 + 0.046691189499470154im;;; 0.9759591344506418 - 0.21128542054786678im;;; 0.9604441706426364 - 0.2762239892126382im]
            numberofports = 1
            numberoffrequencies = 5
            matrixformat = "Other"
            twoportdataorder = "12_21"
            parameter = "S"
            frequencyunit = "GHz"
            format = "MA"
            R = 50.0
            version = 2.0
            @test_throws(
                ErrorException("Unknown matrixformat."),
                JosephsonCircuits.arraytonetworkdata(frequencies,N, numberofports, numberoffrequencies, 
                    matrixformat, twoportdataorder, parameter, frequencyunit, format, R, version)
            )
        end

        begin
            frequencies = 4.0e9:5.0e8:6.0e9
            N = [0.9546262517670427 - 0.296397700700921im;;; 0.8915960960938982 - 0.44358732281729774im;;; 0.9857309246425359 + 0.046691189499470154im;;; 0.9759591344506418 - 0.21128542054786678im;;; 0.9604441706426364 - 0.2762239892126382im]
            numberofports = 1
            numberoffrequencies = 5
            matrixformat = "Full"
            twoportdataorder = "12_21"
            parameter = "K"
            frequencyunit = "GHz"
            format = "MA"
            R = 50.0
            version = 2.0
            @test_throws(
                ErrorException("Unknown format or version"),
                JosephsonCircuits.arraytonetworkdata(frequencies,N, numberofports, numberoffrequencies, 
                    matrixformat, twoportdataorder, parameter, frequencyunit, format, R, version)
            )
        end
    end

    @testset "networkdatatoarray errors" begin

        networkdata = [4.0, 0.9995813511383583, -17.248815971093425, 4.5, 0.9958480363660398, -26.451285931791276, 5.0, 0.9868361175866559, 2.711906450972103, 5.5, 0.9985678550072272, -12.21545548845392, 6.0, 0.9993761539770525, -16.045248853866596]
        numberofports = 1
        numberoffrequencies = 5
        matrixformat = "Full"
        twoportdataorder = "12_21"
        parameter = "K"
        frequencyunit = "ghz"
        format = "ma"
        R = 50.0
        version = 2.0
        @test_throws(
            ErrorException("Error: Unknown parameter."),
            JosephsonCircuits.networkdatatoarray(networkdata, numberofports,
                numberoffrequencies, matrixformat, twoportdataorder, parameter,
                frequencyunit, format, R, version)
            )

    end

    @testset "touchstone_write" begin

        parsedexamples = Any[
            JosephsonCircuits.TouchstoneFile([2.0e6], [0.874020294860635 - 0.18794819544685323im;;;], "mhz", "s", "ma", 50.0, 1.0, 1, "12_21", 1, 0, [50.0], String[], "Full", Tuple{Char, Vector{Int64}}[], ["Example 1:", "1-port S-parameter file, single frequency point", "freq magS11 angS11"], [2.0, 0.894, -12.136], Float64[]), 
            JosephsonCircuits.TouchstoneFile([2.0e6], [0.8743473504516138 - 0.18801852505875893im;;;], "mhz", "s", "db", 50.0, 1.0, 1, "12_21", 1, 0, [50.0], String[], "Full", Tuple{Char, Vector{Int64}}[], ["Example 1a:", "1-port S-parameter file, single frequency point", "freq magS11 angS11"], [2.0, -0.97, -12.136], Float64[]),
            JosephsonCircuits.TouchstoneFile([1.0e8, 2.0e8, 3.0e8, 4.0e8, 5.0e8], [74.06913073179194 - 5.179418175501303im;;; 55.63103127400724 - 22.47639560495472im;;; 37.494337072416684 - 37.49433707241668im;;; 14.084146883576725 - 26.488427785767808im;;; 0.013089304827962698 - 0.7498857713672935im], "mhz", "z", "ma", 75.0, 1.0, 1, "12_21", 5, 0, [75.0], String[], "Full", Tuple{Char, Vector{Int64}}[], ["Example 2:", "1-port Z-parameter file, multiple frequency points", "freq magZ11 angZ11"], [100.0, 0.99, -4.0, 200.0, 0.8, -22.0, 300.0, 0.707, -45.0, 400.0, 0.4, -62.0, 500.0, 0.01, -89.0], Float64[]),
            JosephsonCircuits.TouchstoneFile([2000.0], [0.8538543439842087 - 0.4164525894496235im 0.009676875823986707 + 0.03881182905103986im; -3.286202326825212 + 1.3949101287067074im 0.6598994788032183 - 0.011518588248607119im;;;], "khz", "h", "ma", 1.0, 1.0, 2, "21_12", 1, 0, [1.0, 1.0], String[], "Full", Tuple{Char, Vector{Int64}}[], ["Example 3:", "2-port H-parameter file, single frequency point", " freq magH11 angH11 magH21 angH21 magH12 angH12 magH22 angH22"], [2.0, 0.95, -26.0, 3.57, 157.0, 0.04, 76.0, 0.66, -1.0], Float64[]),
            JosephsonCircuits.TouchstoneFile([1.0e9, 2.0e9, 1.0e10], [0.3926 - 0.1211im -0.0003 - 0.0021im; -0.0003 - 0.0021im 0.3926 - 0.1211im;;; 0.3517 - 0.3054im -0.0096 - 0.0298im; -0.0096 - 0.0298im 0.3517 - 0.3054im;;; 0.3419 + 0.3336im -0.0134 + 0.0379im; -0.0134 + 0.0379im 0.3419 + 0.3336im], "ghz", "s", "ri", 50.0, 1.0, 2, "21_12", 3, 0, [50.0, 50.0], String[], "Full", Tuple{Char, Vector{Int64}}[], ["Example 4:", "2-port S-parameter file, three frequency points", "freq RelS11 ImS11 ReS21 ImS21 ReS12 ImS12 ReS22 ImS22"], [1.0, 0.3926, -0.1211, -0.0003, -0.0021, -0.0003, -0.0021, 0.3926, -0.1211, 2.0, 0.3517, -0.3054, -0.0096, -0.0298, -0.0096, -0.0298, 0.3517, -0.3054, 10.0, 0.3419, 0.3336, -0.0134, 0.0379, -0.0134, 0.0379, 0.3419, 0.3336], Float64[]),
            JosephsonCircuits.TouchstoneFile([5.0e9, 6.0e9, 7.0e9], [-0.5681244079815996 + 0.1929628385351877im 0.29632183851470006 - 0.2686882357291961im 0.16693665375723588 - 0.38539869438327984im 0.09803970583787712 - 0.5208533537179372im; 0.29632183851470006 - 0.2686882357291961im -0.5679895560694177 + 0.1933594171383067im 0.09803970583787712 - 0.5208533537179372im 0.16693665375723588 - 0.38539869438327984im; 0.16693665375723588 - 0.38539869438327984im 0.09803970583787712 - 0.5208533537179372im -0.5681244079815996 + 0.1929628385351877im 0.29632183851470006 - 0.2686882357291961im; 0.09803970583787712 - 0.5208533537179372im 0.16693665375723588 - 0.38539869438327984im 0.29632183851470006 - 0.2686882357291961im -0.5681244079815996 + 0.1929628385351877im;;; -0.495464624294036 + 0.28180632724119187im 0.286081989392916 - 0.2795659051905141im 0.062441313054034775 - 0.40521732739862937im -0.05730515806890161 - 0.5671120866801361im; 0.286081989392916 - 0.2795659051905141im -0.495464624294036 + 0.28180632724119187im -0.05730515806890161 - 0.5671120866801361im 0.062441313054034775 - 0.40521732739862937im; 0.062441313054034775 - 0.40521732739862937im -0.05730515806890161 - 0.5671120866801361im -0.495464624294036 + 0.28180632724119187im 0.286081989392916 - 0.2795659051905141im; -0.05730515806890161 - 0.5671120866801361im 0.062441313054034775 - 0.40521732739862937im 0.286081989392916 - 0.2795659051905141im -0.495464624294036 + 0.28180632724119187im;;; -0.3638265243449566 + 0.3429726813946975im 0.3102719136297667 - 0.32593149527549903im -0.05845471959176759 - 0.3653533163356367im -0.2540535762162701 - 0.565558821354352im; 0.3102719136297667 - 0.32593149527549903im -0.3638265243449566 + 0.3429726813946975im -0.2540535762162701 - 0.565558821354352im -0.05845471959176759 - 0.3653533163356367im; -0.05845471959176759 - 0.3653533163356367im -0.2540535762162701 - 0.565558821354352im -0.3638265243449566 + 0.3429726813946975im 0.3102719136297667 - 0.32593149527549903im; -0.2540535762162701 - 0.565558821354352im -0.05845471959176759 - 0.3653533163356367im 0.3102719136297667 - 0.32593149527549903im -0.3638265243449566 + 0.3429726813946975im], "ghz", "s", "ma", 50.0, 1.0, 4, "12_21", 3, 0, [50.0, 50.0, 50.0, 50.0], String[], "Full", Tuple{Char, Vector{Int64}}[], ["Example 5:", " 4-port S-parameter data, taken at three frequency points", "row 1", "row 2", "row 3", "row 4", "row 1", "row 2", "row 3", "row 4", "row 1", "row 2", "row 3", "row 4"], [5.0, 0.6, 161.24, 0.4, -42.2, 0.42, -66.58, 0.53, -79.34, 0.4, -42.2, 0.6, 161.2, 0.53, -79.34, 0.42, -66.58, 0.42, -66.58, 0.53, -79.34, 0.6, 161.24, 0.4, -42.2, 0.53, -79.34, 0.42, -66.58, 0.4, -42.2, 0.6, 161.24, 6.0, 0.57, 150.37, 0.4, -44.34, 0.41, -81.24, 0.57, -95.77, 0.4, -44.34, 0.57, 150.37, 0.57, -95.77, 0.41, -81.24, 0.41, -81.24, 0.57, -95.77, 0.57, 150.37, 0.4, -44.34, 0.57, -95.77, 0.41, -81.24, 0.4, -44.34, 0.57, 150.37, 7.0, 0.5, 136.69, 0.45, -46.41, 0.37, -99.09, 0.62, -114.19, 0.45, -46.41, 0.5, 136.69, 0.62, -114.19, 0.37, -99.09, 0.37, -99.09, 0.62, -114.19, 0.5, 136.69, 0.45, -46.41, 0.62, -114.19, 0.37, -99.09, 0.45, -46.41, 0.5, 136.69], Float64[]),
            JosephsonCircuits.TouchstoneFile([2.0e9, 2.2e10], [0.8538543439842087 - 0.4164525894496235im 0.009676875823986707 + 0.03881182905103986im; -3.286202326825212 + 1.3949101287067074im 0.6403951793421577 - 0.1596684510957807im;;; -0.48541019662496837 - 0.35267115137548394im 0.10724622203665693 + 0.0899902653561155im; 0.9958577760546714 + 0.835623892592501im 0.048807215938688565 - 0.5578690309313775im], "ghz", "s", "ma", 50.0, 1.0, 2, "21_12", 2, 2, [50.0, 50.0], String[], "Full", Tuple{Char, Vector{Int64}}[], ["Example 8:", "2-port network, S-parameter and noise data", " NOISE PARAMETERS"], [2.0, 0.95, -26.0, 3.57, 157.0, 0.04, 76.0, 0.66, -14.0, 22.0, 0.6, -144.0, 1.3, 40.0, 0.14, 40.0, 0.56, -85.0], [4.0, 0.7, 0.64, 69.0, 0.38, 18.0, 2.7, 0.46, -33.0, 0.4]),
            JosephsonCircuits.TouchstoneFile([5.0e9], [-0.5681244079815996 + 0.1929628385351877im 0.29632183851470006 - 0.2686882357291961im 0.16693665375723588 - 0.38539869438327984im 0.09803970583787712 - 0.5208533537179372im; 0.29632183851470006 - 0.2686882357291961im -0.5679895560694177 + 0.1933594171383067im 0.09803970583787712 - 0.5208533537179372im 0.16693665375723588 - 0.38539869438327984im; 0.16693665375723588 - 0.38539869438327984im 0.09803970583787712 - 0.5208533537179372im -0.5681244079815996 + 0.1929628385351877im 0.29632183851470006 - 0.2686882357291961im; 0.09803970583787712 - 0.5208533537179372im 0.16693665375723588 - 0.38539869438327984im 0.29632183851470006 - 0.2686882357291961im -0.5681244079815996 + 0.1929628385351877im;;;], "ghz", "s", "ma", 50.0, 2.0, 4, "12_21", 1, 0, [50.0, 75.0, 0.01, 0.01], String[], "Full", Tuple{Char, Vector{Int64}}[], ["Example 4 (Version 2.0):", " 4-port S-parameter data", " Default impedance is overridden by [Reference]", " Data cannot be represented using 1.0 syntax", " Note that the [Reference] keyword arguments appear on a separate line", "row 1", "row 2", "row 3", "row 4"], [5.0, 0.6, 161.24, 0.4, -42.2, 0.42, -66.58, 0.53, -79.34, 0.4, -42.2, 0.6, 161.2, 0.53, -79.34, 0.42, -66.58, 0.42, -66.58, 0.53, -79.34, 0.6, 161.24, 0.4, -42.2, 0.53, -79.34, 0.42, -66.58, 0.4, -42.2, 0.6, 161.24], Float64[]),
            JosephsonCircuits.TouchstoneFile([5.0e9], [-0.5681244079815996 + 0.1929628385351877im 0.29632183851470006 - 0.2686882357291961im 0.16693665375723588 - 0.38539869438327984im 0.09803970583787712 - 0.5208533537179372im; 0.29632183851470006 - 0.2686882357291961im -0.5679895560694177 + 0.1933594171383067im 0.09803970583787712 - 0.5208533537179372im 0.16693665375723588 - 0.38539869438327984im; 0.16693665375723588 - 0.38539869438327984im 0.09803970583787712 - 0.5208533537179372im -0.5681244079815996 + 0.1929628385351877im 0.29632183851470006 - 0.2686882357291961im; 0.09803970583787712 - 0.5208533537179372im 0.16693665375723588 - 0.38539869438327984im 0.29632183851470006 - 0.2686882357291961im -0.5681244079815996 + 0.1929628385351877im;;;], "ghz", "s", "ma", 50.0, 2.0, 4, "12_21", 1, 0, [50.0, 75.0, 0.01, 0.01], String[], "Full", Tuple{Char, Vector{Int64}}[], ["Example 5 (Version 2.0):", " 4-port S-parameter data", " Default impedance is overridden by the [Reference] keyword arguments", " Data cannot be represented using 1.0 syntax", "row 1", "row 2", "row 3", "row 4"], [5.0, 0.6, 161.24, 0.4, -42.2, 0.42, -66.58, 0.53, -79.34, 0.4, -42.2, 0.6, 161.2, 0.53, -79.34, 0.42, -66.58, 0.42, -66.58, 0.53, -79.34, 0.6, 161.24, 0.4, -42.2, 0.53, -79.34, 0.42, -66.58, 0.4, -42.2, 0.6, 161.24], Float64[]),
            JosephsonCircuits.TouchstoneFile([5.0e9], [-0.5681244079815996 + 0.1929628385351877im 0.29632183851470006 - 0.2686882357291961im 0.16693665375723588 - 0.38539869438327984im 0.09803970583787712 - 0.5208533537179372im; 0.29632183851470006 - 0.2686882357291961im -0.5679895560694177 + 0.1933594171383067im 0.09803970583787712 - 0.5208533537179372im 0.16693665375723588 - 0.38539869438327984im; 0.16693665375723588 - 0.38539869438327984im 0.09803970583787712 - 0.5208533537179372im -0.5681244079815996 + 0.1929628385351877im 0.29632183851470006 - 0.2686882357291961im; 0.09803970583787712 - 0.5208533537179372im 0.16693665375723588 - 0.38539869438327984im 0.29632183851470006 - 0.2686882357291961im -0.5681244079815996 + 0.1929628385351877im;;;], "ghz", "s", "ma", 50.0, 2.0, 4, "12_21", 1, 0, [50.0, 75.0, 0.01, 0.01], String[], "Lower", Tuple{Char, Vector{Int64}}[], ["Example 6 (Version 2.0):", " 4-port S-parameter data", " Default impedance is overridden by the [Reference] keyword arguments", " Note that [Reference] arguments are split across two lines", " Data cannot be represented using 1.0 syntax", "row 1", "row 2", "row 3", "row 4"], [5.0, 0.6, 161.24, 0.4, -42.2, 0.6, 161.2, 0.42, -66.58, 0.53, -79.34, 0.6, 161.24, 0.53, -79.34, 0.42, -66.58, 0.4, -42.2, 0.6, 161.24], Float64[]),
            JosephsonCircuits.TouchstoneFile([1.0e8, 2.0e8, 3.0e8, 4.0e8, 5.0e8], [74.06913073179194 - 5.179418175501303im;;; 55.63103127400724 - 22.47639560495472im;;; 37.494337072416684 - 37.49433707241668im;;; 14.084146883576725 - 26.488427785767808im;;; 0.013089304827962698 - 0.7498857713672935im], "mhz", "z", "ma", 50.0, 2.0, 1, "12_21", 5, 0, [20.0], String[], "Full", Tuple{Char, Vector{Int64}}[], ["Example 7 (Version 2.0):", "1-port Z-parameter file, multiple frequency points", "freq magZ11 angZ11"], [100.0, 74.25, -4.0, 200.0, 60.0, -22.0, 300.0, 53.025, -45.0, 400.0, 30.0, -62.0, 500.0, 0.75, -89.0], Float64[]),
            JosephsonCircuits.TouchstoneFile([2.0e6], [0.874020294860635 - 0.18794819544685323im;;;], "mhz", "s", "ma", 50.0, 1.0, 1, "12_21", 1, 0, [50.0], String[], "Full", Tuple{Char, Vector{Int64}}[], ["Example 8 (Version 1.0):", "1-port S-parameter file, single frequency point", "freq magS11 angS11"], [2.0, 0.894, -12.136], Float64[]),
            JosephsonCircuits.TouchstoneFile([1.0e8, 2.0e8, 3.0e8, 4.0e8, 5.0e8], [74.06913073179194 - 5.179418175501303im;;; 55.63103127400724 - 22.47639560495472im;;; 37.494337072416684 - 37.49433707241668im;;; 14.084146883576725 - 26.488427785767808im;;; 0.013089304827962698 - 0.7498857713672935im], "mhz", "z", "ma", 75.0, 1.0, 1, "12_21", 5, 0, [75.0], String[], "Full", Tuple{Char, Vector{Int64}}[], ["Example 9 (Version 1.0):", "1-port Z-parameter file, multiple frequency points", "freq magZ11 angZ11"], [100.0, 0.99, -4.0, 200.0, 0.8, -22.0, 300.0, 0.707, -45.0, 400.0, 0.4, -62.0, 500.0, 0.01, -89.0], Float64[]),
            JosephsonCircuits.TouchstoneFile([1.0e8, 2.0e8, 3.0e8, 4.0e8, 5.0e8], [74.06913073179194 - 5.179418175501303im;;; 55.63103127400724 - 22.47639560495472im;;; 37.494337072416684 - 37.49433707241668im;;; 14.084146883576725 - 26.488427785767808im;;; 0.013089304827962698 - 0.7498857713672935im], "mhz", "z", "ma", 50.0, 2.0, 1, "12_21", 5, 0, [20.0], String[], "Full", Tuple{Char, Vector{Int64}}[], ["Example 10 (Version 2.0):", "1-port Z-parameter file, multiple frequency points", "freq magZ11 angZ11"], [100.0, 74.25, -4.0, 200.0, 60.0, -22.0, 300.0, 53.025, -45.0, 400.0, 30.0, -62.0, 500.0, 0.75, -89.0], Float64[]),
            JosephsonCircuits.TouchstoneFile([2000.0], [0.8538543439842087 - 0.4164525894496235im 0.009676875823986707 + 0.03881182905103986im; -3.286202326825212 + 1.3949101287067074im 0.6403951793421577 - 0.1596684510957807im;;;], "khz", "h", "ma", 1.0, 1.0, 2, "21_12", 1, 0, [1.0, 1.0], String[], "Full", Tuple{Char, Vector{Int64}}[], ["Example 11 (Version 1.0):", "2-port H-parameter file, single frequency point", " freq magH11 angH11 magH21 angH21 magH12 angH12 magH22 angH22"], [2.0, 0.95, -26.0, 3.57, 157.0, 0.04, 76.0, 0.66, -14.0], Float64[]),
            JosephsonCircuits.TouchstoneFile([2000.0], [0.8538543439842087 - 0.4164525894496235im 0.009676875823986707 + 0.03881182905103986im; -3.286202326825212 + 1.3949101287067074im 0.6403951793421577 - 0.1596684510957807im;;;], "khz", "h", "ma", 1.0, 2.0, 2, "21_12", 1, 0, [1.0, 1.0], String[], "Full", Tuple{Char, Vector{Int64}}[], ["Example 12 (Version 2.0):", "2-port H-parameter file, single frequency point", " freq magH11 angH11 magH21 angH21 magH12 angH12 magH22 angH22"], [2.0, 0.95, -26.0, 3.57, 157.0, 0.04, 76.0, 0.66, -14.0], Float64[]),
            JosephsonCircuits.TouchstoneFile([1.0e9, 2.0e9, 1.0e10], [0.3926 - 0.1211im -0.0003 - 0.0021im; -0.0003 - 0.0021im 0.3926 - 0.1211im;;; 0.3517 - 0.3054im -0.0096 - 0.0298im; -0.0096 - 0.0298im 0.3517 - 0.3054im;;; 0.3419 + 0.3336im -0.0134 + 0.0379im; -0.0134 + 0.0379im 0.3419 + 0.3336im], "ghz", "s", "ri", 50.0, 1.0, 2, "21_12", 3, 0, [50.0, 50.0], String[], "Full", Tuple{Char, Vector{Int64}}[], ["Example 13 (Version 1.0):", "2-port S-parameter file, three frequency points", "freq ReS11 ImS11 ReS21 ImS21 ReS12 ImS12 ReS22 ImS22"], [1.0, 0.3926, -0.1211, -0.0003, -0.0021, -0.0003, -0.0021, 0.3926, -0.1211, 2.0, 0.3517, -0.3054, -0.0096, -0.0298, -0.0096, -0.0298, 0.3517, -0.3054, 10.0, 0.3419, 0.3336, -0.0134, 0.0379, -0.0134, 0.0379, 0.3419, 0.3336], Float64[]),
            JosephsonCircuits.TouchstoneFile([5.0e9, 6.0e9, 7.0e9], [-0.5681244079815996 + 0.1929628385351877im 0.29632183851470006 - 0.2686882357291961im 0.16693665375723588 - 0.38539869438327984im 0.09803970583787712 - 0.5208533537179372im; 0.29632183851470006 - 0.2686882357291961im -0.5679895560694177 + 0.1933594171383067im 0.09803970583787712 - 0.5208533537179372im 0.16693665375723588 - 0.38539869438327984im; 0.16693665375723588 - 0.38539869438327984im 0.09803970583787712 - 0.5208533537179372im -0.5681244079815996 + 0.1929628385351877im 0.29632183851470006 - 0.2686882357291961im; 0.09803970583787712 - 0.5208533537179372im 0.16693665375723588 - 0.38539869438327984im 0.29632183851470006 - 0.2686882357291961im -0.5681244079815996 + 0.1929628385351877im;;; -0.495464624294036 + 0.28180632724119187im 0.286081989392916 - 0.2795659051905141im 0.062441313054034775 - 0.40521732739862937im -0.05730515806890161 - 0.5671120866801361im; 0.286081989392916 - 0.2795659051905141im -0.495464624294036 + 0.28180632724119187im -0.05730515806890161 - 0.5671120866801361im 0.062441313054034775 - 0.40521732739862937im; 0.062441313054034775 - 0.40521732739862937im -0.05730515806890161 - 0.5671120866801361im -0.495464624294036 + 0.28180632724119187im 0.286081989392916 - 0.2795659051905141im; -0.05730515806890161 - 0.5671120866801361im 0.062441313054034775 - 0.40521732739862937im 0.286081989392916 - 0.2795659051905141im -0.495464624294036 + 0.28180632724119187im;;; -0.3638265243449566 + 0.3429726813946975im 0.3102719136297667 - 0.32593149527549903im -0.05845471959176759 - 0.3653533163356367im -0.2540535762162701 - 0.565558821354352im; 0.3102719136297667 - 0.32593149527549903im -0.3638265243449566 + 0.3429726813946975im -0.2540535762162701 - 0.565558821354352im -0.05845471959176759 - 0.3653533163356367im; -0.05845471959176759 - 0.3653533163356367im -0.2540535762162701 - 0.565558821354352im -0.3638265243449566 + 0.3429726813946975im 0.3102719136297667 - 0.32593149527549903im; -0.2540535762162701 - 0.565558821354352im -0.05845471959176759 - 0.3653533163356367im 0.3102719136297667 - 0.32593149527549903im -0.3638265243449566 + 0.3429726813946975im], "ghz", "s", "ma", 50.0, 1.0, 4, "12_21", 3, 0, [50.0, 50.0, 50.0, 50.0], String[], "Full", Tuple{Char, Vector{Int64}}[], ["Example 14 (Version 1.0):", " 4-port S-parameter data, taken at three frequency points", " note that data points need not be aligned", "row 1", "row 2", "row 3", "row 4", "row 1", "row 2", "row 3", "row 4", "row 1", "row 2", "row 3", "row 4"], [5.0, 0.6, 161.24, 0.4, -42.2, 0.42, -66.58, 0.53, -79.34, 0.4, -42.2, 0.6, 161.2, 0.53, -79.34, 0.42, -66.58, 0.42, -66.58, 0.53, -79.34, 0.6, 161.24, 0.4, -42.2, 0.53, -79.34, 0.42, -66.58, 0.4, -42.2, 0.6, 161.24, 6.0, 0.57, 150.37, 0.4, -44.34, 0.41, -81.24, 0.57, -95.77, 0.4, -44.34, 0.57, 150.37, 0.57, -95.77, 0.41, -81.24, 0.41, -81.24, 0.57, -95.77, 0.57, 150.37, 0.4, -44.34, 0.57, -95.77, 0.41, -81.24, 0.4, -44.34, 0.57, 150.37, 7.0, 0.5, 136.69, 0.45, -46.41, 0.37, -99.09, 0.62, -114.19, 0.45, -46.41, 0.5, 136.69, 0.62, -114.19, 0.37, -99.09, 0.37, -99.09, 0.62, -114.19, 0.5, 136.69, 0.45, -46.41, 0.62, -114.19, 0.37, -99.09, 0.45, -46.41, 0.5, 136.69], Float64[]), 
            JosephsonCircuits.TouchstoneFile([5.0e6], [8.0 + 9.0im 2.0 - 1.0im 3.0 - 2.0im 1.0 + 3.0im 1.0 + 0.1im 0.2 - 0.2im; 2.0 - 1.0im 7.0 + 7.0im 1.8 - 2.0im -1.0 - 1.0im -0.5 + 0.5im 0.2 - 0.1im; 3.0 - 2.0im 1.8 - 2.0im 5.8 + 6.0im 1.2 + 0.8im 0.9 + 0.7im 0.3 - 0.5im; 1.0 + 3.0im -1.0 - 1.0im 1.2 + 0.8im 6.3 + 8.0im 2.0 - 0.5im 1.5 + 0.6im; 1.0 + 0.1im -0.5 + 0.5im 0.9 + 0.7im 2.0 - 0.5im 4.7 - 6.0im -1.0 + 2.0im; 0.2 - 0.2im 0.2 - 0.1im 0.3 - 0.5im 1.5 + 0.6im -1.0 + 2.0im 5.5 - 7.0im;;;], "mhz", "y", "ri", 50.0, 2.0, 6, "12_21", 1, 0, [50.0, 75.0, 75.0, 50.0, 0.01, 0.01], String[], "Full", [('D', [2, 3]), ('D', [6, 5]), ('C', [2, 3]), ('C', [6, 5]), ('S', [4]), ('S', [1])], ["Example 16 (Version 2.0):", "6-port component shown; note that all five ports are used in some", "relationship"], [5.0, 8.0, 9.0, 2.0, -1.0, 3.0, -2.0, 1.0, 3.0, 1.0, 0.1, 0.2, -0.2, 2.0, -1.0, 7.0, 7.0, 1.8, -2.0, -1.0, -1.0, -0.5, 0.5, 0.2, -0.1, 3.0, -2.0, 1.8, -2.0, 5.8, 6.0, 1.2, 0.8, 0.9, 0.7, 0.3, -0.5, 1.0, 3.0, -1.0, -1.0, 1.2, 0.8, 6.3, 8.0, 2.0, -0.5, 1.5, 0.6, 1.0, 0.1, -0.5, 0.5, 0.9, 0.7, 2.0, -0.5, 4.7, -6.0, -1.0, 2.0, 0.2, -0.2, 0.2, -0.1, 0.3, -0.5, 1.5, 0.6, -1.0, 2.0, 5.5, -7.0], Float64[]),
            JosephsonCircuits.TouchstoneFile([2.0e9, 2.2e10], [0.8538543439842087 - 0.4164525894496235im 0.009676875823986707 + 0.03881182905103986im; -3.286202326825212 + 1.3949101287067074im 0.6403951793421577 - 0.1596684510957807im;;; -0.48541019662496837 - 0.35267115137548394im 0.10724622203665693 + 0.0899902653561155im; 0.9958577760546714 + 0.835623892592501im 0.048807215938688565 - 0.5578690309313775im], "ghz", "s", "ma", 50.0, 2.0, 2, "21_12", 2, 2, [50.0, 25.0], String[], "Full", Tuple{Char, Vector{Int64}}[], ["Example 17 (Version 2.0):", "2-port network, S-parameter and noise data", "Default MA format, GHz frequencies, 50 ohm reference, S-parameters"], [2.0, 0.95, -26.0, 3.57, 157.0, 0.04, 76.0, 0.66, -14.0, 22.0, 0.6, -144.0, 1.3, 40.0, 0.14, 40.0, 0.56, -85.0], [4.0, 0.7, 0.64, 69.0, 19.0, 18.0, 2.7, 0.46, -33.0, 20.0]),
            JosephsonCircuits.TouchstoneFile([2.0e9, 2.2e10], [0.8538543439842087 - 0.4164525894496235im 0.009676875823986707 + 0.03881182905103986im; -3.286202326825212 + 1.3949101287067074im 0.6403951793421577 - 0.1596684510957807im;;; -0.48541019662496837 - 0.35267115137548394im 0.10724622203665693 + 0.0899902653561155im; 0.9958577760546714 + 0.835623892592501im 0.048807215938688565 - 0.5578690309313775im], "ghz", "s", "ma", 50.0, 1.0, 2, "21_12", 2, 2, [50.0, 50.0], String[], "Full", Tuple{Char, Vector{Int64}}[], ["Example 18 (Version 1.0):", "2-port network, S-parameter and noise data", "Default MA format, GHz frequencies, 50 ohm reference, S-parameters", " NETWORK PARAMETERS", " NOISE PARAMETERS"], [2.0, 0.95, -26.0, 3.57, 157.0, 0.04, 76.0, 0.66, -14.0, 22.0, 0.6, -144.0, 1.3, 40.0, 0.14, 40.0, 0.56, -85.0], [4.0, 0.7, 0.64, 69.0, 0.38, 18.0, 2.7, 0.46, -33.0, 0.4]),
            JosephsonCircuits.TouchstoneFile([2.0e9, 2.2e10], [0.8538543439842087 - 0.4164525894496235im 0.009676875823986707 + 0.03881182905103986im; -3.286202326825212 + 1.3949101287067074im 0.6403951793421577 - 0.1596684510957807im;;; -0.48541019662496837 - 0.35267115137548394im 0.10724622203665693 + 0.0899902653561155im; 0.9958577760546714 + 0.835623892592501im 0.048807215938688565 - 0.5578690309313775im], "ghz", "s", "ma", 50.0, 2.0, 2, "21_12", 2, 2, [50.0, 25.0], String[], "Full", Tuple{Char, Vector{Int64}}[], ["Example 19 (Version 2.0):", "2-port network, S-parameter and noise data", "Default MA format, GHz frequencies, 50 ohm reference, S-parameters", " NETWORK PARAMETERS", " NOISE PARAMETERS"], [2.0, 0.95, -26.0, 3.57, 157.0, 0.04, 76.0, 0.66, -14.0, 22.0, 0.6, -144.0, 1.3, 40.0, 0.14, 40.0, 0.56, -85.0], [4.0, 0.7, 0.64, 69.0, 19.0, 18.0, 2.7, 0.46, -33.0, 20.0]),
            JosephsonCircuits.TouchstoneFile([2.0e9, 2.2e10], [0.8538543439842087 - 0.4164525894496235im 0.009676875823986707 + 0.03881182905103986im; -3.286202326825212 + 1.3949101287067074im 0.6403951793421577 - 0.1596684510957807im;;; -0.48541019662496837 - 0.35267115137548394im 0.10724622203665693 + 0.0899902653561155im; 0.9958577760546714 + 0.835623892592501im 0.048807215938688565 - 0.5578690309313775im], "ghz", "s", "ma", 50.0, 2.0, 2, "21_12", 2, 2, [50.0, 25.0], String["This is information."], "Full", Tuple{Char, Vector{Int64}}[], ["Example 19 (Version 2.0):", "2-port network, S-parameter and noise data", "Default MA format, GHz frequencies, 50 ohm reference, S-parameters", " NETWORK PARAMETERS", " NOISE PARAMETERS"], [2.0, 0.95, -26.0, 3.57, 157.0, 0.04, 76.0, 0.66, -14.0, 22.0, 0.6, -144.0, 1.3, 40.0, 0.14, 40.0, 0.56, -85.0], [4.0, 0.7, 0.64, 69.0, 19.0, 18.0, 2.7, 0.46, -33.0, 20.0]),
        ]

        writtenexamples = String[
            "!Example 1:\n!1-port S-parameter file, single frequency point\n!freq magS11 angS11\n# mhz s ma R 50.0\n! freq mags11 angs11 \n2.0 0.894 -12.136\n",
            "!Example 1a:\n!1-port S-parameter file, single frequency point\n!freq magS11 angS11\n# mhz s db R 50.0\n! freq logmags11 angs11 \n2.0 -0.97 -12.136\n",
            "!Example 2:\n!1-port Z-parameter file, multiple frequency points\n!freq magZ11 angZ11\n# mhz z ma R 75.0\n! freq magz11 angz11 \n100.0 0.99 -4.0\n200.0 0.8 -22.0\n300.0 0.707 -45.0\n400.0 0.4 -62.0\n500.0 0.01 -89.0\n",
            "!Example 3:\n!2-port H-parameter file, single frequency point\n! freq magH11 angH11 magH21 angH21 magH12 angH12 magH22 angH22\n# khz h ma R 1.0\n! freq magh11 angh11 magh21 angh21 magh12 angh12 magh22 angh22 \n2.0 0.95 -26.0 3.57 157.0 0.04 76.0 0.66 -1.0\n",
            "!Example 4:\n!2-port S-parameter file, three frequency points\n!freq RelS11 ImS11 ReS21 ImS21 ReS12 ImS12 ReS22 ImS22\n# ghz s ri R 50.0\n! freq Res11 Ims11 Res21 Ims21 Res12 Ims12 Res22 Ims22 \n1.0 0.3926 -0.1211 -0.0003 -0.0021 -0.0003 -0.0021 0.3926 -0.1211\n2.0 0.3517 -0.3054 -0.0096 -0.0298 -0.0096 -0.0298 0.3517 -0.3054\n10.0 0.3419 0.3336 -0.0134 0.0379 -0.0134 0.0379 0.3419 0.3336\n",
            "!Example 5:\n! 4-port S-parameter data, taken at three frequency points\n!row 1\n!row 2\n!row 3\n!row 4\n!row 1\n!row 2\n!row 3\n!row 4\n!row 1\n!row 2\n!row 3\n!row 4\n# ghz s ma R 50.0\n! freq mags11 angs11 mags12 angs12 mags13 angs13 mags14 angs14 mags21 angs21 mags22 angs22 mags23 angs23 mags24 angs24 mags31 angs31 mags32 angs32 mags33 angs33 mags34 angs34 mags41 angs41 mags42 angs42 mags43 angs43 mags44 angs44 \n5.0 0.6 161.24 0.4 -42.2 0.42 -66.58 0.53 -79.34\n 0.4 -42.2 0.6 161.2 0.53 -79.34\n 0.42 -66.58\n 0.42 -66.58 0.53 -79.34 0.6 161.24 0.4 -42.2\n 0.53 -79.34 0.42 -66.58 0.4 -42.2\n 0.6 161.24\n6.0 0.57 150.37 0.4 -44.34 0.41 -81.24 0.57 -95.77\n 0.4 -44.34 0.57 150.37 0.57 -95.77\n 0.41 -81.24\n 0.41 -81.24 0.57 -95.77 0.57 150.37 0.4 -44.34\n 0.57 -95.77 0.41 -81.24 0.4 -44.34\n 0.57 150.37\n7.0 0.5 136.69 0.45 -46.41 0.37 -99.09 0.62 -114.19\n 0.45 -46.41 0.5 136.69 0.62 -114.19\n 0.37 -99.09\n 0.37 -99.09 0.62 -114.19 0.5 136.69 0.45 -46.41\n 0.62 -114.19 0.37 -99.09 0.45 -46.41\n 0.5 136.69\n",
            "!Example 8:\n!2-port network, S-parameter and noise data\n! NOISE PARAMETERS\n# ghz s ma R 50.0\n! freq mags11 angs11 mags21 angs21 mags12 angs12 mags22 angs22 \n2.0 0.95 -26.0 3.57 157.0 0.04 76.0 0.66 -14.0\n22.0 0.6 -144.0 1.3 40.0 0.14 40.0 0.56 -85.0\n[Noise Data]\n4.0 0.7 0.64 69.0 0.38\n18.0 2.7 0.46 -33.0 0.4\n",
            "!Example 4 (Version 2.0):\n! 4-port S-parameter data\n! Default impedance is overridden by [Reference]\n! Data cannot be represented using 1.0 syntax\n! Note that the [Reference] keyword arguments appear on a separate line\n!row 1\n!row 2\n!row 3\n!row 4\n[Version] 2.0\n# ghz s ma R 50.0\n[Number of Ports] 4\n[Number of Frequencies] 1\n[Reference] 50.0 75.0 0.01 0.01\n[Network Data]\n! freq mags11 angs11 mags12 angs12 mags13 angs13 mags14 angs14 mags21 angs21 mags22 angs22 mags23 angs23 mags24 angs24 mags31 angs31 mags32 angs32 mags33 angs33 mags34 angs34 mags41 angs41 mags42 angs42 mags43 angs43 mags44 angs44 \n5.0 0.6 161.24 0.4 -42.2 0.42 -66.58 0.53 -79.34\n 0.4 -42.2 0.6 161.2 0.53 -79.34\n 0.42 -66.58\n 0.42 -66.58 0.53 -79.34 0.6 161.24 0.4 -42.2\n 0.53 -79.34 0.42 -66.58 0.4 -42.2\n 0.6 161.24\n[End]",
            "!Example 5 (Version 2.0):\n! 4-port S-parameter data\n! Default impedance is overridden by the [Reference] keyword arguments\n! Data cannot be represented using 1.0 syntax\n!row 1\n!row 2\n!row 3\n!row 4\n[Version] 2.0\n# ghz s ma R 50.0\n[Number of Ports] 4\n[Number of Frequencies] 1\n[Reference] 50.0 75.0 0.01 0.01\n[Network Data]\n! freq mags11 angs11 mags12 angs12 mags13 angs13 mags14 angs14 mags21 angs21 mags22 angs22 mags23 angs23 mags24 angs24 mags31 angs31 mags32 angs32 mags33 angs33 mags34 angs34 mags41 angs41 mags42 angs42 mags43 angs43 mags44 angs44 \n5.0 0.6 161.24 0.4 -42.2 0.42 -66.58 0.53 -79.34\n 0.4 -42.2 0.6 161.2 0.53 -79.34\n 0.42 -66.58\n 0.42 -66.58 0.53 -79.34 0.6 161.24 0.4 -42.2\n 0.53 -79.34 0.42 -66.58 0.4 -42.2\n 0.6 161.24\n[End]",
            "!Example 6 (Version 2.0):\n! 4-port S-parameter data\n! Default impedance is overridden by the [Reference] keyword arguments\n! Note that [Reference] arguments are split across two lines\n! Data cannot be represented using 1.0 syntax\n!row 1\n!row 2\n!row 3\n!row 4\n[Version] 2.0\n# ghz s ma R 50.0\n[Number of Ports] 4\n[Number of Frequencies] 1\n[Reference] 50.0 75.0 0.01 0.01\n[Matrix Format] Lower\n[Network Data]\n! freq mags11 angs11 mags21 angs21 mags22 angs22 mags31 angs31 mags32 angs32 mags33 angs33 mags41 angs41 mags42 angs42 mags43 angs43 mags44 angs44 \n5.0 0.6 161.24 0.4 -42.2 0.6 161.2 0.42 -66.58\n 0.53 -79.34 0.6 161.24 0.53 -79.34\n 0.42 -66.58\n 0.4 -42.2 0.6 161.24\n[End]",
            "!Example 7 (Version 2.0):\n!1-port Z-parameter file, multiple frequency points\n!freq magZ11 angZ11\n[Version] 2.0\n# mhz z ma R 50.0\n[Number of Ports] 1\n[Number of Frequencies] 5\n[Reference] 20.0\n[Network Data]\n! freq magz11 angz11 \n100.0 74.25 -4.0\n200.0 60.0 -22.0\n300.0 53.025 -45.0\n400.0 30.0 -62.0\n500.0 0.75 -89.0\n[End]",
            "!Example 8 (Version 1.0):\n!1-port S-parameter file, single frequency point\n!freq magS11 angS11\n# mhz s ma R 50.0\n! freq mags11 angs11 \n2.0 0.894 -12.136\n",
            "!Example 9 (Version 1.0):\n!1-port Z-parameter file, multiple frequency points\n!freq magZ11 angZ11\n# mhz z ma R 75.0\n! freq magz11 angz11 \n100.0 0.99 -4.0\n200.0 0.8 -22.0\n300.0 0.707 -45.0\n400.0 0.4 -62.0\n500.0 0.01 -89.0\n",
            "!Example 10 (Version 2.0):\n!1-port Z-parameter file, multiple frequency points\n!freq magZ11 angZ11\n[Version] 2.0\n# mhz z ma R 50.0\n[Number of Ports] 1\n[Number of Frequencies] 5\n[Reference] 20.0\n[Network Data]\n! freq magz11 angz11 \n100.0 74.25 -4.0\n200.0 60.0 -22.0\n300.0 53.025 -45.0\n400.0 30.0 -62.0\n500.0 0.75 -89.0\n[End]",
            "!Example 11 (Version 1.0):\n!2-port H-parameter file, single frequency point\n! freq magH11 angH11 magH21 angH21 magH12 angH12 magH22 angH22\n# khz h ma R 1.0\n! freq magh11 angh11 magh21 angh21 magh12 angh12 magh22 angh22 \n2.0 0.95 -26.0 3.57 157.0 0.04 76.0 0.66 -14.0\n",
            "!Example 12 (Version 2.0):\n!2-port H-parameter file, single frequency point\n! freq magH11 angH11 magH21 angH21 magH12 angH12 magH22 angH22\n[Version] 2.0\n# khz h ma R 1.0\n[Number of Ports] 2\n[Two-Port Data Order] 21_12\n[Number of Frequencies] 1\n[Reference] 1.0 1.0\n[Network Data]\n! freq magh11 angh11 magh21 angh21 magh12 angh12 magh22 angh22 \n2.0 0.95 -26.0 3.57 157.0 0.04 76.0 0.66 -14.0\n[End]",
            "!Example 13 (Version 1.0):\n!2-port S-parameter file, three frequency points\n!freq ReS11 ImS11 ReS21 ImS21 ReS12 ImS12 ReS22 ImS22\n# ghz s ri R 50.0\n! freq Res11 Ims11 Res21 Ims21 Res12 Ims12 Res22 Ims22 \n1.0 0.3926 -0.1211 -0.0003 -0.0021 -0.0003 -0.0021 0.3926 -0.1211\n2.0 0.3517 -0.3054 -0.0096 -0.0298 -0.0096 -0.0298 0.3517 -0.3054\n10.0 0.3419 0.3336 -0.0134 0.0379 -0.0134 0.0379 0.3419 0.3336\n",
            "!Example 14 (Version 1.0):\n! 4-port S-parameter data, taken at three frequency points\n! note that data points need not be aligned\n!row 1\n!row 2\n!row 3\n!row 4\n!row 1\n!row 2\n!row 3\n!row 4\n!row 1\n!row 2\n!row 3\n!row 4\n# ghz s ma R 50.0\n! freq mags11 angs11 mags12 angs12 mags13 angs13 mags14 angs14 mags21 angs21 mags22 angs22 mags23 angs23 mags24 angs24 mags31 angs31 mags32 angs32 mags33 angs33 mags34 angs34 mags41 angs41 mags42 angs42 mags43 angs43 mags44 angs44 \n5.0 0.6 161.24 0.4 -42.2 0.42 -66.58 0.53 -79.34\n 0.4 -42.2 0.6 161.2 0.53 -79.34\n 0.42 -66.58\n 0.42 -66.58 0.53 -79.34 0.6 161.24 0.4 -42.2\n 0.53 -79.34 0.42 -66.58 0.4 -42.2\n 0.6 161.24\n6.0 0.57 150.37 0.4 -44.34 0.41 -81.24 0.57 -95.77\n 0.4 -44.34 0.57 150.37 0.57 -95.77\n 0.41 -81.24\n 0.41 -81.24 0.57 -95.77 0.57 150.37 0.4 -44.34\n 0.57 -95.77 0.41 -81.24 0.4 -44.34\n 0.57 150.37\n7.0 0.5 136.69 0.45 -46.41 0.37 -99.09 0.62 -114.19\n 0.45 -46.41 0.5 136.69 0.62 -114.19\n 0.37 -99.09\n 0.37 -99.09 0.62 -114.19 0.5 136.69 0.45 -46.41\n 0.62 -114.19 0.37 -99.09 0.45 -46.41\n 0.5 136.69\n",
            # the file below doesn't have the mixed mode order data. need to support writing that.
            "!Example 16 (Version 2.0):\n!6-port component shown; note that all five ports are used in some\n!relationship\n[Version] 2.0\n# mhz y ri R 50.0\n[Number of Ports] 6\n[Number of Frequencies] 1\n[Reference] 50.0 75.0 75.0 50.0 0.01 0.01\n[Network Data]\n! freq Rey11 Imy11 Rey12 Imy12 Rey13 Imy13 Rey14 Imy14 Rey15 Imy15 Rey16 Imy16 Rey21 Imy21 Rey22 Imy22 Rey23 Imy23 Rey24 Imy24 Rey25 Imy25 Rey26 Imy26 Rey31 Imy31 Rey32 Imy32 Rey33 Imy33 Rey34 Imy34 Rey35 Imy35 Rey36 Imy36 Rey41 Imy41 Rey42 Imy42 Rey43 Imy43 Rey44 Imy44 Rey45 Imy45 Rey46 Imy46 Rey51 Imy51 Rey52 Imy52 Rey53 Imy53 Rey54 Imy54 Rey55 Imy55 Rey56 Imy56 Rey61 Imy61 Rey62 Imy62 Rey63 Imy63 Rey64 Imy64 Rey65 Imy65 Rey66 Imy66 \n5.0 8.0 9.0 2.0 -1.0 3.0 -2.0 1.0 3.0\n 1.0 0.1 0.2 -0.2\n 2.0 -1.0 7.0 7.0 1.8 -2.0 -1.0 -1.0\n -0.5 0.5 0.2 -0.1\n 3.0 -2.0 1.8 -2.0 5.8 6.0 1.2 0.8\n 0.9 0.7 0.3 -0.5\n 1.0 3.0 -1.0 -1.0 1.2 0.8 6.3 8.0\n 2.0 -0.5 1.5 0.6\n 1.0 0.1 -0.5 0.5 0.9 0.7 2.0 -0.5\n 4.7 -6.0 -1.0 2.0\n 0.2 -0.2 0.2 -0.1 0.3 -0.5 1.5 0.6\n -1.0 2.0 5.5 -7.0\n[End]",
            "!Example 17 (Version 2.0):\n!2-port network, S-parameter and noise data\n!Default MA format, GHz frequencies, 50 ohm reference, S-parameters\n[Version] 2.0\n# ghz s ma R 50.0\n[Number of Ports] 2\n[Two-Port Data Order] 21_12\n[Number of Frequencies] 2\n[Number of Noise Frequencies] 2\n[Reference] 50.0 25.0\n[Network Data]\n! freq mags11 angs11 mags21 angs21 mags12 angs12 mags22 angs22 \n2.0 0.95 -26.0 3.57 157.0 0.04 76.0 0.66 -14.0\n22.0 0.6 -144.0 1.3 40.0 0.14 40.0 0.56 -85.0\n[Noise Data]\n4.0 0.7 0.64 69.0 19.0\n18.0 2.7 0.46 -33.0 20.0\n[End]",
            "!Example 18 (Version 1.0):\n!2-port network, S-parameter and noise data\n!Default MA format, GHz frequencies, 50 ohm reference, S-parameters\n! NETWORK PARAMETERS\n! NOISE PARAMETERS\n# ghz s ma R 50.0\n! freq mags11 angs11 mags21 angs21 mags12 angs12 mags22 angs22 \n2.0 0.95 -26.0 3.57 157.0 0.04 76.0 0.66 -14.0\n22.0 0.6 -144.0 1.3 40.0 0.14 40.0 0.56 -85.0\n[Noise Data]\n4.0 0.7 0.64 69.0 0.38\n18.0 2.7 0.46 -33.0 0.4\n",
            "!Example 19 (Version 2.0):\n!2-port network, S-parameter and noise data\n!Default MA format, GHz frequencies, 50 ohm reference, S-parameters\n! NETWORK PARAMETERS\n! NOISE PARAMETERS\n[Version] 2.0\n# ghz s ma R 50.0\n[Number of Ports] 2\n[Two-Port Data Order] 21_12\n[Number of Frequencies] 2\n[Number of Noise Frequencies] 2\n[Reference] 50.0 25.0\n[Network Data]\n! freq mags11 angs11 mags21 angs21 mags12 angs12 mags22 angs22 \n2.0 0.95 -26.0 3.57 157.0 0.04 76.0 0.66 -14.0\n22.0 0.6 -144.0 1.3 40.0 0.14 40.0 0.56 -85.0\n[Noise Data]\n4.0 0.7 0.64 69.0 19.0\n18.0 2.7 0.46 -33.0 20.0\n[End]",
            "!Example 19 (Version 2.0):\n!2-port network, S-parameter and noise data\n!Default MA format, GHz frequencies, 50 ohm reference, S-parameters\n! NETWORK PARAMETERS\n! NOISE PARAMETERS\n[Version] 2.0\n# ghz s ma R 50.0\n[Number of Ports] 2\n[Begin Information] This is information.\n[End Information][Two-Port Data Order] 21_12\n[Number of Frequencies] 2\n[Number of Noise Frequencies] 2\n[Reference] 50.0 25.0\n[Network Data]\n! freq mags11 angs11 mags21 angs21 mags12 angs12 mags22 angs22 \n2.0 0.95 -26.0 3.57 157.0 0.04 76.0 0.66 -14.0\n22.0 0.6 -144.0 1.3 40.0 0.14 40.0 0.56 -85.0\n[Noise Data]\n4.0 0.7 0.64 69.0 19.0\n18.0 2.7 0.46 -33.0 20.0\n[End]",
        ]

        for i in 1:length(parsedexamples)
            io = IOBuffer();
            JosephsonCircuits.touchstone_write(io, parsedexamples[i]);
            # @show String(take!(io));
            @test String(take!(io)) == writtenexamples[i]
        end
    end

    @testset "touchstone_write errors" begin

        message = "Number of mixed mode order descriptors must match the number of ports."
        io = IOBuffer();
        @test_throws(
            ErrorException(message),
            JosephsonCircuits.touchstone_write(io,[2.0e6], zeros(Complex{Float64},1,1,1),mixedmodeorder = [('D',[1,2]),('C',[1,2])])
        )

        message = "Version 1.0 or 1.1 files do not support mixed-mode order (common and differential modes. Save as a Version 2.0 file or delete mixed-mode order data."
        io = IOBuffer();
        @test_throws(
            ErrorException(message),
            JosephsonCircuits.touchstone_write(io,[2.0e6], zeros(Complex{Float64},2,2,1),mixedmodeorder = [('D',[1,2]),('C',[1,2])],version=1.0)
        )

        # message = "Unknown format"
        # io = IOBuffer();
        # @test_throws message JosephsonCircuits.touchstone_write(io,[2.0e6], zeros(Complex{Float64},2,2,1),format="DA")

        message = "Unknown parameter"
        io = IOBuffer();
        @test_throws(
            ErrorException(message),
            JosephsonCircuits.touchstone_write(io,JosephsonCircuits.TouchstoneFile([2.0e6], [0.874020294860635 - 0.18794819544685323im;;;], "mhz", "k", "ma", 50.0, 1.0, 1, "12_21", 1, 0, [50.0], String[], "Full", Tuple{Char, Vector{Int64}}[], ["Example 1:", "1-port S-parameter file, single frequency point", "freq magS11 angS11"], [2.0, 0.894, -12.136], Float64[]))
        )

        # Unknown format v1
        message = "Unknown format"
        io = IOBuffer();
        @test_throws(
            ErrorException(message),
            JosephsonCircuits.touchstone_write(io,JosephsonCircuits.TouchstoneFile([2.0e6], [0.874020294860635 - 0.18794819544685323im;;;], "mhz", "s", "da", 50.0, 1.0, 1, "12_21", 1, 0, [50.0], String[], "Full", Tuple{Char, Vector{Int64}}[], ["Example 1:", "1-port S-parameter file, single frequency point", "freq magS11 angS11"], [2.0, 0.894, -12.136], Float64[]))
        )

        # Unknown format v2
        message = "Unknown format"
        io = IOBuffer();
        @test_throws(
            ErrorException(message),
            JosephsonCircuits.touchstone_write(io,JosephsonCircuits.TouchstoneFile([2.0e6], [0.874020294860635 - 0.18794819544685323im;;;], "mhz", "s", "da", 50.0, 2.0, 1, "12_21", 1, 0, [50.0], String[], "Full", Tuple{Char, Vector{Int64}}[], ["Example 1:", "1-port S-parameter file, single frequency point", "freq magS11 angS11"], [2.0, 0.894, -12.136], Float64[]))
        )

    end

    @testset "touchstone_save RI" begin

        #find the temporary directory
        path = tempdir()

        #generate unique filenames
        filename = joinpath(path,"JosephsonCircuits-"* string(UUIDs.uuid1()) * ".ts")

        frequencies = [1.0e9, 2.0e9, 10.0e9]
        N = [0.3926 - 0.1211im -0.0003 - 0.0021im; -0.0003 - 0.0021im 0.3926 - 0.1211im;;; 0.3517 - 0.3054im -0.0096 - 0.0298im; -0.0096 - 0.0298im 0.3517 - 0.3054im;;; 0.3419 + 0.3336im -0.0134 + 0.0379im; -0.0134 + 0.0379im 0.3419 + 0.3336im]
        version = 2.0
        R = 50.0
        format = "RI"
        frequencyunit = "Hz"
        comments = ["Example 4:","2-port S-parameter file, three frequency points"]

        JosephsonCircuits.touchstone_save(filename, frequencies, N;
            version = version, R = R,format = format, frequencyunit = frequencyunit,
            comments = comments)

        out=JosephsonCircuits.touchstone_load(filename)

        # clean up the temporary file
        rm(filename)

        @test out.f == frequencies
        @test isapprox(out.N, N)
        @test out.version == version
        @test out.R == R
        @test out.format == lowercase(format)
        @test out.frequencyunit == lowercase(frequencyunit)

    end

    @testset "touchstone_save MA" begin

        #find the temporary directory
        path = tempdir()

        #generate unique filenames
        filename = joinpath(path,"JosephsonCircuits-"* string(UUIDs.uuid1()) * ".ts")

        frequencies = [1.0e9, 2.0e9, 10.0e9]
        N = [0.3926 - 0.1211im -0.0003 - 0.0021im; -0.0003 - 0.0021im 0.3926 - 0.1211im;;; 0.3517 - 0.3054im -0.0096 - 0.0298im; -0.0096 - 0.0298im 0.3517 - 0.3054im;;; 0.3419 + 0.3336im -0.0134 + 0.0379im; -0.0134 + 0.0379im 0.3419 + 0.3336im]
        version = 2.0
        R = 50.0
        format = "MA"
        frequencyunit = "Hz"
        comments = ["Example 4:","2-port S-parameter file, three frequency points"]

        JosephsonCircuits.touchstone_save(filename, frequencies, N;
            version = version, R = R,format = format, frequencyunit = frequencyunit,
            comments = comments)

        out=JosephsonCircuits.touchstone_load(filename)

        # clean up the temporary file
        rm(filename)

        @test out.f == frequencies
        @test isapprox(out.N, N)
        @test out.version == version
        @test out.R == R
        @test out.format == lowercase(format)
        @test out.frequencyunit == lowercase(frequencyunit)

    end

    @testset "touchstone_save dB" begin

        #find the temporary directory
        path = tempdir()

        #generate unique filenames
        filename = joinpath(path,"JosephsonCircuits-"* string(UUIDs.uuid1()) * ".ts")

        frequencies = [1.0e9, 2.0e9, 10.0e9]
        N = [0.3926 - 0.1211im -0.0003 - 0.0021im; -0.0003 - 0.0021im 0.3926 - 0.1211im;;; 0.3517 - 0.3054im -0.0096 - 0.0298im; -0.0096 - 0.0298im 0.3517 - 0.3054im;;; 0.3419 + 0.3336im -0.0134 + 0.0379im; -0.0134 + 0.0379im 0.3419 + 0.3336im]
        version = 2.0
        R = 50.0
        format = "dB"
        frequencyunit = "Hz"
        comments = ["Example 4:","2-port S-parameter file, three frequency points"]

        JosephsonCircuits.touchstone_save(filename, frequencies, N;
            version = version, R = R,format = format, frequencyunit = frequencyunit,
            comments = comments)

        out=JosephsonCircuits.touchstone_load(filename)

        # clean up the temporary file
        rm(filename)

        @test out.f == frequencies
        @test isapprox(out.N, N)
        @test out.version == version
        @test out.R == R
        @test out.format == lowercase(format)
        @test out.frequencyunit == lowercase(frequencyunit)

    end

    @testset "touchstone_save errors" begin

        #find the temporary directory
        path = tempdir()

        #generate unique filenames
        filename = joinpath(path,"JosephsonCircuits-"* string(UUIDs.uuid1()) * ".ts")

        frequencies = [1.0e9, 2.0e9, 10.0e9]
        N = [0.3926 - 0.1211im -0.0003 - 0.0021im; -0.0003 - 0.0021im 0.3926 - 0.1211im;;; 0.3517 - 0.3054im -0.0096 - 0.0298im; -0.0096 - 0.0298im 0.3517 - 0.3054im;;; 0.3419 + 0.3336im -0.0134 + 0.0379im; -0.0134 + 0.0379im 0.3419 + 0.3336im]
        version = 2.0
        R = 50.0
        format = "RI"
        frequencyunit = "Hz"
        comments = ["Example 4:","2-port S-parameter file, three frequency points"]

        # error if network data arrays not square
        message = "Network data arrays are not square."
        @test_throws(
            ErrorException(message),
            JosephsonCircuits.touchstone_save(filename, frequencies,
                N[1:end-1,:,:]; version = version, R = R,format = format,
                frequencyunit = frequencyunit,comments = comments)
        )

        # unknown format error
        message = "Unknown parameter"
        @test_throws(
            ErrorException(message),
            JosephsonCircuits.touchstone_save(filename, frequencies, N;
                version = version, R = R,format = format, parameter = "K",
                frequencyunit = frequencyunit, comments = comments)
        )

        # unknown format error
        message = "Unknown format"
        @test_throws(
            ErrorException(message),
            JosephsonCircuits.touchstone_save(filename, frequencies, N;
                version = version, R = R,format = "Unknown",
                frequencyunit = frequencyunit, comments = comments)
        )
        if isfile(filename)
            rm(filename)
        end

        #file extensions
        filename = joinpath(path,"JosephsonCircuits-"* string(UUIDs.uuid1()) * ".s4p")
        message = "Extension of .s4p is not the recommended extension of .ts or .s2p for a file with 2 ports."
        @test_warn message JosephsonCircuits.touchstone_save(filename, frequencies, N;
            version = version, R = R,format = format, frequencyunit = frequencyunit,
            comments = comments)
        rm(filename)

        # incorrect file extension v1 file
        filename = joinpath(path,"JosephsonCircuits-"* string(UUIDs.uuid1()) * ".s4p")
        message = "Extension of .s4p is not the recommended extension of .s2p for a version 1.0 file with 2 ports."
        @test_warn message JosephsonCircuits.touchstone_save(filename, frequencies, N;
            version = 1.0, R = R,format = format, frequencyunit = frequencyunit,
            comments = comments)
        rm(filename)

        # no file extension
        filename = joinpath(path,"JosephsonCircuits-"* string(UUIDs.uuid1()))
        message = "Adding extension of .s2p"
        @test_warn message JosephsonCircuits.touchstone_save(filename, frequencies, N;
            version = version, R = R,format = format, frequencyunit = frequencyunit,
            comments = comments)
        rm(filename * ".s2p")

    end

    @testset "matrixindices" begin
        @test_throws(
            ErrorException("Unknown two port data order string."),
            JosephsonCircuits.matrixindices(2,"Full","21_31")
        )
        @test_throws(
            ErrorException("Unknown matrix format."),
            JosephsonCircuits.matrixindices(2,"Fake","21_12")
        )
        @test_throws(
            ErrorException("Two port data order = 21_12 is only allowed if the number of ports is two."),
            JosephsonCircuits.matrixindices(4,"Full","21_12")
        )
    end

    @testset "frequencyscale" begin
        @test_throws(
            ErrorException("Unknown frequency unit THz"),
            JosephsonCircuits.frequencyscale("THz")
        )
    end

    @testset "parsetwoportdataorder" begin
        @test_throws(
            ErrorException("Unknown [Two-Port Data Order] parameter: [two-port data order] 21_34"),
            JosephsonCircuits.parsetwoportdataorder("[two-port data order] 21_34")
        )
    end

    @testset "parsenumberoffrequencies" begin
        @test_throws(
            ErrorException("Number of frequencies must be an integer greater than zero: [number of frequencies] 0"),
            JosephsonCircuits.parsenumberoffrequencies("[number of frequencies] 0")
        )
    end

    @testset "parsenumberofnoisefrequencies" begin
        @test_throws(
            ErrorException("Number of noise frequencies must be an integer greater than zero: [number of noise frequencies] 0"),
            JosephsonCircuits.parsenumberofnoisefrequencies("[number of noise frequencies] 0")
        )
    end

    @testset "parsereference!" begin
        io = IOBuffer("[Reference] 50.0 \n60.0 75.0 1.0\n[Number of Frequencies] 1")
        numberofports = 3
        comments = String[]
        reference = Float64[]
        line = JosephsonCircuits.stripcommentslowercase!(comments,readline(io))
        @test_throws(
            ErrorException("Too many values on [Reference] line: 60.0 75.0 1.0"),
            JosephsonCircuits.parsereference!(reference, comments, line, numberofports, io)
        )
    end

    @testset "parsematrixformat" begin
        @test_throws(
            ErrorException("Unknown format: unknown"),
            JosephsonCircuits.parsematrixformat("[matrix format] unknown")
        )
    end

    @testset "parsenetworkdata!" begin
        networkdata = Float64[]
        comments = String[]
        io = IOBuffer("2 .95 -26 3.57 157 .04 76 .66 -14\n# GHz S MA R 50\n22 .60 -144 1.30 40 .14 40 .56 -85\n[Noise Data]\n4 .7 .64 69 19\n18 2.7 .46 -33 20\n[End]")
        @test_throws(
            ErrorException("Second option line in network data."),
            JosephsonCircuits.parsenetworkdata!(networkdata,comments,io)
        )
    end

    @testset "parsenoisedata!"  begin

    networkdata = Float64[]
    noisedata = Float64[]
    comments = String[]
    io = IOBuffer("2 .95 -26 3.57 157 .04 76 .66 -14\n22 .60 -144 1.30 40 .14 40 .56 -85\n[Noise Data]\n4 .7 .64 69 19 1.2\n18 2.7 .46 -33 20\n[End]")
    JosephsonCircuits.parsenetworkdata!(networkdata,comments,io)
    @test_throws(
        ErrorException("Noise data lines must have 5 entries."),
        JosephsonCircuits.parsenoisedata!(noisedata,comments,io)
    )

    networkdata = Float64[]
    noisedata = Float64[]
    comments = String[]
    io = IOBuffer("2 .95 -26 3.57 157 .04 76 .66 -14\n22 .60 -144 1.30 40 .14 40 .56 -85\n[Noise Data]\n4 .7 .64 69 19\n[Noise Data]\n18 2.7 .46 -33 20\n[End]")
    JosephsonCircuits.parsenetworkdata!(networkdata,comments,io)
    @test_throws(
        ErrorException("Only one [Noise Data] keyword allowed."),
        JosephsonCircuits.parsenoisedata!(noisedata,comments,io)
    )

    networkdata = Float64[]
    noisedata = Float64[]
    comments = String[]
    io = IOBuffer("2 .95 -26 3.57 157 .04 76 .66 -14\n22 .60 -144 1.30 40 .14 40 .56 -85\n[Noise Data]\n4 .7 .64 69 19\n# GHz S MA R 50\n18 2.7 .46 -33 20\n[End]")
    JosephsonCircuits.parsenetworkdata!(networkdata,comments,io)
    @test_throws(
        ErrorException("Second option line in noise data."),
        JosephsonCircuits.parsenoisedata!(noisedata,comments,io)
    )

    networkdata = Float64[]
    noisedata = Float64[]
    comments = String[]
    io = IOBuffer("2 .95 -26 3.57 157 .04 76 .66 -14\n22 .60 -144 1.30 40 .14 40 .56 -85\n[Noise Data]\n18 2.7 .46 -33 20\n4 .7 .64 69 19\n[End]")
    JosephsonCircuits.parsenetworkdata!(networkdata,comments,io)
    @test_throws(
        ErrorException("Frequencies descending in noise data"),
        JosephsonCircuits.parsenoisedata!(noisedata,comments,io)
    )

    end

end