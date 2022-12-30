using JosephsonCircuits
using Test
import UUIDs

@testset verbose=true "touchstone" begin

    @testset "touchstone_parse" begin

        examples = [
            "!Example 1:\n!1-port S-parameter file, single frequency point\n# MHz S MA R 50\n!freq magS11 angS11\n2.000 0.894 -12.136",
            "!Example 2:\n!1-port Z-parameter file, multiple frequency points\n# MHz Z MA R 75\n!freq magZ11 angZ11\n100 0.99 -4\n200 0.80 -22\n300 0.707 -45\n400 0.40 -62\n500 0.01 -89",
            "!Example 3:\n!2-port H-parameter file, single frequency point\n# KHz H MA R 1\n! freq magH11 angH11 magH21 angH21 magH12 angH12 magH22 angH22\n2 .95 -26 3.57 157 .04 76 .66 -1",
            "!Example 4:\n!2-port S-parameter file, three frequency points\n# GHZ S RI R 50.0\n!freq RelS11 ImS11 ReS21 ImS21 ReS12 ImS12 ReS22 ImS22\n1.0000 0.3926 -0.1211 -0.0003 -0.0021 -0.0003 -0.0021 0.3926 -0.1211\n2.0000 0.3517 -0.3054 -0.0096 -0.0298 -0.0096 -0.0298 0.3517 -0.3054\n10.000 0.3419 0.3336 -0.0134 0.0379 -0.0134 0.0379 0.3419 0.3336",
            "!Example 5:\n! 4-port S-parameter data, taken at three frequency points\n# GHZ S MA R 50\n5.00000 0.60 161.24 0.40 -42.20 0.42 -66.58 0.53 -79.34 !row 1\n0.40 -42.20 0.60 161.20 0.53 -79.34 0.42 -66.58 !row 2\n0.42 -66.58 0.53 -79.34 0.60 161.24 0.40 -42.20 !row 3\n0.53 -79.34 0.42 -66.58 0.40 -42.20 0.60 161.24 !row 4\n6.00000 0.57 150.37 0.40 -44.34 0.41 -81.24 0.57 -95.77 !row 1\n0.40 -44.34 0.57 150.37 0.57 -95.77 0.41 -81.24 !row 2\n0.41 -81.24 0.57 -95.77 0.57 150.37 0.40 -44.34 !row 3\n0.57 -95.77 0.41 -81.24 0.40 -44.34 0.57 150.37 !row 4\n7.00000 0.50 136.69 0.45 -46.41 0.37 -99.09 0.62 -114.19 !row 1\n0.45 -46.41 0.50 136.69 0.62 -114.19 0.37 -99.09 !row 2\n0.37 -99.09 0.62 -114.19 0.50 136.69 0.45 -46.41 !row 3\n0.62 -114.19 0.37 -99.09 0.45 -46.41 0.50 136.69 !row 4",
            "!Example 8:\n!2-port network, S-parameter and noise data\n# GHZ S MA R 50\n2 .95 -26 3.57 157 .04 76 .66 -14\n22 .60 -144 1.30 40 .14 40 .56 -85\n! NOISE PARAMETERS\n4 .7 .64 69 .38\n18 2.7 .46 -33 .40",
            "!Example 4 (Version 2.0):\n! 4-port S-parameter data\n! Default impedance is overridden by [Reference]\n! Data cannot be represented using 1.0 syntax\n! Note that the [Reference] keyword arguments appear on a separate line\n[Version] 2.0\n# GHz S MA R 50\n[Number of Ports] 4\n[Reference]\n50 75 0.01 0.01\n[Number of Frequencies] 1\n[Network Data]\n5.00000 0.60 161.24 0.40 -42.20 0.42 -66.58 0.53 -79.34 !row 1\n0.40 -42.20 0.60 161.20 0.53 -79.34 0.42 -66.58 !row 2\n0.42 -66.58 0.53 -79.34 0.60 161.24 0.40 -42.20 !row 3\n0.53 -79.34 0.42 -66.58 0.40 -42.20 0.60 161.24 !row 4",
            "!Example 5 (Version 2.0):\n! 4-port S-parameter data\n! Default impedance is overridden by the [Reference] keyword arguments\n! Data cannot be represented using 1.0 syntax\n[Version] 2.0\n# GHz S MA R 50\n[Number of Ports] 4\n[Number of Frequencies] 1\n[Reference] 50 75 0.01 0.01\n[Matrix Format] Full\n[Network Data]\n5.00000 0.60 161.24 0.40 -42.20 0.42 -66.58 0.53 -79.34 !row 1\n0.40 -42.20 0.60 161.20 0.53 -79.34 0.42 -66.58 !row 2\n0.42 -66.58 0.53 -79.34 0.60 161.24 0.40 -42.20 !row 3\n0.53 -79.34 0.42 -66.58 0.40 -42.20 0.60 161.24 !row 4\n[End]",
            # "!Example 6 (Version 2.0):\n! 4-port S-parameter data\n! Default impedance is overridden by the [Reference] keyword arguments\n! Note that [Reference] arguments are split across two lines\n! Data cannot be represented using 1.0 syntax\n[Version] 2.0\n# GHz S MA R 50\n[Number of Ports] 4\n[Number of Frequencies] 1\n[Reference] 50 75\n0.01 0.01\n[Matrix Format] Lower\n[Network Data]\n5.00000 0.60 161.24 !row 1\n0.40 -42.20 0.60 161.20 !row 2\n0.42 -66.58 0.53 -79.34 0.60 161.24 !row 3\n0.53 -79.34 0.42 -66.58 0.40 -42.20 0.60 161.24 !row 4\n[End]",
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
            JosephsonCircuits.TouchstoneFile([1.0e8, 2.0e8, 3.0e8, 4.0e8, 5.0e8], [74.06913073179194 - 5.179418175501303im;;; 55.63103127400724 - 22.47639560495472im;;; 37.494337072416684 - 37.49433707241668im;;; 14.084146883576725 - 26.488427785767808im;;; 0.013089304827962698 - 0.7498857713672935im], "mhz", "z", "ma", 75.0, 1.0, 1, "12_21", 5, 0, [75.0], String[], "Full", Tuple{Char, Vector{Int64}}[], ["Example 2:", "1-port Z-parameter file, multiple frequency points", "freq magZ11 angZ11"], [100.0, 0.99, -4.0, 200.0, 0.8, -22.0, 300.0, 0.707, -45.0, 400.0, 0.4, -62.0, 500.0, 0.01, -89.0], Float64[]),
            JosephsonCircuits.TouchstoneFile([2000.0], [0.8538543439842087 - 0.4164525894496235im 0.009676875823986707 + 0.03881182905103986im; -3.286202326825212 + 1.3949101287067074im 0.6598994788032183 - 0.011518588248607119im;;;], "khz", "h", "ma", 1.0, 1.0, 2, "21_12", 1, 0, [1.0, 1.0], String[], "Full", Tuple{Char, Vector{Int64}}[], ["Example 3:", "2-port H-parameter file, single frequency point", " freq magH11 angH11 magH21 angH21 magH12 angH12 magH22 angH22"], [2.0, 0.95, -26.0, 3.57, 157.0, 0.04, 76.0, 0.66, -1.0], Float64[]),
            JosephsonCircuits.TouchstoneFile([1.0e9, 2.0e9, 1.0e10], [0.3926 - 0.1211im -0.0003 - 0.0021im; -0.0003 - 0.0021im 0.3926 - 0.1211im;;; 0.3517 - 0.3054im -0.0096 - 0.0298im; -0.0096 - 0.0298im 0.3517 - 0.3054im;;; 0.3419 + 0.3336im -0.0134 + 0.0379im; -0.0134 + 0.0379im 0.3419 + 0.3336im], "ghz", "s", "ri", 50.0, 1.0, 2, "21_12", 3, 0, [50.0, 50.0], String[], "Full", Tuple{Char, Vector{Int64}}[], ["Example 4:", "2-port S-parameter file, three frequency points", "freq RelS11 ImS11 ReS21 ImS21 ReS12 ImS12 ReS22 ImS22"], [1.0, 0.3926, -0.1211, -0.0003, -0.0021, -0.0003, -0.0021, 0.3926, -0.1211, 2.0, 0.3517, -0.3054, -0.0096, -0.0298, -0.0096, -0.0298, 0.3517, -0.3054, 10.0, 0.3419, 0.3336, -0.0134, 0.0379, -0.0134, 0.0379, 0.3419, 0.3336], Float64[]),
            JosephsonCircuits.TouchstoneFile([5.0e9, 6.0e9, 7.0e9], [-0.5681244079815996 + 0.1929628385351877im 0.29632183851470006 - 0.2686882357291961im 0.16693665375723588 - 0.38539869438327984im 0.09803970583787712 - 0.5208533537179372im; 0.29632183851470006 - 0.2686882357291961im -0.5679895560694177 + 0.1933594171383067im 0.09803970583787712 - 0.5208533537179372im 0.16693665375723588 - 0.38539869438327984im; 0.16693665375723588 - 0.38539869438327984im 0.09803970583787712 - 0.5208533537179372im -0.5681244079815996 + 0.1929628385351877im 0.29632183851470006 - 0.2686882357291961im; 0.09803970583787712 - 0.5208533537179372im 0.16693665375723588 - 0.38539869438327984im 0.29632183851470006 - 0.2686882357291961im -0.5681244079815996 + 0.1929628385351877im;;; -0.495464624294036 + 0.28180632724119187im 0.286081989392916 - 0.2795659051905141im 0.062441313054034775 - 0.40521732739862937im -0.05730515806890161 - 0.5671120866801361im; 0.286081989392916 - 0.2795659051905141im -0.495464624294036 + 0.28180632724119187im -0.05730515806890161 - 0.5671120866801361im 0.062441313054034775 - 0.40521732739862937im; 0.062441313054034775 - 0.40521732739862937im -0.05730515806890161 - 0.5671120866801361im -0.495464624294036 + 0.28180632724119187im 0.286081989392916 - 0.2795659051905141im; -0.05730515806890161 - 0.5671120866801361im 0.062441313054034775 - 0.40521732739862937im 0.286081989392916 - 0.2795659051905141im -0.495464624294036 + 0.28180632724119187im;;; -0.3638265243449566 + 0.3429726813946975im 0.3102719136297667 - 0.32593149527549903im -0.05845471959176759 - 0.3653533163356367im -0.2540535762162701 - 0.565558821354352im; 0.3102719136297667 - 0.32593149527549903im -0.3638265243449566 + 0.3429726813946975im -0.2540535762162701 - 0.565558821354352im -0.05845471959176759 - 0.3653533163356367im; -0.05845471959176759 - 0.3653533163356367im -0.2540535762162701 - 0.565558821354352im -0.3638265243449566 + 0.3429726813946975im 0.3102719136297667 - 0.32593149527549903im; -0.2540535762162701 - 0.565558821354352im -0.05845471959176759 - 0.3653533163356367im 0.3102719136297667 - 0.32593149527549903im -0.3638265243449566 + 0.3429726813946975im], "ghz", "s", "ma", 50.0, 1.0, 4, "12_21", 3, 0, [50.0, 50.0, 50.0, 50.0], String[], "Full", Tuple{Char, Vector{Int64}}[], ["Example 5:", " 4-port S-parameter data, taken at three frequency points", "row 1", "row 2", "row 3", "row 4", "row 1", "row 2", "row 3", "row 4", "row 1", "row 2", "row 3", "row 4"], [5.0, 0.6, 161.24, 0.4, -42.2, 0.42, -66.58, 0.53, -79.34, 0.4, -42.2, 0.6, 161.2, 0.53, -79.34, 0.42, -66.58, 0.42, -66.58, 0.53, -79.34, 0.6, 161.24, 0.4, -42.2, 0.53, -79.34, 0.42, -66.58, 0.4, -42.2, 0.6, 161.24, 6.0, 0.57, 150.37, 0.4, -44.34, 0.41, -81.24, 0.57, -95.77, 0.4, -44.34, 0.57, 150.37, 0.57, -95.77, 0.41, -81.24, 0.41, -81.24, 0.57, -95.77, 0.57, 150.37, 0.4, -44.34, 0.57, -95.77, 0.41, -81.24, 0.4, -44.34, 0.57, 150.37, 7.0, 0.5, 136.69, 0.45, -46.41, 0.37, -99.09, 0.62, -114.19, 0.45, -46.41, 0.5, 136.69, 0.62, -114.19, 0.37, -99.09, 0.37, -99.09, 0.62, -114.19, 0.5, 136.69, 0.45, -46.41, 0.62, -114.19, 0.37, -99.09, 0.45, -46.41, 0.5, 136.69], Float64[]),
            JosephsonCircuits.TouchstoneFile([2.0e9, 2.2e10], [0.8538543439842087 - 0.4164525894496235im 0.009676875823986707 + 0.03881182905103986im; -3.286202326825212 + 1.3949101287067074im 0.6403951793421577 - 0.1596684510957807im;;; -0.48541019662496837 - 0.35267115137548394im 0.10724622203665693 + 0.0899902653561155im; 0.9958577760546714 + 0.835623892592501im 0.048807215938688565 - 0.5578690309313775im], "ghz", "s", "ma", 50.0, 1.0, 2, "21_12", 2, 2, [50.0, 50.0], String[], "Full", Tuple{Char, Vector{Int64}}[], ["Example 8:", "2-port network, S-parameter and noise data", " NOISE PARAMETERS"], [2.0, 0.95, -26.0, 3.57, 157.0, 0.04, 76.0, 0.66, -14.0, 22.0, 0.6, -144.0, 1.3, 40.0, 0.14, 40.0, 0.56, -85.0], [4.0, 0.7, 0.64, 69.0, 0.38, 18.0, 2.7, 0.46, -33.0, 0.4]),
            JosephsonCircuits.TouchstoneFile([5.0e9], [-0.5681244079815996 + 0.1929628385351877im 0.29632183851470006 - 0.2686882357291961im 0.16693665375723588 - 0.38539869438327984im 0.09803970583787712 - 0.5208533537179372im; 0.29632183851470006 - 0.2686882357291961im -0.5679895560694177 + 0.1933594171383067im 0.09803970583787712 - 0.5208533537179372im 0.16693665375723588 - 0.38539869438327984im; 0.16693665375723588 - 0.38539869438327984im 0.09803970583787712 - 0.5208533537179372im -0.5681244079815996 + 0.1929628385351877im 0.29632183851470006 - 0.2686882357291961im; 0.09803970583787712 - 0.5208533537179372im 0.16693665375723588 - 0.38539869438327984im 0.29632183851470006 - 0.2686882357291961im -0.5681244079815996 + 0.1929628385351877im;;;], "ghz", "s", "ma", 50.0, 2.0, 4, "12_21", 1, 0, [50.0, 75.0, 0.01, 0.01], String[], "Full", Tuple{Char, Vector{Int64}}[], ["Example 4 (Version 2.0):", " 4-port S-parameter data", " Default impedance is overridden by [Reference]", " Data cannot be represented using 1.0 syntax", " Note that the [Reference] keyword arguments appear on a separate line", "row 1", "row 2", "row 3", "row 4"], [5.0, 0.6, 161.24, 0.4, -42.2, 0.42, -66.58, 0.53, -79.34, 0.4, -42.2, 0.6, 161.2, 0.53, -79.34, 0.42, -66.58, 0.42, -66.58, 0.53, -79.34, 0.6, 161.24, 0.4, -42.2, 0.53, -79.34, 0.42, -66.58, 0.4, -42.2, 0.6, 161.24], Float64[]),
            JosephsonCircuits.TouchstoneFile([5.0e9], [-0.5681244079815996 + 0.1929628385351877im 0.29632183851470006 - 0.2686882357291961im 0.16693665375723588 - 0.38539869438327984im 0.09803970583787712 - 0.5208533537179372im; 0.29632183851470006 - 0.2686882357291961im -0.5679895560694177 + 0.1933594171383067im 0.09803970583787712 - 0.5208533537179372im 0.16693665375723588 - 0.38539869438327984im; 0.16693665375723588 - 0.38539869438327984im 0.09803970583787712 - 0.5208533537179372im -0.5681244079815996 + 0.1929628385351877im 0.29632183851470006 - 0.2686882357291961im; 0.09803970583787712 - 0.5208533537179372im 0.16693665375723588 - 0.38539869438327984im 0.29632183851470006 - 0.2686882357291961im -0.5681244079815996 + 0.1929628385351877im;;;], "ghz", "s", "ma", 50.0, 2.0, 4, "12_21", 1, 0, [50.0, 75.0, 0.01, 0.01], String[], "Full", Tuple{Char, Vector{Int64}}[], ["Example 5 (Version 2.0):", " 4-port S-parameter data", " Default impedance is overridden by the [Reference] keyword arguments", " Data cannot be represented using 1.0 syntax", "row 1", "row 2", "row 3", "row 4"], [5.0, 0.6, 161.24, 0.4, -42.2, 0.42, -66.58, 0.53, -79.34, 0.4, -42.2, 0.6, 161.2, 0.53, -79.34, 0.42, -66.58, 0.42, -66.58, 0.53, -79.34, 0.6, 161.24, 0.4, -42.2, 0.53, -79.34, 0.42, -66.58, 0.4, -42.2, 0.6, 161.24], Float64[]),
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

    @testset "touchstone_save" begin

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
        @test out.N == N
        @test out.version == version
        @test out.R == R
        @test out.format == lowercase(format)
        @test out.frequencyunit == lowercase(frequencyunit)

    end

end