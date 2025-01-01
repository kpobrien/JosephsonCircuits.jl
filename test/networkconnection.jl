using JosephsonCircuits
using LinearAlgebra
using Test
import StaticArrays

@testset verbose=true "network connection" begin

    # one network
    @testset "connectS one network errors" begin
        begin
            Sa = Complex{Float64}[1 2;3 4]
            @test_throws(
                ArgumentError("Port `k` is smaller than one."),
                JosephsonCircuits.connectS(Sa,0,1)
            )

            @test_throws(
                ArgumentError("Port `l` is smaller than one."),
                JosephsonCircuits.connectS(Sa,1,0)
            )

            @test_throws(
                ArgumentError("Port `k` is larger than number of ports in `Sa`."),
                JosephsonCircuits.connectS(Sa,3,1)
            )

            @test_throws(
                ArgumentError("Port `l` is larger than number of ports in `Sa`."),
                JosephsonCircuits.connectS(Sa,1,3)
            )

            @test_throws(
                ArgumentError("`k` and `l` cannot be equal because a port cannot be merged with itself."),
                JosephsonCircuits.connectS(Sa,1,1)
            )
        end

        begin
            Sa = Complex{Float64}[1 2 3;4 5 6]
            @test_throws(
                DimensionMismatch("Lengths of first two dimensions of `Sa` must be equal."),
                JosephsonCircuits.connectS(Sa,1,2)
            )
        end

        # one network, in place
        begin
            Sa = rand(Complex{Float64},3,3)
            Sout = zeros(Complex{Float64},2,2)
            @test_throws(
                DimensionMismatch("Length of first two dimensions must be 2 smaller for `Sout` than `Sa` because we are merging two ports."),
                JosephsonCircuits.connectS!(Sout,Sa,1,2)
            )
        end
    end

    @testset "connectS! one network errors" begin
        begin
            Sa = rand(Complex{Float64},4,4)
            Sout = zeros(Complex{Float64},2,2,1)
            @test_throws(
                DimensionMismatch("`Sout` and `Sa` must have the same number of dimensions."),
                JosephsonCircuits.connectS!(Sout,Sa,1,2)
            )
        end

        begin
            Sa = rand(Complex{Float64},4)
            Sout = zeros(Complex{Float64},2)
            @test_throws(
                DimensionMismatch("`Sout` and `Sa` must have at least two dimensions."),
                JosephsonCircuits.connectS!(Sout,Sa,1,2)
            )
        end

        begin
            Sa = rand(Complex{Float64},4,4)
            Sout = zeros(Complex{Float64},2,3)
            @test_throws(
                DimensionMismatch("Lengths of first two dimensions of `Sout` must be equal."),
                JosephsonCircuits.connectS!(Sout,Sa,1,2)
            )
        end

        begin
            Sa = rand(Complex{Float64},3,3,3)
            Sout = zeros(Complex{Float64},1,1,4)
            @test_throws(
                DimensionMismatch("Non-port axis lengths of `Sa` and `Sout` must be equal."),
                JosephsonCircuits.connectS!(Sout,Sa,1,2)
            )
        end
    end

    # two networks
    @testset "connectS two networks errors" begin
        begin
            Sa = Complex{Float64}[1 2;3 4]
            Sb = Complex{Float64}[5 6;7 8]
            @test_throws(
                ArgumentError("Port `k` is smaller than one."),
                JosephsonCircuits.connectS(Sa,Sb,0,1)
            )

            @test_throws(
                ArgumentError("Port `l` is smaller than one."),
                JosephsonCircuits.connectS(Sa,Sb,1,0)
            )

            @test_throws(
                ArgumentError("Port `k` is larger than number of ports in `Sa`."),
                JosephsonCircuits.connectS(Sa,Sb,3,1)
            )

            @test_throws(
                ArgumentError("Port `l` is larger than number of ports in `Sb`."),
                JosephsonCircuits.connectS(Sa,Sb,1,3)
            )
        end

        begin
            Sa = Complex{Float64}[1 2 3;4 5 6]
            Sb = Complex{Float64}[5 6;7 8]
            @test_throws(
                DimensionMismatch("Lengths of first two dimensions of `Sa` must be equal."),
                JosephsonCircuits.connectS(Sa,Sb,1,2)
            )
        end

        begin
            Sa = Complex{Float64}[1 2;4 5]
            Sb = Complex{Float64}[5 6 7;8 9 10]
            @test_throws(
                DimensionMismatch("Lengths of first two dimensions of `Sb` must be equal."),
                JosephsonCircuits.connectS(Sa,Sb,1,2)
            )
        end
    end

    # two networks, in place
    @testset "connectS! two networks errors" begin
        begin
            Sa = rand(Complex{Float64},3,3)
            Sb = rand(Complex{Float64},3,3)
            Sout = zeros(Complex{Float64},2,2)
            @test_throws(
                DimensionMismatch("First two dimensions of `Sout` must be `m+n-2`."),
                JosephsonCircuits.connectS!(Sout,Sa,Sb,1,2)
            )
        end

        begin
            Sa = rand(Complex{Float64},3,3)
            Sb = rand(Complex{Float64},3,3,1)
            Sout = zeros(Complex{Float64},4,4)
            @test_throws(
                DimensionMismatch("`Sa` and `Sb` must have the same number of dimensions."),
                JosephsonCircuits.connectS!(Sout,Sa,Sb,1,2)
            )
        end

        begin
            Sa = rand(Complex{Float64},4,4)
            Sb = rand(Complex{Float64},4,4)
            Sout = zeros(Complex{Float64},6,6,1)
            @test_throws(
                DimensionMismatch("`Sout`, `Sa`, and `Sb` must have the same number of dimensions."),
                JosephsonCircuits.connectS!(Sout,Sa,Sb,1,2)
            )
        end

        begin
            Sa = rand(Complex{Float64},4)
            Sb = rand(Complex{Float64},4)
            Sout = zeros(Complex{Float64},6)
            @test_throws(
                DimensionMismatch("`Sout`, `Sa`, and `Sb` must have atleast two dimensions."),
                JosephsonCircuits.connectS!(Sout,Sa,Sb,1,2)
            )
        end

        begin
            Sa = rand(Complex{Float64},4,4)
            Sb = rand(Complex{Float64},4,4)
            Sout = zeros(Complex{Float64},6,7)
            @test_throws(
                DimensionMismatch("Lengths of first two dimensions of `Sout` must be equal."),
                JosephsonCircuits.connectS!(Sout,Sa,Sb,1,2)
            )
        end

        begin
            Sa = rand(Complex{Float64},3,3,3)
            Sb = rand(Complex{Float64},3,3,3)
            Sout = zeros(Complex{Float64},4,4,4)
            @test_throws(
                DimensionMismatch("Non-port axis lengths of `Sa`, `Sb`, and `Sout` must be equal."),
                JosephsonCircuits.connectS!(Sout,Sa,Sb,1,2)
            )
        end
    end

    # check for consistency between the one and two network functions
    @testset "connectS consistency" begin

        begin
            Sa = rand(Complex{Float64},3,3,3)
            Sb = rand(Complex{Float64},3,3,3)
            Sboth = zeros(Complex{Float64},6,6,3)
            Sboth[1:size(Sa,1),1:size(Sa,1),:,:] .= Sa
            Sboth[size(Sa,1)+1:end,size(Sa,1)+1:end,:,:] .= Sb
            Sout1 = JosephsonCircuits.connectS(Sboth,1,2+size(Sa,1))
            Sout2 = JosephsonCircuits.connectS(Sa,Sb,1,2)
            @test isapprox(Sout1,Sout2)
        end

        begin
            Sa = rand(Complex{Float64},3,3)
            Sb = rand(Complex{Float64},3,3)
            Sboth = zeros(Complex{Float64},6,6)
            Sboth[1:size(Sa,1),1:size(Sa,1),:,:] .= Sa
            Sboth[size(Sa,1)+1:end,size(Sa,1)+1:end,:,:] .= Sb
            Sout1 = JosephsonCircuits.connectS(Sboth,2,1+size(Sa,1))
            Sout2 = JosephsonCircuits.connectS(Sa,Sb,2,1)
            @test isapprox(Sout1,Sout2)
        end

        begin
            Sa = rand(Complex{Float64},3,3,3)
            Sb = rand(Complex{Float64},3,3,3)
            Sboth = zeros(Complex{Float64},6,6,3)
            Sboth[1:size(Sa,1),1:size(Sa,1),:,:] .= Sa
            Sboth[size(Sa,1)+1:end,size(Sa,1)+1:end,:,:] .= Sb
            Sboth = [StaticArrays.MMatrix{size(Sboth,1),size(Sboth,2)}(Sboth[:,:,i]) for i=1:size(Sboth,3)]
            Sa = [StaticArrays.MMatrix{size(Sa,1),size(Sa,2)}(Sa[:,:,i]) for i=1:size(Sa,3)]
            Sb = [StaticArrays.MMatrix{size(Sb,1),size(Sb,2)}(Sb[:,:,i]) for i=1:size(Sb,3)]
            Sout1 = JosephsonCircuits.connectS.(Sboth,1,2+size(Sa,1))
            Sout2 = JosephsonCircuits.connectS.(Sa,Sb,1,2)
            @test isapprox(Sout1,Sout2)
        end

        begin
            Sa = rand(Complex{Float64},3,3)
            Sb = rand(Complex{Float64},3,3)
            Sboth = zeros(Complex{Float64},6,6)
            Sboth[1:size(Sa,1),1:size(Sa,1),:,:] .= Sa
            Sboth[size(Sa,1)+1:end,size(Sa,1)+1:end,:,:] .= Sb
            # convert to MMatrix
            Sboth = StaticArrays.MMatrix{size(Sboth,1),size(Sboth,2)}(Sboth)
            Sa = StaticArrays.MMatrix{size(Sa,1),size(Sa,2)}(Sa)
            Sb = StaticArrays.MMatrix{size(Sb,1),size(Sb,2)}(Sb)
            Sout1 = JosephsonCircuits.connectS(Sboth,2,1+size(Sa,1))
            Sout2 = JosephsonCircuits.connectS(Sa,Sb,2,1)
            @test isapprox(Sout1,Sout2)
        end
    end

    # one network
    @testset "connectSports one network" begin
        portsa = [(:S1,1),(:S1,2)]
        @test_throws(
            ArgumentError("Port `k` is smaller than one."),
            JosephsonCircuits.connectSports(portsa,0,1)
        )

        @test_throws(
            ArgumentError("Port `l` is smaller than one."),
            JosephsonCircuits.connectSports(portsa,1,0)
        )

        @test_throws(
            ArgumentError("Port `k` is larger than number of ports."),
            JosephsonCircuits.connectSports(portsa,3,1)
        )

        @test_throws(
            ArgumentError("Port `l` is larger than number of ports."),
            JosephsonCircuits.connectSports(portsa,1,3)
        )

        @test_throws(
            ArgumentError("`k` and `l` cannot be equal because a port cannot be merged with itself."),
            JosephsonCircuits.connectSports(portsa,1,1)
        )
    end

    # two networks
    @testset "connectSports two networks" begin
        portsa = [(:S1,1),(:S1,2)]
        portsb = [(:S2,1),(:S2,2)]
        @test_throws(
            ArgumentError("Port `k` is smaller than one."),
            JosephsonCircuits.connectSports(portsa,portsb,0,1)
        )

        @test_throws(
            ArgumentError("Port `l` is smaller than one."),
            JosephsonCircuits.connectSports(portsa,portsb,1,0)
        )

        @test_throws(
            ArgumentError("Port `k` is larger than number of ports in `portsa`."),
            JosephsonCircuits.connectSports(portsa,portsb,3,1)
        )

        @test_throws(
            ArgumentError("Port `l` is larger than number of ports in `portsb`."),
            JosephsonCircuits.connectSports(portsa,portsb,1,3)
        )
    end

    @testset "connectS with list of connections" begin
        # define an open
        Sopen = ones(Complex{Float64},1,1)

        # and a short
        Sshort = -ones(Complex{Float64},1,1)

        # and a match
        Smatch = zeros(Complex{Float64},1,1)

        # a splitter
        Ssplitter = Complex{Float64}[-1/3 2/3 2/3;2/3 -1/3 2/3;2/3 2/3 -1/3]

        S1 = rand(Complex{Float64},3,3)
        S2 = rand(Complex{Float64},2,2)

        # with strings
        networks = [("S1",S1),("S2",S2),("S3",Ssplitter),("S4",Sopen)]
        connections = [("S1","S1",1,2),("S1","S2",3,1),("S3","S2",2,2),("S3","S4",3,1)]
        Sout1, ports = JosephsonCircuits.connectS(networks,connections)

        Sout2 = begin
            S = JosephsonCircuits.connectS(Ssplitter,Sopen,3,1)
            S = JosephsonCircuits.connectS(S,S2,2,2)
            S = JosephsonCircuits.connectS(S1,S,3,2)
            S = JosephsonCircuits.connectS(S,1,2)
        end
        @test isapprox(Sout1[1],Sout2)

    end

    @testset "connectS with list of connections, many frequencies" begin
        N = 100

        # define an open
        Sopen = ones(Complex{Float64},1,1,N)

        # and a short
        Sshort = -ones(Complex{Float64},1,1,N)

        # and a match
        Smatch = zeros(Complex{Float64},1,1,N)

        # a splitter
        Ssplitter = stack([Complex{Float64}[-1/3 2/3 2/3;2/3 -1/3 2/3;2/3 2/3 -1/3] for i in 1:N])
        
        S1 = rand(Complex{Float64},3,3,N)
        S2 = rand(Complex{Float64},2,2,N)

        networks = [("S1",S1),("S2",S2),("S3",Ssplitter),("S4",Sopen)]
        connections = [("S1","S1",1,2),("S1","S2",3,1),("S3","S2",2,2),("S3","S4",3,1)]
        networkdata, ports = JosephsonCircuits.connectS(networks,connections)
        Sout1 = networkdata[1]

        Sout2 = begin
            S = JosephsonCircuits.connectS(Ssplitter,Sopen,3,1)
            S = JosephsonCircuits.connectS(S,S2,2,2)
            S = JosephsonCircuits.connectS(S1,S,3,2)
            S = JosephsonCircuits.connectS(S,1,2)
        end

        @test isapprox(Sout1,Sout2)

        Sout3 = JosephsonCircuits.solveS(networks,connections)

        @test isapprox(Sout1,Sout3[1])

    end

    @testset "connectS with list of connections, small_splitters" begin
        N = 100

        port1 = rand(Complex{Float64},2,2,N)
        port2 = rand(Complex{Float64},2,2,N)
        port3 = rand(Complex{Float64},2,2,N)
        port4 = rand(Complex{Float64},2,2,N)
        port5 = rand(Complex{Float64},2,2,N)

        networks = [
            ("port1",port1),("port2",port2),("port3",port3),("port4",port4),
            ("port5",port5),
            ]

        connections = [
            [("port1", 1),("port2", 1),("port3", 1),("port4", 1),("port5", 1)],
        ]

        out1 = JosephsonCircuits.connectS(networks, connections;
            small_splitters=false);

        out2 = JosephsonCircuits.connectS(networks, connections;
            small_splitters=true);
        @test isapprox(out1[1],out2[1])

        # test solveS
        out3 = JosephsonCircuits.solveS(networks, connections;
            small_splitters=false,
            factorization=JosephsonCircuits.LUfactorization())
        @test isapprox(out1[1][1],out3[1])

        out4 = JosephsonCircuits.solveS(networks, connections;
            small_splitters=false,
            factorization=JosephsonCircuits.KLUfactorization())
        @test isapprox(out1[1][1],out4[1])

        out5 = JosephsonCircuits.solveS(networks, connections;
            small_splitters=true,
            factorization=JosephsonCircuits.LUfactorization())
        @test isapprox(out1[1][1],out5[1])

        out6 = JosephsonCircuits.solveS(networks, connections;
            small_splitters=true,
            factorization=JosephsonCircuits.KLUfactorization())
        @test isapprox(out1[1][1],out6[1])

    end

    @testset "make_connection! errors" begin
        begin
            networks = [("S1",Complex{Float64}[0.0 1.0;1.0 0.0]),("S2",Complex{Float64}[0.5 0.5;0.5 0.5])];
            connections = [("S1","S2",1,2)];
            userinput = ones(Bool,length(networks))
            storage = Dict{Int,typeof(networks[1][2])}()
            g, fconnectionlist, fweightlist, ports, networkdata = JosephsonCircuits.connectS_initialize(networks,connections)
            # corrupt the vector of ports
            ports = [[("S1", 1), ("S1", 2)],[("S2", 1), ("S3", 2)]]
            @test_throws(
                ArgumentError("""Destination port ("S2", 2) not found in the ports [("S2", 1), ("S3", 2)] of the destination node 2."""),
                JosephsonCircuits.make_connection!(g, fconnectionlist, fweightlist, ports, networkdata,1,1,1,userinput,storage)
            )
        end

        begin
            networks = [("S1",Complex{Float64}[0.0 1.0;1.0 0.0]),("S2",Complex{Float64}[0.5 0.5;0.5 0.5])];
            connections = [("S1","S2",1,2)];
            userinput = ones(Bool,length(networks))
            storage = Dict{Int,typeof(networks[1][2])}()
            g, fconnectionlist, fweightlist, ports, networkdata = JosephsonCircuits.connectS_initialize(networks,connections)
            # corrupt the vector of ports
            ports = [[("S3", 1), ("S1", 2)],[("S2", 1), ("S2", 2)]]
            @test_throws(
                ArgumentError("""Source port ("S1", 1) not found in the ports [("S3", 1), ("S1", 2)] of the source node 1."""),
                JosephsonCircuits.make_connection!(g, fconnectionlist, fweightlist, ports, networkdata,1,1,1,userinput,storage)
            )
        end

    end

    @testset "connectS_initialize errors" begin

        begin
            networks = [("S1",Complex{Float64}[0.0 1.0;1.0 0.0]),("S2",Complex{Float64}[0.5 0.5;0.5 0.5])];
            connections = [("S3","S2",1,2)];
            @test_throws(
                ArgumentError("Source (network name, port number) (S3, 1) not found for connection (S3,S2,1,2)."),
                JosephsonCircuits.connectS_initialize(networks,connections)
            )
        end

        begin
            networks = [("S1",Complex{Float64}[0.0 1.0;1.0 0.0]),("S2",Complex{Float64}[0.5 0.5;0.5 0.5])];
            connections = [("S1","S3",1,2)];
            @test_throws(
                ArgumentError("Destination (network name, port number) (S3, 2) not found for connection (S1,S3,1,2)."),
                JosephsonCircuits.connectS_initialize(networks,connections)
            )
        end

        begin
            networks = [("S1",Complex{Float64}[0.0 1.0 1.0;1.0 0.0 0.0]),("S2",Complex{Float64}[0.5 0.5;0.5 0.5])];
            connections = [("S1","S2",1,2)];
            @test_throws(
                ArgumentError("The sizes of the first two dimensions (2,3) of the scattering matrix S1 must be the same."),
                JosephsonCircuits.connectS_initialize(networks,connections)
            )
        end

        begin
            networks = [(:S1,Complex{Float64}[0.0 1.0;1.0 0.0]),(:S2,Complex{Float64}[0.5 0.5;0.5 0.5]),(:S1,Complex{Float64}[0.0 1.0;1.0 0.0])];
            connections = [(:S1,:S2,1,2)];
            @test_throws(
                ArgumentError("Duplicate network names detected [(networkname,count)]: [(:S1, 2)]."),
                JosephsonCircuits.connectS_initialize(networks,connections)
            )
        end

        begin
            networks = [(:S1,Complex{Float64}[0.0 1.0;1.0 0.0],[(:S1,1),(:S1,2)]),(:S2,Complex{Float64}[0.5 0.5;0.5 0.5],[(:S2,1),(:S2,2)]),(:S3,Complex{Float64}[0.0 1.0;1.0 0.0],[(:S3,1),(:S3,2)])];
            connections = [(:S1,:S2,1,2),(:S1,:S3,1,2)];
            @test_throws(
                ArgumentError("Duplicate connections detected [(networkname,port),counts]: [((:S1, 1), 2)]."),
                JosephsonCircuits.connectS_initialize(networks,connections)
            )
        end

        begin
            networks = [(:S1,Complex{Float64}[0.0 1.0;1.0 0.0;;;0.0 1.0;1.0 0.0]),(:S2,Complex{Float64}[0.5 0.5;0.5 0.5]),(:S1,Complex{Float64}[0.0 1.0;1.0 0.0])];
            connections = [(:S1,:S2,1,2)];
            @test_throws(
                ArgumentError("The sizes of the third and higher dimensions of the scattering matrices must be the same. Size of S1 is (2, 2, 2) and size of S2 is (2, 2)."),
                JosephsonCircuits.connectS_initialize(networks,connections)
            )
        end

        begin
            networks = [(:S1,[0 1;1 0]),(:S2,[0.5 0.5;0.5 0.5])];
            connections = [(:S1,:S2,1,2)];
            @test_throws(
                ArgumentError("The element types of the scattering matrices must be the same. Element type of S1 is Int64 and element type of S2 is Float64."),
                JosephsonCircuits.connectS_initialize(networks,connections)
            )
        end

        begin
            networks = [(:S1,Complex{Float64}[0.0 1.0;1.0 0.0]),(:S2,Complex{Float64}[0.5 0.5;0.5 0.5],[(:S3,5),(:S3,5)])];
            connections = [(:S1,:S3,1,5)];
            @test_throws(
                ArgumentError("Duplicate port (:S3, 5) in network S2."),
                JosephsonCircuits.connectS_initialize(networks,connections)
            )
        end

    end

    @testset "parse_connections_sparse errors" begin

        begin
            networks = [(:S1,Complex{Float64}[0.0 1.0;1.0 0.0],[(:S1,1),(:S1,2)]),(:S2,Complex{Float64}[0.5 0.5;0.5 0.5],[(:S2,1),(:S2,2)])];
            connections = [(:S3,:S2,1,2)];
            @test_throws(
                ArgumentError("Source (network name, port number) (S3, 1) not found for connection (S3,S2,1,2)."),
                JosephsonCircuits.parse_connections_sparse(networks,connections)
            )
        end

        begin
            networks = [(:S1,Complex{Float64}[0.0 1.0;1.0 0.0],[(:S1,1),(:S1,2)]),(:S2,Complex{Float64}[0.5 0.5;0.5 0.5],[(:S2,1),(:S2,2)])];
            connections = [(:S1,:S3,1,2)];
            @test_throws(
                ArgumentError("Destination (network name, port number) (S3, 2) not found for connection (S1,S3,1,2)."),
                JosephsonCircuits.parse_connections_sparse(networks,connections)
            )
        end

        begin
            networks = [(:S1,Complex{Float64}[0.0 1.0 0.0;1.0 0.0 0.0],[(:S1,1),(:S1,2)]),(:S2,Complex{Float64}[0.5 0.5;0.5 0.5],[(:S2,1),(:S2,2)])];
            connections = [(:S1,:S2,1,2)];
            @test_throws(
                ArgumentError("The sizes of the first two dimensions (2,3) of the scattering matrix S1 must be the same."),
                JosephsonCircuits.parse_connections_sparse(networks,connections)
            )
        end

        begin
            networks = [(:S1,Complex{Float64}[0.0 1.0;1.0 0.0],[(:S1,1),(:S1,2)]),(:S2,Complex{Float64}[0.5 0.5;0.5 0.5],[(:S2,1),(:S2,2)]),(:S1,Complex{Float64}[0.0 1.0;1.0 0.0],[(:S1,1),(:S1,2)])];
            connections = [(:S1,:S2,1,2)];
            @test_throws(
                ArgumentError("Duplicate network names detected [(networkname,count)]: [(:S1, 2)]."),
                JosephsonCircuits.parse_connections_sparse(networks,connections)
            )
        end

        begin
            networks = [(:S1,Complex{Float64}[0.0 1.0;1.0 0.0],[(:S1,1),(:S1,2)]),(:S2,Complex{Float64}[0.5 0.5;0.5 0.5],[(:S2,1),(:S2,2)]),(:S3,Complex{Float64}[0.0 1.0;1.0 0.0],[(:S3,1),(:S3,2)])];
            connections = [(:S1,:S2,1,2),(:S1,:S3,1,2)];
            @test_throws(
                ArgumentError("Duplicate connections detected [(networkname,port),counts]: [((:S1, 1), 2)]."),
                JosephsonCircuits.parse_connections_sparse(networks,connections)
            )
        end

        begin
            networks = [(:S1,Complex{Float64}[0.0 1.0;1.0 0.0;;;0.0 1.0;1.0 0.0],[(:S1,1),(:S1,2)]),(:S2,Complex{Float64}[0.5 0.5;0.5 0.5],[(:S2,1),(:S2,2)])];
            connections = [(:S1,:S2,1,2)];
            @test_throws(
                ArgumentError("The sizes of the third and higher dimensions of the scattering matrices must be the same. Size of S1 is (2, 2, 2) and size of S2 is (2, 2)."),
                JosephsonCircuits.parse_connections_sparse(networks,connections)
            )
        end

        begin
            networks = [(:S1,[0 1;1 0],[(:S1,1),(:S1,2)]),(:S2,[0.5 0.5;0.5 0.5],[(:S2,1),(:S2,2)])];
            connections = [(:S1,:S2,1,2)];
            @test_throws(
                ArgumentError("The element types of the scattering matrices must be the same. Element type of S1 is Int64 and element type of S2 is Float64."),
                JosephsonCircuits.parse_connections_sparse(networks,connections)
            )
        end

        begin
            networks = [(:S1,Complex{Float64}[0.0 1.0;1.0 0.0],[(:S1,1),(:S1,2)]),(:S2,Complex{Float64}[0.5 0.5;0.5 0.5],[(:S3,5),(:S3,5)])];
            connections = [(:S1,:S3,1,5)];
            @test_throws(
                ArgumentError("Duplicate port (:S3, 5) in network S2."),
                JosephsonCircuits.parse_connections_sparse(networks,connections)
            )
        end
    end

    @testset "add_splitters errors" begin
        begin
            networks = [(:S1,Complex{Float64}[0.0 1.0;1.0 0.0],[(:S1,1),(:S1,2)]),(:S2,Complex{Float64}[0.5 0.5;0.5 0.5],[(:S2,1),(:S2,2)])];
            connections = [[(:S1,1)]];
            @test_throws(
                ArgumentError("Invalid connection [(:S1, 1)] with only network and port."),
                JosephsonCircuits.add_splitters(networks,connections)
            )
        end

        begin
            networks = [(:S1,Complex{Float64}[0.0 1.0 0.0;1.0 0.0 0.0],[(:S1,1),(:S1,2)]),(:S2,Complex{Float64}[0.5 0.5;0.5 0.5],[(:S2,1),(:S2,2)])];
            connections = [[(:S1,1),(:S2,2)]];
            @test_throws(
                ArgumentError("The sizes of the first two dimensions (2,3) of the scattering matrix S1 must be the same."),
                JosephsonCircuits.add_splitters(networks,connections)
            )
        end

        begin
            networks = [(:S1,Complex{Float64}[0.0 1.0;1.0 0.0;;;0.0 1.0;1.0 0.0],[(:S1,1),(:S1,2)]),(:S2,Complex{Float64}[0.5 0.5;0.5 0.5],[(:S2,1),(:S2,2)])];
            connections = [[(:S1,1),(:S2,2)]];
            @test_throws(
                ArgumentError("The sizes of the third and higher dimensions of the scattering matrices must be the same. Size of S1 is (2, 2, 2) and size of S2 is (2, 2)."),
                JosephsonCircuits.add_splitters(networks,connections)
            )
        end

        begin
            networks = [(:S1,[0 1;1 0],[(:S1,1),(:S1,2)]),(:S2,[0.5 0.5;0.5 0.5],[(:S2,1),(:S2,2)])];
            connections = [[(:S1,1),(:S2,2)]];
            @test_throws(
                ArgumentError("The element types of the scattering matrices must be the same. Element type of S1 is Int64 and element type of S2 is Float64."),
                JosephsonCircuits.add_splitters(networks,connections)
            )
        end

    end

    @testset "connectS solveS comparison" begin

        networks =[("S1",rand(Complex{Float64},4,4,10)),("S2",rand(Complex{Float64},3,3,10))];
        connections = [[("S1",1),("S1",2),("S1",3)],[("S1",4),("S2",2)]];
        out1 = JosephsonCircuits.connectS(networks,connections)
        out2 = JosephsonCircuits.solveS(networks,connections)
        @test isapprox(out1[1][1],out2[1])
    end

    @testset "connectS! solveS! in-place updates" begin
        S1 = rand(Complex{Float64},4,4,10)
        S2 = rand(Complex{Float64},3,3,10)
        networks =[("S1",S1),("S2",S2)];
        connections = [[("S1",1),("S1",2),("S1",3)],[("S1",4),("S2",2)]];
        init1 = JosephsonCircuits.connectS_initialize(networks,connections)
        init2 = JosephsonCircuits.solveS_initialize(networks,connections)

        S1a = JosephsonCircuits.connectS!(init1...)[1][1]
        S2a = copy(JosephsonCircuits.solveS!(init2...)[1])
        @test isapprox(S1a,S2a)

        # update the S1 scattering matrix. check that the two solvers give
        # the same results and that it's different from the initial solution
        S1 .= rand(Complex{Float64},4,4,10)
        S1b = JosephsonCircuits.connectS!(init1...)[1][1]
        S2b = copy(JosephsonCircuits.solveS!(init2...)[1])
        @test isapprox(S1b,S2b)
        @test !isapprox(S1a,S1b)
        @test !isapprox(S2a,S2b)
    end

    @testset "connectS solveS splitters" begin
        # tests from https://github.com/scikit-rf/scikit-rf/issues/1221

        # first test of connecting two splitters
        S_splitter = JosephsonCircuits.S_splitter!(zeros(Complex{Float64},3,3))

        networks = [("S1_splitter",S_splitter),("S2_splitter",S_splitter)]
        connections = [
            [("S1_splitter",2),("S2_splitter",2)],
            [("S1_splitter",3),("S2_splitter",3)],
        ]
        sol1 = JosephsonCircuits.solveS(networks,connections;factorization=JosephsonCircuits.LUfactorization())
        sol2 = JosephsonCircuits.connectS(networks,connections)
        sol3 = Complex{Float64}[0 1;1 0]

        @test isapprox(sol1[1],sol3)
        @test isapprox(sol2[1][1],sol3)

        # second test of connecting two splitters with through lines
        S_thru = JosephsonCircuits.S_splitter!(zeros(Complex{Float64},2,2))
        networks = [
            ("S1",S_thru),("S2",S_thru),("S3",S_thru),
            ("S4",S_thru),("S5",S_thru),("S6",S_thru),
        ]
        connections = [
            [("S1",2),("S2",1),("S3",1)],
            [("S2",2),("S4",1)],
            [("S3",2),("S5",1)],
            [("S4",2),("S5",2),("S6",1)],
        ]
        sol2 = JosephsonCircuits.connectS(networks,connections)
        sol3 = Complex{Float64}[0 1;1 0]
        @test isapprox(sol2[1][1],sol3)
    end

    # @testset "connectS solveS mirror" begin
    #     # test from https://github.com/scikit-rf/scikit-rf/issues/1221

    #     S1 = rand(Complex{Float64},3,3)
    #     S1_mirror = inv(S1)

    #     networks = [("S1",S1),("S1_mirror",S1_mirror)]
    #     connections = [[("S1",2),("S1_mirror",2)],[("S1",3),("S1_mirror",3)]]
        
    #     ## solveS sometimes gives singular matrix errors, so don't test on this
    #     ## network.
    #     # sol1 = JosephsonCircuits.solveS(networks,connections;factorization=JosephsonCircuits.LUfactorization())
    #     sol2 = JosephsonCircuits.connectS(networks,connections)
    #     sol3 = Complex{Float64}[0 1;1 0]

    #     # @test isapprox(sol1[1],sol3)
    #     @test isapprox(sol2[1][1],sol3)
    # end

end