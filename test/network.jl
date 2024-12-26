using JosephsonCircuits
using LinearAlgebra
using Test
import StaticArrays

@testset verbose=true "network" begin

    # one network
    @testset "connectS" begin
        Sa = Float64[1 2;3 4]
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

    @testset "connectS" begin
        Sa = Float64[1 2 3;4 5 6]
        @test_throws(
            DimensionMismatch("Lengths of first two dimensions of `Sa` must be equal."),
            JosephsonCircuits.connectS(Sa,1,2)
        )
    end

    # one network, in place
    @testset "connectS!" begin
        Sa = rand(Complex{Float64},3,3)
        Sout = zeros(Complex{Float64},2,2)
        @test_throws(
            DimensionMismatch("Length of first two dimensions must be 2 smaller for `Sout` than `Sa` because we are merging two ports."),
            JosephsonCircuits.connectS!(Sout,Sa,1,2)
        )
    end

    @testset "connectS!" begin
        Sa = rand(Complex{Float64},4,4)
        Sout = zeros(Complex{Float64},2,2,1)
        @test_throws(
            DimensionMismatch("`Sout` and `Sa` must have the same number of dimensions."),
            JosephsonCircuits.connectS!(Sout,Sa,1,2)
        )
    end

    @testset "connectS!" begin
        Sa = rand(Complex{Float64},4)
        Sout = zeros(Complex{Float64},2)
        @test_throws(
            DimensionMismatch("`Sout` and `Sa` must have at least two dimensions."),
            JosephsonCircuits.connectS!(Sout,Sa,1,2)
        )
    end

    @testset "connectS!" begin
        Sa = rand(Complex{Float64},4,4)
        Sout = zeros(Complex{Float64},2,3)
        @test_throws(
            DimensionMismatch("Lengths of first two dimensions of `Sout` must be equal."),
            JosephsonCircuits.connectS!(Sout,Sa,1,2)
        )
    end

    @testset "connectS!" begin
        Sa = rand(Complex{Float64},3,3,3)
        Sout = zeros(Complex{Float64},1,1,4)
        @test_throws(
            DimensionMismatch("Non-port axis lengths of `Sa` and `Sout` must be equal."),
            JosephsonCircuits.connectS!(Sout,Sa,1,2)
        )
    end

    # two networks
    @testset "connectS" begin
        Sa = Float64[1 2;3 4]
        Sb = Float64[5 6;7 8]
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

    @testset "connectS" begin
        Sa = Float64[1 2 3;4 5 6]
        Sb = Float64[5 6;7 8]
        @test_throws(
            DimensionMismatch("Lengths of first two dimensions of `Sa` must be equal."),
            JosephsonCircuits.connectS(Sa,Sb,1,2)
        )
    end

    @testset "connectS" begin
        Sa = Float64[1 2;4 5]
        Sb = Float64[5 6 7;8 9 10]
        @test_throws(
            DimensionMismatch("Lengths of first two dimensions of `Sb` must be equal."),
            JosephsonCircuits.connectS(Sa,Sb,1,2)
        )
    end

    # two networks, in place
    @testset "connectS!" begin
        Sa = rand(Complex{Float64},3,3)
        Sb = rand(Complex{Float64},3,3)
        Sout = zeros(Complex{Float64},2,2)
        @test_throws(
            DimensionMismatch("First two dimensions of `Sout` must be `m+n-2`."),
            JosephsonCircuits.connectS!(Sout,Sa,Sb,1,2)
        )
    end

    @testset "connectS!" begin
        Sa = rand(Complex{Float64},3,3)
        Sb = rand(Complex{Float64},3,3,1)
        Sout = zeros(Complex{Float64},4,4)
        @test_throws(
            DimensionMismatch("`Sa` and `Sb` must have the same number of dimensions."),
            JosephsonCircuits.connectS!(Sout,Sa,Sb,1,2)
        )
    end

    @testset "connectS!" begin
        Sa = rand(Complex{Float64},4,4)
        Sb = rand(Complex{Float64},4,4)
        Sout = zeros(Complex{Float64},6,6,1)
        @test_throws(
            DimensionMismatch("`Sout`, `Sa`, and `Sb` must have the same number of dimensions."),
            JosephsonCircuits.connectS!(Sout,Sa,Sb,1,2)
        )
    end

    @testset "connectS!" begin
        Sa = rand(Complex{Float64},4)
        Sb = rand(Complex{Float64},4)
        Sout = zeros(Complex{Float64},6)
        @test_throws(
            DimensionMismatch("`Sout`, `Sa`, and `Sb` must have atleast two dimensions."),
            JosephsonCircuits.connectS!(Sout,Sa,Sb,1,2)
        )
    end

    @testset "connectS!" begin
        Sa = rand(Complex{Float64},4,4)
        Sb = rand(Complex{Float64},4,4)
        Sout = zeros(Complex{Float64},6,7)
        @test_throws(
            DimensionMismatch("Lengths of first two dimensions of `Sout` must be equal."),
            JosephsonCircuits.connectS!(Sout,Sa,Sb,1,2)
        )
    end

    @testset "connectS!" begin
        Sa = rand(Complex{Float64},3,3,3)
        Sb = rand(Complex{Float64},3,3,3)
        Sout = zeros(Complex{Float64},4,4,4)
        @test_throws(
            DimensionMismatch("Non-port axis lengths of `Sa`, `Sb`, and `Sout` must be equal."),
            JosephsonCircuits.connectS!(Sout,Sa,Sb,1,2)
        )
    end

    # check for consistency between the one and two network functions
    @testset "connectS consistency" begin
        Sa = rand(Complex{Float64},3,3,3)
        Sb = rand(Complex{Float64},3,3,3)
        Sboth = zeros(Complex{Float64},6,6,3)
        Sboth[1:size(Sa,1),1:size(Sa,1),:,:] .= Sa
        Sboth[size(Sa,1)+1:end,size(Sa,1)+1:end,:,:] .= Sb
        Sout1 = JosephsonCircuits.connectS(Sboth,1,2+size(Sa,1))
        Sout2 = JosephsonCircuits.connectS(Sa,Sb,1,2)
        @test isapprox(Sout1,Sout2)
    end

    @testset "connectS consistency" begin
        Sa = rand(Complex{Float64},3,3)
        Sb = rand(Complex{Float64},3,3)
        Sboth = zeros(Complex{Float64},6,6)
        Sboth[1:size(Sa,1),1:size(Sa,1),:,:] .= Sa
        Sboth[size(Sa,1)+1:end,size(Sa,1)+1:end,:,:] .= Sb
        Sout1 = JosephsonCircuits.connectS(Sboth,2,1+size(Sa,1))
        Sout2 = JosephsonCircuits.connectS(Sa,Sb,2,1)
        @test isapprox(Sout1,Sout2)
    end

    @testset "connectS consistency StaticArrays" begin
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

    @testset "connectS consistency StaticArrays" begin
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

    # one network
    @testset "connectSports" begin
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
    @testset "connectSports" begin
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

        # with symbols
        networks = [(:S1,S1),(:S2,S2),(:S3,Ssplitter),(:S4,Sopen)]
        connections = [(:S1,:S1,1,2),(:S1,:S2,3,1),(:S3,:S2,2,2),(:S3,:S4,3,1)]
        Sout1, ports = JosephsonCircuits.connectS(networks,connections)

        Sout2 = begin
            S = JosephsonCircuits.connectS(Ssplitter,Sopen,3,1)
            S = JosephsonCircuits.connectS(S,S2,2,2)
            S = JosephsonCircuits.connectS(S1,S,3,2)
            S = JosephsonCircuits.connectS(S,1,2)
        end
        @test isapprox(Sout1[1],Sout2)

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

        networks = [(:S1,[0.0 1.0;1.0 0.0]),(:S2,[0.5 0.5;0.5 0.5])];
        connections = [(:S1,:S2,1,2)];
        userinput = ones(Bool,length(networks))
        storage = Dict{Int,typeof(networks[1][2])}()
        g, fconnectionlist, fweightlist, ports, networkdata = JosephsonCircuits.connectS_initialize(networks,connections)
        # corrupt the vector of ports
        ports = [[(:S1, 1), (:S1, 2)],[(:S2, 1), (:S3, 2)]]
        @test_throws(
            ArgumentError("Destination port (:S2, 2) not found in the ports [(:S2, 1), (:S3, 2)] of the destination node 2."),
            JosephsonCircuits.make_connection!(g, fconnectionlist, fweightlist, ports, networkdata,1,1,1,userinput,storage)
        )

        networks = [(:S1,[0.0 1.0;1.0 0.0]),(:S2,[0.5 0.5;0.5 0.5])];
        connections = [(:S1,:S2,1,2)];
        userinput = ones(Bool,length(networks))
        storage = Dict{Int,typeof(networks[1][2])}()
        g, fconnectionlist, fweightlist, ports, networkdata = JosephsonCircuits.connectS_initialize(networks,connections)
        # corrupt the vector of ports
        ports = [[(:S3, 1), (:S1, 2)],[(:S2, 1), (:S2, 2)]]
        @test_throws(
            ArgumentError("Source port (:S1, 1) not found in the ports [(:S3, 1), (:S1, 2)] of the source node 1."),
            JosephsonCircuits.make_connection!(g, fconnectionlist, fweightlist, ports, networkdata,1,1,1,userinput,storage)
        )

    end

    @testset "connectS_initialize errors" begin

        networks = [(:S1,[0.0 1.0;1.0 0.0]),(:S2,[0.5 0.5;0.5 0.5])];
        connections = [(:S3,:S2,1,2)];
        @test_throws(
            ArgumentError("Source (network name, port number) (:S3, 1) not found for connection (S3,S2,1,2)."),
            JosephsonCircuits.connectS_initialize(networks,connections)
        )

        networks = [(:S1,[0.0 1.0;1.0 0.0]),(:S2,[0.5 0.5;0.5 0.5])];
        connections = [(:S1,:S3,1,2)];
        @test_throws(
            ArgumentError("Destination (network name, port number) (:S3, 2) not found for connection (S1,S3,1,2)."),
            JosephsonCircuits.connectS_initialize(networks,connections)
        )

        networks = [(:S1,[0.0 1.0 1.0;1.0 0.0 0.0]),(:S2,[0.5 0.5;0.5 0.5])];
        connections = [(:S1,:S2,1,2)];
        @test_throws(
            ArgumentError("The sizes of the first two dimensions (2,3) of the scattering matrix S1 must be the same."),
            JosephsonCircuits.connectS_initialize(networks,connections)
        )

        networks = [(:S1,[0.0 1.0;1.0 0.0]),(:S2,[0.5 0.5;0.5 0.5]),(:S1,[0.0 1.0;1.0 0.0])];
        connections = [(:S1,:S2,1,2)];
        @test_throws(
            ArgumentError("Duplicate network names detected [(networkname,count)]: [(:S1, 2)]."),
            JosephsonCircuits.connectS_initialize(networks,connections)
        )

        networks = [(:S1,[0.0 1.0;1.0 0.0],[(:S1,1),(:S1,2)]),(:S2,[0.5 0.5;0.5 0.5],[(:S2,1),(:S2,2)]),(:S3,[0.0 1.0;1.0 0.0],[(:S3,1),(:S3,2)])];
        connections = [(:S1,:S2,1,2),(:S1,:S3,1,2)];
        @test_throws(
            ArgumentError("Duplicate connections detected [(networkname,port),counts]: [((:S1, 1), 2)]."),
            JosephsonCircuits.connectS_initialize(networks,connections)
        )

        networks = [(:S1,[0.0 1.0;1.0 0.0;;;0.0 1.0;1.0 0.0]),(:S2,[0.5 0.5;0.5 0.5]),(:S1,[0.0 1.0;1.0 0.0])];
        connections = [(:S1,:S2,1,2)];
        @test_throws(
            ArgumentError("The sizes of the third and higher dimensions of the scattering matrices must be the same. Size of S1 is (2, 2, 2) and size of S2 is (2, 2)."),
            JosephsonCircuits.connectS_initialize(networks,connections)
        )

        networks = [(:S1,[0 1;1 0]),(:S2,[0.5 0.5;0.5 0.5])];
        connections = [(:S1,:S2,1,2)];
        @test_throws(
            ArgumentError("The element types of the scattering matrices must be the same. Element type of S1 is Int64 and element type of S2 is Float64."),
            JosephsonCircuits.connectS_initialize(networks,connections)
        )

        networks = [(:S1,[0.0 1.0;1.0 0.0]),(:S2,[0.5 0.5;0.5 0.5],[(:S3,5),(:S3,5)])];
        connections = [(:S1,:S3,1,5)];
        @test_throws(
            ArgumentError("Duplicate port (:S3, 5) in network S2."),
            JosephsonCircuits.connectS_initialize(networks,connections)
        )

    end

    @testset "parse_connections_sparse errors" begin

        networks = [(:S1,[0.0 1.0;1.0 0.0],[(:S1,1),(:S1,2)]),(:S2,[0.5 0.5;0.5 0.5],[(:S2,1),(:S2,2)])];
        connections = [(:S3,:S2,1,2)];
        @test_throws(
            ArgumentError("Source (network name, port number) (:S3, 1) not found for connection (S3,S2,1,2)."),
            JosephsonCircuits.parse_connections_sparse(networks,connections)
        )

        networks = [(:S1,[0.0 1.0;1.0 0.0],[(:S1,1),(:S1,2)]),(:S2,[0.5 0.5;0.5 0.5],[(:S2,1),(:S2,2)])];
        connections = [(:S1,:S3,1,2)];
        @test_throws(
            ArgumentError("Destination (network name, port number) (:S3, 2) not found for connection (S1,S3,1,2)."),
            JosephsonCircuits.parse_connections_sparse(networks,connections)
        )

        networks = [(:S1,[0.0 1.0 0.0;1.0 0.0 0.0],[(:S1,1),(:S1,2)]),(:S2,[0.5 0.5;0.5 0.5],[(:S2,1),(:S2,2)])];
        connections = [(:S1,:S2,1,2)];
        @test_throws(
            ArgumentError("The sizes of the first two dimensions (2,3) of the scattering matrix S1 must be the same."),
            JosephsonCircuits.parse_connections_sparse(networks,connections)
        )

        networks = [(:S1,[0.0 1.0;1.0 0.0],[(:S1,1),(:S1,2)]),(:S2,[0.5 0.5;0.5 0.5],[(:S2,1),(:S2,2)]),(:S1,[0.0 1.0;1.0 0.0],[(:S1,1),(:S1,2)])];
        connections = [(:S1,:S2,1,2)];
        @test_throws(
            ArgumentError("Duplicate network names detected [(networkname,count)]: [(:S1, 2)]."),
            JosephsonCircuits.parse_connections_sparse(networks,connections)
        )

        networks = [(:S1,[0.0 1.0;1.0 0.0],[(:S1,1),(:S1,2)]),(:S2,[0.5 0.5;0.5 0.5],[(:S2,1),(:S2,2)]),(:S3,[0.0 1.0;1.0 0.0],[(:S3,1),(:S3,2)])];
        connections = [(:S1,:S2,1,2),(:S1,:S3,1,2)];
        @test_throws(
            ArgumentError("Duplicate connections detected [(networkname,port),counts]: [((:S1, 1), 2)]."),
            JosephsonCircuits.parse_connections_sparse(networks,connections)
        )

        networks = [(:S1,[0.0 1.0;1.0 0.0;;;0.0 1.0;1.0 0.0],[(:S1,1),(:S1,2)]),(:S2,[0.5 0.5;0.5 0.5],[(:S2,1),(:S2,2)])];
        connections = [(:S1,:S2,1,2)];
        @test_throws(
            ArgumentError("The sizes of the third and higher dimensions of the scattering matrices must be the same. Size of S1 is (2, 2, 2) and size of S2 is (2, 2)."),
            JosephsonCircuits.parse_connections_sparse(networks,connections)
        )

        networks = [(:S1,[0 1;1 0],[(:S1,1),(:S1,2)]),(:S2,[0.5 0.5;0.5 0.5],[(:S2,1),(:S2,2)])];
        connections = [(:S1,:S2,1,2)];
        @test_throws(
            ArgumentError("The element types of the scattering matrices must be the same. Element type of S1 is Int64 and element type of S2 is Float64."),
            JosephsonCircuits.parse_connections_sparse(networks,connections)
        )

        networks = [(:S1,[0.0 1.0;1.0 0.0],[(:S1,1),(:S1,2)]),(:S2,[0.5 0.5;0.5 0.5],[(:S3,5),(:S3,5)])];
        connections = [(:S1,:S3,1,5)];
        @test_throws(
            ArgumentError("Duplicate port (:S3, 5) in network S2."),
            JosephsonCircuits.parse_connections_sparse(networks,connections)
        )

    end

    @testset "add_splitters errors" begin
        networks = [(:S1,[0.0 1.0;1.0 0.0],[(:S1,1),(:S1,2)]),(:S2,[0.5 0.5;0.5 0.5],[(:S2,1),(:S2,2)])];
        connections = [[(:S1,1)]];
        @test_throws(
            ArgumentError("Invalid connection [(:S1, 1)] with only network and port."),
            JosephsonCircuits.add_splitters(networks,connections)
        )

        networks = [(:S1,[0.0 1.0 0.0;1.0 0.0 0.0],[(:S1,1),(:S1,2)]),(:S2,[0.5 0.5;0.5 0.5],[(:S2,1),(:S2,2)])];
        connections = [[(:S1,1),(:S2,2)]];
        @test_throws(
            ArgumentError("The sizes of the first two dimensions (2,3) of the scattering matrix S1 must be the same."),
            JosephsonCircuits.add_splitters(networks,connections)
        )

        networks = [(:S1,[0.0 1.0;1.0 0.0;;;0.0 1.0;1.0 0.0],[(:S1,1),(:S1,2)]),(:S2,[0.5 0.5;0.5 0.5],[(:S2,1),(:S2,2)])];
        connections = [[(:S1,1),(:S2,2)]];
        @test_throws(
            ArgumentError("The sizes of the third and higher dimensions of the scattering matrices must be the same. Size of S1 is (2, 2, 2) and size of S2 is (2, 2)."),
            JosephsonCircuits.add_splitters(networks,connections)
        )

        networks = [(:S1,[0 1;1 0],[(:S1,1),(:S1,2)]),(:S2,[0.5 0.5;0.5 0.5],[(:S2,1),(:S2,2)])];
        connections = [[(:S1,1),(:S2,2)]];
        @test_throws(
            ArgumentError("The element types of the scattering matrices must be the same. Element type of S1 is Int64 and element type of S2 is Float64."),
            JosephsonCircuits.add_splitters(networks,connections)
        )

    end

    @testset "S_splitter! errors" begin
        @test_throws(
            ArgumentError("The sizes of the first two dimensions (3,2) of the scattering matrix must be the same."),
            JosephsonCircuits.S_splitter!(ones(3,2))
        )
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

    @testset "connectS solveS mirror" begin
        # test from https://github.com/scikit-rf/scikit-rf/issues/1221

        S1 = rand(Complex{Float64},3,3)
        S1_mirror = inv(S1)

        networks = [("S1",S1),("S1_mirror",S1_mirror)]
        connections = [[("S1",2),("S1_mirror",2)],[("S1",3),("S1_mirror",3)]]
        
        ## solveS sometimes gives singular matrix errors, so don't test on this
        ## network.
        # sol1 = JosephsonCircuits.solveS(networks,connections;factorization=JosephsonCircuits.LUfactorization())
        sol2 = JosephsonCircuits.connectS(networks,connections)
        sol3 = Complex{Float64}[0 1;1 0]

        # @test isapprox(sol1[1],sol3)
        @test isapprox(sol2[1][1],sol3)
    end

    @testset "StoZ, StoY, StoA, StoB, StoABCD consistency" begin
        # the different functions we want to test
        for f in [
                (JosephsonCircuits.ZtoS,JosephsonCircuits.StoZ),
                (JosephsonCircuits.YtoS,JosephsonCircuits.StoY),
                (JosephsonCircuits.AtoS,JosephsonCircuits.StoA),
                (JosephsonCircuits.BtoS,JosephsonCircuits.StoB),
                (JosephsonCircuits.ABCDtoS,JosephsonCircuits.StoABCD),
            ]
            # single matrix input
            for portimpedances in [
                    rand(Complex{Float64}), rand(Complex{Float64},2),
                    (StaticArrays.@MVector rand(Complex{Float64},2))
                ]
                for arg1 in [rand(Complex{Float64},2,2), (StaticArrays.@MMatrix rand(Complex{Float64},2,2))]
                    arg2 = f[1](arg1,portimpedances=portimpedances)
                    arg3 = f[2](arg2,portimpedances=portimpedances)
                    @test isapprox(arg1,arg3)
                end
            end
            # array input
            for portimpedances in [rand(Complex{Float64}), rand(Complex{Float64},2,10)]
                for arg1 in [rand(Complex{Float64},2,2,10)]
                    arg2 = f[1](arg1,portimpedances=portimpedances)
                    arg3 = f[2](arg2,portimpedances=portimpedances)
                    @test isapprox(arg1,arg3)
                end
            end
            # vector of matrices
            for portimpedances in [rand(Complex{Float64}), rand(Complex{Float64},2), (StaticArrays.@MVector rand(Complex{Float64},2))]
                for arg1 in [
                        [rand(Complex{Float64},2,2) for i in 1:10],
                        [(StaticArrays.@MMatrix rand(Complex{Float64},2,2)) for i in 1:10],
                    ]
                    arg2 = [f[1](arg1[i],portimpedances=portimpedances) for i in 1:10]
                    arg3 = [f[2](arg2[i],portimpedances=portimpedances) for i in 1:10]
                    @test isapprox(arg1,arg3)
                end
            end
        end
    end

    @testset "StoT, AtoB, ZtoA, YtoA, YtoB, ZtoB, ZtoY consistency" begin
        # the different functions we want to test
        for f in [
                (JosephsonCircuits.StoT,JosephsonCircuits.TtoS),
                (JosephsonCircuits.AtoB,JosephsonCircuits.BtoA),
                (JosephsonCircuits.ZtoA,JosephsonCircuits.AtoZ),
                (JosephsonCircuits.YtoA,JosephsonCircuits.AtoY),
                (JosephsonCircuits.YtoB,JosephsonCircuits.BtoY),
                (JosephsonCircuits.ZtoB,JosephsonCircuits.BtoZ),
                (JosephsonCircuits.ZtoY,JosephsonCircuits.YtoZ),
            ]
            # single matrix input
            for arg1 in [rand(Complex{Float64},2,2), (StaticArrays.@MMatrix rand(Complex{Float64},2,2))]
                arg2 = f[1](arg1)
                arg3 = f[2](arg2)
                @test isapprox(arg1,arg3)
            end
            # array input
            for arg1 in [rand(Complex{Float64},2,2,10)]
                arg2 = f[1](arg1)
                arg3 = f[2](arg2)
                @test isapprox(arg1,arg3)
            end
            # vector of matrices
            for arg1 in [
                    [rand(Complex{Float64},2,2) for i in 1:10],
                    [(StaticArrays.@MMatrix rand(Complex{Float64},2,2)) for i in 1:10],
                ]
                arg2 = [f[1](arg1[i]) for i in 1:10]
                arg3 = [f[2](arg2[i]) for i in 1:10]
                @test isapprox(arg1,arg3)
            end
        end
    end

    @testset "ZC_basis_coupled_lines" begin
        Cg = 1.6670474181399462e-10
        Cm = 9.320861870486729e-12
        Ls = 4.167618545349866e-7
        Lm = 2.330215467621682e-8
        C = JosephsonCircuits.LinearAlgebra.Symmetric([Cg -Cm;-Cm Cg])
        L = JosephsonCircuits.LinearAlgebra.Symmetric([Ls Lm;Lm Ls])
        b = JosephsonCircuits.ZC_basis_coupled_tlines(L,C)
        ZC = b.ZC
        TI = b.TI
        TV = b.TV
        theta = b.theta
        U = b.U
        lambda = b.lambda
        S = b.S
        @test isapprox(U'*C*U,theta^2)
        @test isapprox(S'*(theta*U'*L*U*theta)*S,lambda^2)
        @test isapprox(inv(TI)*C*L*TI,lambda^2)
        @test isapprox(TV'*TI,I)
    end

    @testset "A_coupled_tlines" begin

        Zeven = 51.0
        Zodd = 49.0
        neven = 1.1
        nodd = 1.08
        l = 3.5e-3
        c = JosephsonCircuits.speed_of_light
        omega = 2*pi*(1:10)*1e9

        L, C = JosephsonCircuits.even_odd_to_maxwell(Zeven, Zodd, neven, nodd)
        A1 = JosephsonCircuits.A_coupled_tlines(L,C,l,omega)
        A2 = stack(JosephsonCircuits.ABCD_coupled_tline.(Zeven,Zodd,neven*omega/c*l,nodd*omega/c*l))
        @test isapprox(A1,A2)

    end

    @testset "Z_canonical_coupled_line_circuits" begin
        Zeven = 52.0
        Zodd = 48.0
        neven = 2.2
        nodd = 2.1
        l = 2*3.0e-3
        c = 2.998e8
        omega = 2*pi*5e9

        function S_canonical_coupled_line_circuits(i::Int, Ze, Zo, thetae, thetao)

            # define an open
            Sopen = ones(Complex{Float64},1,1)

            # and a short
            Sshort = -ones(Complex{Float64},1,1)
            
            # and a coupled transmission line
            S = JosephsonCircuits.ZtoS(JosephsonCircuits.Z_coupled_tline(Zeven,Zodd,thetae,thetao));

            if i == 1
                S = JosephsonCircuits.connectS(S,Sopen,3,1)
                S = JosephsonCircuits.connectS(S,Sshort,1,1)
            elseif i == 2
                S = JosephsonCircuits.connectS(S,Sshort,4,1)
                S = JosephsonCircuits.connectS(S,Sshort,1,1)
            elseif i == 3
                S = JosephsonCircuits.connectS(S,Sopen,4,1)
                S = JosephsonCircuits.connectS(S,Sopen,1,1)
            elseif i == 4
                S = JosephsonCircuits.connectS(S,Sopen,4,1)
                S = JosephsonCircuits.connectS(S,Sshort,3,1)
            elseif i == 5
                S = JosephsonCircuits.connectS(S,Sopen,3,1)
                S = JosephsonCircuits.connectS(S,Sopen,1,1)
            elseif i == 6
                S = JosephsonCircuits.connectS(S,Sshort,3,1)
                S = JosephsonCircuits.connectS(S,Sshort,1,1)
            elseif i == 7
                S = JosephsonCircuits.connectS(S,4,3)
            elseif i == 8
                S = JosephsonCircuits.connectS(S,Sopen,4,1)
                S = JosephsonCircuits.connectS(S,Sshort,1,1)
            elseif i == 9
                S = JosephsonCircuits.connectS(S,Sshort,4,1)
                S = JosephsonCircuits.connectS(S,Sshort,3,1)
            elseif i == 10
                S = JosephsonCircuits.connectS(S,Sopen,4,1)
                S = JosephsonCircuits.connectS(S,Sopen,3,1)
            else
                throw(ArgumentError("Canonical coupled line circuit number must be 1-10."))
            end
                
            return S
        end

        for i in 1:10
            S1 = S_canonical_coupled_line_circuits(i, Zeven, Zodd, neven*omega/c*l, nodd*omega/c*l)
            S2 = JosephsonCircuits.ZtoS(JosephsonCircuits.Z_canonical_coupled_line_circuits(i, Zeven, Zodd, neven*omega/c*l, nodd*omega/c*l));
            @test isapprox(S1,S2)
        end

        @test_throws(
            ArgumentError("Canonical coupled line circuit number must be 1-10."),
            JosephsonCircuits.Z_canonical_coupled_line_circuits(11, Zeven, Zodd, neven*omega/c*l, nodd*omega/c*l),
        )

    end

    @testset "canonical_coupled_line_circuits" begin

        Zeven = 52.0
        Zodd = 48.0
        neven = 2.2
        nodd = 2.1

        # 3
        outa = (L1 = 1.1963832214440388e-7, L2 = 1.1963832214440388e-7, M = -3.780393078912389e-9, C1 = 1.4112327104537203e-10, C2 = 1.4112327104537203e-10, Cm = 2.4055103019097465e-12)
        outb = JosephsonCircuits.canonical_coupled_line_circuits(3, Zeven, Zodd, neven, nodd)
        @test all([isapprox(outai,outbi) for (outai,outbi) in zip(outa,outb)])

        # 8
        outa = (L1 = 1.1941849319292907e-7, L2 = 3.5891496643321166e-7, M = -8.333511920257754e-9, C1 = 1.4352878134728178e-10)
        outb = JosephsonCircuits.canonical_coupled_line_circuits(8, Zeven, Zodd, neven, nodd)
        @test all([isapprox(outai,outbi) for (outai,outbi) in zip(outa,outb)])

        # 9
        outa = (L1 = 3.5891496643321166e-7, L2 = 3.5891496643321166e-7, M = 2.268235847347433e-8)
        outb = JosephsonCircuits.canonical_coupled_line_circuits(9, Zeven, Zodd, neven, nodd)
        @test all([isapprox(outai,outbi) for (outai,outbi) in zip(outa,outb)])

        # 10
        outa = (L1 = 1.1963832214440388e-7, L2 = 1.1963832214440388e-7, M = 7.560786157824777e-9, C1 = 1.4112327104537203e-10, C2 = 1.4112327104537203e-10, Cm = 2.4055103019097465e-12)
        outb = JosephsonCircuits.canonical_coupled_line_circuits(10, Zeven, Zodd, neven, nodd)
        @test all([isapprox(outai,outbi) for (outai,outbi) in zip(outa,outb)])


        @test_throws(
            ArgumentError("Canonical coupled line circuit number must be 1-10."),
            JosephsonCircuits.canonical_coupled_line_circuits(11, Zeven, Zodd, neven, nodd),
        )

    end

end