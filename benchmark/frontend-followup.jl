# Run with the same options as frontend.jl; this also measures that suite.
include(joinpath(@__DIR__, "frontend.jl"))

function interfaceinput(n)
    components = Pair{Symbol,Any}[]
    connections = Any[]
    pins = Pair{Symbol,Tuple{Symbol,Int}}[]
    ports = Pair{Symbol,Tuple{Symbol,typeof(Ground)}}[]
    for i in 1:n
        r, k, w = Symbol(:r,i), Symbol(:k,i), Symbol(:w,i)
        push!(components, r => Resistor(50.0))
        push!(connections, [(r,2), Ground])
        push!(pins, k => (r,1))
        push!(ports, w => (k, Ground))
    end
    return (components, connections, Interface(pins, ports))
end
buildinterface(x) = Circuit(x...)
function followup(label)
    # Stress the dense root table when nearly every wire is ground.
    n4096 = benchsize(4096)
    input = Any[(Symbol(:r, i), 0, 0, Resistor(50.0)) for i in 1:n4096]
    c = record(label, "allground$(n4096)", "construct", construct, input)
    record(label, "allground$(n4096)", "elaborate", elaborate, c)
    for n in (2, 16, benchsize(1024))
        name = "interface$n"
        c = record(label, name, "construct", buildinterface, interfaceinput(n))
        connections = [PortRef(:a, Symbol(:w,i)) => PortRef(:b, Symbol(:w,i)) for i in 1:n]
        input = ([:a => c, :b => c], connections)
        top = record(label, "ports$n", "construct", construct, input)
        record(label, "ports$n", "parse", JosephsonCircuits.parsecircuitlevel, top)
        record(label, "ports$n", "elaborate", elaborate, top)
    end
end
followup(BENCH_LABEL)
