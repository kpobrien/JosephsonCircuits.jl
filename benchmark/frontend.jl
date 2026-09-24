# Run in a fresh process for each revision, with the same Julia and depot:
# julia --startup-file=no --project=. benchmark/frontend.jl LABEL
# Inputs are prepared before timing; package loading and input generation
# are excluded. First-call figures include frontend JIT compilation.
using JosephsonCircuits
using Statistics
include(joinpath(@__DIR__, "options.jl"))

const N4096 = benchsize(4096)
const N1024 = benchsize(1024)
const N128 = benchsize(128)

function jpa()
    Any[(:p, 1, 0, Port(1)), (:cc, 1, 2, Capacitor(1e-13)),
        (:jj, 2, 0, JosephsonJunction(1e-9)), (:cj, 2, 0, Capacitor(1e-12))]
end
function ladder(n)
    entries = Any[(:p, 1, 0, Port(1))]
    for i in 1:n
        push!(entries, (Symbol(:l, i), i, i+1, Inductor(1e-10)))
        push!(entries, (Symbol(:c, i), i+1, 0, Capacitor(1e-13)))
    end
    entries
end
function groupinput(n)
    components = Pair{Symbol,Any}[:p => Port(1), :gnd => Ground()]
    connections = Vector{Tuple{Symbol,Int}}[]
    ground = Tuple{Symbol,Int}[(:p, 2), (:gnd, 1)]
    previous = (:p, 1)
    for i in 1:n
        l, c = Symbol(:l, i), Symbol(:c, i)
        push!(components, l => Inductor(1e-10), c => Capacitor(1e-13))
        push!(connections, [previous, (l, 1)], [(l, 2), (c, 1)])
        push!(ground, (c, 2))
        previous = (l, 2)
    end
    push!(connections, ground)
    return (components, connections)
end
construct(input::AbstractVector) = Circuit(input)
construct(input::Tuple) = Circuit(input...)
function legacy(n)
    entries = Tuple{String,String,String,Any}[("P1", "1", "0", 1),
        ("R1", "0", "1", 50.0)]
    for i in 1:n
        push!(entries, ("L$i", string(i), string(i+1), 1e-10))
        push!(entries, ("C$i", string(i+1), "0", 1e-13))
    end
    entries
end
function hierarchy(n)
    cell = Circuit([(:l, 1, 2, Inductor(1e-10)),
        (:c, 2, 0, Capacitor(1e-13))];
        pins = [:in => (:l, 1), :out => (:l, 2)])
    Any[(:p, 1, 0, Port(1));
        [(Symbol(:cell, i), i, i+1, cell) for i in 1:n]]
end

function record(label, case, stage, f, x)
    GC.gc()
    first = @timed f(x)
    f(x) # one additional warm-up, outside the sample set
    samples = [@timed(f(x)) for _ in 1:BENCH_SAMPLES]
    coldcompile = hasproperty(first, :compile_time) ? first.compile_time : NaN
    println(join((label, case, stage, first.time, coldcompile,
        median(s.time for s in samples), median(s.bytes for s in samples),
        median(s.gcstats.malloc + s.gcstats.realloc + s.gcstats.poolalloc +
            s.gcstats.bigalloc for s in samples),
        Base.summarysize(first.value)), '\t'))
    return first.value
end

function main(label)
    # Prepare every raw input before the first measured constructor.
    cases = [("jpa", jpa()), ("ladder$(N128)", ladder(N128)),
        ("ladder$(N4096)", ladder(N4096)), ("groups$(N4096)", groupinput(N4096)),
        ("legacy$(N4096)", legacy(N4096))]
    println("label\tcase\tstage\tfirst_s\tfirst_compile_s\twarm_s\tbytes\tallocations\tresult_bytes")
    for (name, input) in cases
        c = record(label, name, "construct", construct, input)
        record(label, name, "parse", JosephsonCircuits.parsecircuitlevel, c)
        e = record(label, name, "elaborate", elaborate, c)
        record(label, name, "lower", compile, e)
    end
    input = hierarchy(N1024)
    c = record(label, "hierarchy$(N1024)", "construct", construct, input)
    record(label, "hierarchy$(N1024)", "parse", JosephsonCircuits.parsecircuitlevel, c)
    e = record(label, "hierarchy$(N1024)", "elaborate", elaborate, c)
    record(label, "hierarchy$(N1024)", "lower", compile, e)
end
main(BENCH_LABEL)
