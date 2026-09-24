# Frontend benchmark. Run in a fresh process with one Julia/BLAS thread:
# julia --startup-file=no --project=. benchmark/circuit-preparation.jl label
# Package loading and input construction are outside the timers.
using JosephsonCircuits, LinearAlgebra, SparseArrays, Statistics, Random
include(joinpath(@__DIR__, "options.jl"))
const JC = JosephsonCircuits
const label = BENCH_LABEL
const NLADDER = benchsize(4096)
const NCOUPLED = benchsize(512)
function measure(case, f, input)
    @nospecialize f input
    firstcall, compile, samples = benchtime(f, input)
    println(join((label,case,firstcall.time,compile,
        median(s.time for s in samples),median(s.bytes for s in samples)), '\t'))
    flush(stdout)
    return firstcall.value
end
function ladder(n; coupled=false)
    e = Any[(:p,1,0,Port(1))]
    for i in 1:n
        push!(e,(Symbol(:l,i),i,i+1,Inductor(1e-10)),
            (Symbol(:c,i),i+1,0,Capacitor(1e-13)))
        if coupled && iseven(i)
            push!(e,(Symbol(:k,i),Symbol(:l,i-1),Symbol(:l,i),MutualInductor(0.2)))
        end
    end
    return Circuit(e)
end
println("label\tcase\tfirst_s\tcompile_s\twarm_s\tbytes")
c = ladder(NLADDER)
function front(c)
    cc = JC.compile(c)
    return cc, numericmatrices(cc, Dict{Symbol,Any}(); Nmodes = 3)
end
bound(cc) = JC.bindvalues(cc, JC.componentvaluestonumber(cc.componentvalues, Dict{Any,Any}()))
cc,nm = measure("ladder$(NLADDER)/prepare",front,c)
b = bound(cc)
makeplan(x) = JC.circuitmatrixplan(x[1]; Nmodes = 3)
plan = measure("ladder$(NLADDER)/plan",makeplan,(cc,b))
assemble(x) = JC.assemblematrices(x...)
measure("ladder$(NLADDER)/assemble",assemble,(plan,b))
work = JC.CircuitMatrixWorkspace(plan,nm)
refill(x) = JC.assemblematrices!(x...)
measure("ladder$(NLADDER)/refill",refill,(nm,plan,b,work))
numeric(x) = numericmatrices(x...;Nmodes=3)
measure("ladder$(NLADDER)/numericmatrices",numeric,(cc,b.values))
cc2,nm2 = front(ladder(NCOUPLED;coupled=true)); b2=bound(cc2)
plan2=measure("coupled$(NCOUPLED)/plan",makeplan,(cc2,b2))
measure("coupled$(NCOUPLED)/assemble",assemble,(plan2,b2))
work2 = JC.CircuitMatrixWorkspace(plan2,nm2)
measure("coupled$(NCOUPLED)/refill",refill,(nm2,plan2,b2,work2))
rng=MersenneTwister(11926)
A=Matrix(Diagonal(-collect(range(1.,10.;length=48))))
B=randn(rng,48,4); C=1e-3*randn(rng,4,48); D=zeros(4,4)
p=JC.RationalScatteringProvider(A,B,C,D)
ws=collect(range(0.,20.;length=benchsize(256))); out=zeros(ComplexF64,4,4,length(ws))
evaluate(x)=JC.evaluateprovider!(x...)
measure("rational48x4x$(benchsize(256))/evaluate",evaluate,(out,p,ws))
rf=JC.resolventfactors(A,B)
CZ=C*rf.Z
rw=JC.ResolventWorkspace(rf); rdest=zeros(ComplexF64,4,4)
transfer(x)=JC.transferat!(x[1],x[2],x[3],x[4],3.0,48,x[5])
measure("rational48x4/transfer",transfer,(rdest,rf,CZ,D,rw))

measure("ladder8/prepare",front,ladder(8))
