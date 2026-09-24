# Sensitivity benchmark. Run in a fresh process with one Julia/BLAS thread:
# julia --startup-file=no --project=. benchmark/sensitivities.jl label
# The pump solve is outside the timers. Each derivative of the operating
# point is set up per call, so a call with a few components measures that
# setup and a call with every component measures the derivatives
# themselves.
using JosephsonCircuits, LinearAlgebra, SparseArrays, Statistics
include(joinpath(@__DIR__, "options.jl"))
const JC = JosephsonCircuits
const label = BENCH_LABEL
const NCHAIN = benchsize(128)
function measure(case, f, input)
    @nospecialize f input
    firstcall, compile, samples = benchtime(f, input)
    println(join((label,case,firstcall.time,compile,
        median(s.time for s in samples),median(s.bytes for s in samples)), '\t'))
    flush(stdout)
    return firstcall.value
end
# a chain of junctions with a port at each end
function chain(n)
    e = Any[("P1","1","0",Port(1;Z0=:R))]
    for i in 1:n
        push!(e,("Lj$i","$i","$(i+1)",JosephsonJunction(:Lj)),
            ("C$i","$i","0",Capacitor(:Cg)))
    end
    push!(e,("C$(n+1)","$(n+1)","0",Capacitor(:Cg)),("R2","$(n+1)","0",Resistor(:R)))
    return compile(Circuit(e))
end
println("label\tcase\tfirst_s\tcompile_s\twarm_s\tbytes")
c = chain(NCHAIN)
defs = Dict(:Lj=>100e-12,:Cg=>40e-15,:R=>50.0)
nl = hbnlsolve((2*pi*7e9,),(16,),[(mode=(1,),port=1,current=1e-8)],c,defs;
    returnoperatingpoint=true,keyedarrays=false)
op = nl.operatingpoint
nm = numericmatrices(c,defs;Nmodes=op.Nmodes)
few = ["Lj1","C1"]
every = [n for n in c.componentnames if startswith(n,"Lj") || startswith(n,"C")]
indices(names) = [JC.componentindex(c,n) for n in names]
residual(x) = JC.calcresidualsensitivity(op,c,nm,x)
dfew = measure("chain$(NCHAIN)/residual, 2 components",residual,indices(few))
devery = measure("chain$(NCHAIN)/residual, $(length(every)) components",residual,indices(every))
ws = 2*pi*collect(range(4e9,6e9;length=4))
sweep(x) = hblinsolve(ws,c,defs;Nmodulationharmonics=(2,),nonlinear=nl,
    nbatches=1,keyedarrays=false,sensitivitynames=x.names,
    returnSsensitivity=true,sensitivityresidual=x.dFr,
    sensitivitymode=x.mode).Ssensitivity
measure("chain$(NCHAIN)/forward sweep, 2 components",sweep,
    (names=few,dFr=dfew,mode=:forward))
measure("chain$(NCHAIN)/reverse sweep, $(length(every)) components",sweep,
    (names=every,dFr=devery,mode=:reverse))
