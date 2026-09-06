# Circuits shared across test files, so a change to a canonical test
# device is made in one place. Files include this defensively
# (`isdefined(Main, ...) || include(...)`) so each test file still runs
# standalone.

# a JPA: one port, one junction, a coupling capacitor -- symbolic values
# with a definitions dictionary, the standard user-facing form
function testjpacircuit()
    circuit = Tuple{String,String,String,Union{Complex{Float64},Symbol,Int64}}[]
    push!(circuit,("P1","1","0",1))
    push!(circuit,("R1","1","0",:Rleft))
    push!(circuit,("C1","1","2",:Cc))
    push!(circuit,("Lj1","2","0",:Lj))
    push!(circuit,("C2","2","0",:Cj))
    circuitdefs = Dict(
        :Lj => 1000.0e-12,
        :Cc => 100.0e-15,
        :Cj => 1000.0e-15,
        :Rleft => 50.0,
    )
    return circuit, circuitdefs
end

# the same JPA with numeric literal values and an empty definitions
# dictionary, which exercises the numeric-value parsing path
function testjpacircuitnumeric()
    circuit = [("P1","1","0",1), ("R1","1","0",50.0), ("C1","1","2",100e-15),
               ("Lj1","2","0",1000e-12), ("C2","2","0",1000e-15)]
    return circuit, Dict{Any,Any}()
end

# a four junction transmission line chain with a port at each end
function testchaincircuit()
    circuit = Tuple{String,String,String,Union{Complex{Float64},Symbol,Int64}}[]
    push!(circuit, ("P1","1","0",1)); push!(circuit, ("R1","1","0",:R))
    for i in 1:4
        push!(circuit, ("Lj$(i)","$(i)","$(i+1)",:Lj))
        push!(circuit, ("C$(i)","$(i)","0",:Cg))
    end
    push!(circuit, ("C5","5","0",:Cg)); push!(circuit, ("R2","5","0",:R))
    circuitdefs = Dict{Symbol,Complex{Float64}}(
        :Lj => 100e-12, :Cg => 40e-15, :R => 50.0)
    return circuit, circuitdefs
end

# --- the cross check matrix -------------------------------------------------
#
# The circuits the systematic comparisons of test/crosscheck.jl run over: one
# of each class the package handles, each small enough that the whole matrix
# costs about as much as one large test file used to. A case records the
# circuit, its pump and sources, the signal frequencies, the harmonic counts
# and the solver options it needs, and what is true of it (whether it has a
# pump solve, whether it is lossless, whether it is linear), so that the
# comparisons can skip what does not apply. The weak signal `Is` of the two
# tone comparison is 1e-7 of the largest source, so that the nonlinear
# response to it is linear to that order.
function crosscheckcases()
    GHz = 2*pi*1e9
    Z0 = 50.0
    cases = Any[]

    # a lossless linear two port: a resonator between two coupling
    # capacitors. No pump solve; the linearized sweep alone.
    push!(cases, (name = "linear lossless resonator",
        circuit = Circuit([(:p1, 1, 0, Port(1; Z0 = Z0)),
            (:cc1, 1, 2, Capacitor(100e-15)), (:l1, 2, 0, Inductor(1e-9)),
            (:c1, 2, 0, Capacitor(1e-12)), (:cc2, 2, 3, Capacitor(100e-15)),
            (:p2, 3, 0, Port(2; Z0 = Z0))]),
        defs = Dict{Any,Any}(), ws = GHz .* [4.0, 5.0, 6.0], wp = (GHz*5.0,),
        sources = [], signalport = 1, Is = 0.0,
        Npump = (1,), Nmod = (1,), kw = (;),
        pumped = false, lossless = true, linear = true))

    # the same with a resistor inside: a passive lossy two port, whose noise
    # covariance Bosma's theorem gives in closed form
    push!(cases, (name = "linear lossy resonator",
        circuit = Circuit([(:p1, 1, 0, Port(1; Z0 = Z0)),
            (:cc1, 1, 2, Capacitor(100e-15)), (:l1, 2, 0, Inductor(1e-9)),
            (:c1, 2, 0, Capacitor(1e-12)), (:r1, 2, 0, Resistor(2000.0)),
            (:cc2, 2, 3, Capacitor(100e-15)), (:p2, 3, 0, Port(2; Z0 = Z0))]),
        defs = Dict{Any,Any}(), ws = GHz .* [4.0, 5.0, 6.0], wp = (GHz*5.0,),
        sources = [], signalport = 1, Is = 0.0,
        Npump = (1,), Nmod = (1,), kw = (;),
        pumped = false, lossless = false, linear = true))

    # the four wave mixing JPA: one port, one junction, a coupling capacitor
    jpa(extra...) = Circuit([(:p1, 1, 0, Port(1; Z0 = Z0)),
        (:cc, 1, 2, Capacitor(100e-15)), (:jj, 2, 0, JosephsonJunction(1000e-12)),
        (:cj, 2, 0, Capacitor(1000e-15)), extra...])
    Ip = 0.00565e-6
    push!(cases, (name = "4WM JPA",
        circuit = jpa(), defs = Dict{Any,Any}(),
        ws = GHz .* [4.6, 4.75, 4.9], wp = (GHz*4.75001,),
        sources = [(mode = (1,), port = 1, current = Ip)], signalport = 1,
        Is = 1e-7*Ip, Npump = (4,), Nmod = (4,), kw = (;),
        pumped = true, lossless = true, linear = false))
    push!(cases, (name = "4WM JPA with loss",
        circuit = jpa((:rl, 2, 0, Resistor(2.0e5))), defs = Dict{Any,Any}(),
        ws = GHz .* [4.6, 4.75, 4.9], wp = (GHz*4.75001,),
        sources = [(mode = (1,), port = 1, current = Ip)], signalport = 1,
        Is = 1e-7*Ip, Npump = (4,), Nmod = (4,), kw = (;),
        pumped = true, lossless = false, linear = false))

    # the four cell junction chain with a port at one end and a resistive
    # termination, a noise channel, at the other
    chain, chaindefs = testchaincircuit()
    push!(cases, (name = "4WM junction chain",
        circuit = chain, defs = chaindefs,
        ws = GHz .* [4.03, 5.47, 7.11], wp = (GHz*6.0,),
        sources = [(mode = (1,), port = 1, current = 0.6e-6)], signalport = 1,
        Is = 1e-7*0.6e-6, Npump = (4,), Nmod = (4,), kw = (;),
        pumped = true, lossless = false, linear = false))

    # three wave mixing through a direct current bias: a junction resonator
    # fed through an inductor, so that the direct current of the port bias
    # reaches the junction (a coupling capacitor would block it) and sets
    # the odd order nonlinearity; pumped through the port at twice the
    # resonance
    # the junction's critical current is 0.66 uA, so the bias of 0.3 uA
    # is at 0.45 of it and the pump well below
    Lj = 500e-12; Lb = 300e-12; Cr = 1.0e-12
    wr = 1/sqrt((Lb + Lj)*Cr)
    push!(cases, (name = "3WM with a DC bias",
        circuit = Circuit([(:p1, 1, 0, Port(1; Z0 = Z0)),
            (:lb, 1, 2, Inductor(Lb)), (:jj, 2, 0, JosephsonJunction(Lj)),
            (:cr, 2, 0, Capacitor(Cr))]),
        defs = Dict{Any,Any}(),
        ws = wr .* [0.9713, 1.0000, 1.0291], wp = (2*wr*1.0017,),
        sources = [(mode = (0,), port = 1, current = 0.3e-6),
            (mode = (1,), port = 1, current = 0.1e-6)], signalport = 1,
        Is = 1e-7*0.3e-6, Npump = (4,), Nmod = (4,),
        kw = (dc = true, threewavemixing = true, fourwavemixing = true),
        pumped = true, lossless = true, linear = false))

    # a flux biased SQUID: two junctions in a loop closed by an inductor,
    # the loop inductor coupled to a flux line which a high impedance port
    # drives with a direct current; the signal port pumps at twice the
    # resonance (biased), or the flux line does (flux pumped)
    squid(pumpport) = (
        circuit = Circuit([(:p1, 1, 0, Port(1; Z0 = Z0)),
            (:cc, 1, 2, Capacitor(60e-15)),
            (:jj1, 2, 0, JosephsonJunction(400e-12)), (:cj1, 2, 0, Capacitor(0.6e-12)),
            (:ll, 2, 3, Inductor(40e-12)),
            (:jj2, 3, 0, JosephsonJunction(400e-12)), (:cj2, 3, 0, Capacitor(0.6e-12)),
            (:ldc, 4, 0, Inductor(200e-12)), (:k1, :ll, :ldc, MutualInductor(0.9)),
            (:p2, 4, 0, Port(2; Z0 = 1000.0))]),
        pumpport = pumpport)
    wsq = 1/sqrt(200e-12*1.2e-12)
    # the mutual inductance is 0.9*sqrt(40*200) pH = 80 pH, so 8 uA on
    # the flux line threads 0.3 flux quanta: a single stable state
    for (name, pumpport, Ip) in (("flux biased SQUID", 1, 0.3e-6),
            ("flux pumped SQUID", 2, 4.0e-6))
        s = squid(pumpport)
        push!(cases, (name = name, circuit = s.circuit, defs = Dict{Any,Any}(),
            ws = wsq .* [0.9713, 1.0000, 1.0291], wp = (2*wsq*1.0017,),
            sources = [(mode = (0,), port = 2, current = 8.0e-6),
                (mode = (1,), port = pumpport, current = Ip)], signalport = 1,
            Is = 1e-7*8.0e-6, Npump = (4,), Nmod = (4,),
            kw = (dc = true, threewavemixing = true, fourwavemixing = true),
            pumped = true, lossless = true, linear = false))
    end

    # a JPA whose junction is fed through a transformer: mutual inductors
    # promote coupled branches to auxiliary unknowns
    push!(cases, (name = "mutual inductor JPA",
        circuit = Circuit([(:p1, 1, 0, Port(1; Z0 = Z0)),
            (:cc, 1, 2, Capacitor(100e-15)), (:l1, 2, 0, Inductor(300e-12)),
            (:l2, 3, 0, Inductor(300e-12)), (:k1, :l1, :l2, MutualInductor(0.9)),
            (:jj, 3, 0, JosephsonJunction(500e-12)), (:cj, 3, 0, Capacitor(1000e-15))]),
        defs = Dict{Any,Any}(),
        ws = GHz .* [4.5, 5.0, 5.5], wp = (GHz*5.0001,),
        sources = [(mode = (1,), port = 1, current = 1.5e-6)], signalport = 1,
        Is = 1e-7*1.5e-6, Npump = (4,), Nmod = (4,), kw = (;),
        pumped = true, lossless = true, linear = false))

    # the JPA behind a dissipative scattering block, a 3 dB attenuator
    att = ScatteringParameters(JosephsonCircuits.ABCDtoS(
        JosephsonCircuits.ABCD_attenuator_T(Z0, 3.0)); nports = 2, zref = Z0)
    push!(cases, (name = "scattering block JPA",
        circuit = Circuit([(:p1, 1, 0, Port(1; Z0 = Z0)), (:att, 1, 2, att),
            (:cc, 2, 3, Capacitor(100e-15)), (:jj, 3, 0, JosephsonJunction(1000e-12)),
            (:cj, 3, 0, Capacitor(1000e-15))]),
        defs = Dict{Any,Any}(),
        ws = GHz .* [4.6, 4.75, 4.9], wp = (GHz*4.75001,),
        sources = [(mode = (1,), port = 1, current = 1.4*Ip)], signalport = 1,
        Is = 1e-7*Ip, Npump = (4,), Nmod = (4,), kw = (;),
        pumped = true, lossless = false, linear = false))

    # the JPA fed through a lossless transmission line
    push!(cases, (name = "transmission line JPA",
        circuit = Circuit([(:p1, 1, 0, Port(1; Z0 = Z0)),
            (:tl, 1, 2, TransmissionLine(Z0, 5e-3; noise = Lossless())),
            (:cc, 2, 3, Capacitor(100e-15)), (:jj, 3, 0, JosephsonJunction(1000e-12)),
            (:cj, 3, 0, Capacitor(1000e-15))]),
        defs = Dict{Any,Any}(),
        ws = GHz .* [4.6, 4.75, 4.9], wp = (GHz*4.75001,),
        sources = [(mode = (1,), port = 1, current = Ip)], signalport = 1,
        Is = 1e-7*Ip, Npump = (4,), Nmod = (4,), kw = (;),
        pumped = true, lossless = true, linear = false))

    return cases
end

# the keyword names of `hbnlsolve` for the mode selection `hbsolve` spells
# as `dc`, `threewavemixing` and `fourwavemixing`
function nonlinearkw(kw::NamedTuple)
    return (dc = get(kw, :dc, false), even = get(kw, :threewavemixing, false),
        odd = get(kw, :fourwavemixing, true))
end
