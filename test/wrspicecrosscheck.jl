using JosephsonCircuits
using Test
using XicTools_jll

# The transient against WRSPICE, where its executable is available: the
# same problem stepped by the package and run through the simulator, and
# the traces, the phases and the demodulation compared. The tolerances
# are those of WRSPICE's own accuracy at the chosen dphimax and of the
# small subgap loss its jj model keeps, which the package's junctions do
# not have.
@testset "the transient against WRSPICE" begin
    JC = JosephsonCircuits
    if !XicTools_jll.is_available()
        @info "skipping the WRSPICE cross check; XicTools_jll provides no wrspice for this platform"
    else
        w = 2pi*5e9
        drive = t -> 1e-7*sin(w*t)*(1 - exp(-t/1e-9))
        ts = (0.0, 5e-9)

        @testset "a driven JPA" begin
            circuit = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)),
                (:cc, 1, 2, Capacitor(100e-15)),
                (:jj, 2, 0, JosephsonJunction(1000e-12)),
                (:cj, 2, 0, Capacitor(1000e-15))])
            p = transientproblem(circuit; sources = [TransientSource(1, drive)])
            native = transientsolve(p, ts; dt = 1e-12, record = :phases)
            spice = transientsolve(p, ts; dt = 1e-12, record = :phases,
                method = WRspice())
            @test spice.times == native.times
            @test isnothing(spice.flux)
            @test all(isnan, spice.finalflux)
            # the port voltage trace and the outgoing wave
            @test maximum(abs, native.voltage .- spice.voltage) <
                0.02*maximum(abs, native.voltage)
            @test maximum(abs, native.outgoing .- spice.outgoing) <
                0.02*maximum(abs, native.outgoing)
            # the incident wave of a matched port is set by the drive
            # alone, so the two solvers agree to roundoff
            @test maximum(abs, native.incident .- spice.incident) <
                1e-12*maximum(abs, native.incident)
            # the junction phase from the phase node, in the package's
            # orientation: close to the package's record and far from
            # its negation
            @test size(spice.phases) == size(native.phases)
            @test maximum(abs, native.phases .- spice.phases) <
                0.05*maximum(abs, native.phases)
            @test maximum(abs, native.phases .+ spice.phases) >
                maximum(abs, native.phases)
            # the demodulated reflection agrees
            dn = transientdemodulate(native, 1, w/(2pi))
            ds = transientdemodulate(spice, 1, w/(2pi))
            @test abs(dn - ds) < 0.02*abs(dn)
            # the I/Q measurement reads the WRSPICE traces unchanged
            plan = transientiqplan(spice.times, [w/(2pi)]; duration = 1e-9)
            zn = transientiq(plan, native.outgoing)
            zs = transientiq(plan, spice.outgoing)
            @test maximum(abs, zn .- zs) < 0.02*maximum(abs, zn)
            # without a record of the phases none are kept
            @test isnothing(transientsolve(p, ts; dt = 1e-12,
                method = WRspice()).phases)
            # a decimated save is the same run printed less often
            sd = transientsolve(p, ts; dt = 1e-12, saveevery = 5,
                method = WRspice())
            @test sd.times == spice.times[1:5:end]
            @test maximum(abs, sd.voltage .- spice.voltage[:, 1:5:end]) <
                1e-4*maximum(abs, spice.voltage)
        end

        @testset "a mismatched transmission line in front of a junction" begin
            # reflections bounce on the 100 ps line at 30 ohms between
            # the 50 ohm port and the junction; the package steps the
            # line by the method of characteristics and WRSPICE by its
            # lossless line element
            circuit = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)),
                (:line, 1, 2, TransmissionLine(30.0, 0.1e-9*3e8; vp = 3e8)),
                (:cc, 2, 3, Capacitor(100e-15)),
                (:jj, 3, 0, JosephsonJunction(1000e-12)),
                (:cj, 3, 0, Capacitor(1000e-15))])
            p = transientproblem(circuit; sources = [TransientSource(1, drive)])
            native = transientsolve(p, ts; dt = 1e-12, method = GaussLegendre())
            spice = transientsolve(p, ts; dt = 1e-12, method = WRspice())
            @test maximum(abs, native.voltage .- spice.voltage) <
                0.03*maximum(abs, native.voltage)
            @test maximum(abs, native.outgoing .- spice.outgoing) <
                0.03*maximum(abs, native.outgoing)
            dn = transientdemodulate(native, 1, w/(2pi))
            ds = transientdemodulate(spice, 1, w/(2pi))
            @test abs(dn - ds) < 0.03*abs(dn)
        end

        @testset "junction orientations, a bias and a named source" begin
            # one junction written against the branch orientation, one
            # junction between two unearthed nodes, and a ramped current
            # through a named source
            circuit = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)),
                (:cc, 1, 2, Capacitor(100e-15)),
                (:jj1, 2, 3, JosephsonJunction(500e-12)),
                (:cj1, 2, 3, Capacitor(300e-15)),
                (:jj2, 3, 0, JosephsonJunction(500e-12)),
                (:cj2, 3, 0, Capacitor(1000e-15)),
                (:ib, 2, 0, CurrentSource(0.0))])
            p = transientproblem(circuit; sources = [TransientSource(1, drive),
                TransientSource(:ib, t -> 2e-8*(1 - exp(-t/1e-9)))])
            native = transientsolve(p, ts; dt = 1e-12, record = :phases)
            spice = transientsolve(p, ts; dt = 1e-12, record = :phases,
                method = WRspice())
            @test size(spice.phases, 1) == 2
            scale = maximum(abs, native.phases)
            @test maximum(abs, native.phases .- spice.phases) < 0.05*scale
            for r in axes(native.phases, 1)
                @test maximum(abs, native.phases[r, :] .+ spice.phases[r, :]) >
                    maximum(abs, native.phases[r, :])
            end
            @test maximum(abs, native.voltage .- spice.voltage) <
                0.02*maximum(abs, native.voltage)
        end
    end
end
