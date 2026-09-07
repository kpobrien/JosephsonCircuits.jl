# Photon normalized temporal modes of a real power wave trace: the
# measurement weights of a set of modes over the positive frequency bins of
# a record, built with FFTW once, then applied and pulled back as matrix
# products on the backend.

"""
    TransientQuantumPlan

A finite-record, positive-frequency temporal-mode measurement. Rows of its
output are `(X1,P1,X2,P2,...)`, with `[X,P]=im` and vacuum variance `1/2`.
`coefficients[:,j]` specifies mode j in the positive-frequency Fourier basis;
`gram` accounts for overlapping modes on the same port. `vacuum` and
`commutator` are the corresponding real covariance and commutator matrices.
Numeric measurement weights live on `backend`; metadata remain on the host.
Storage is O(samples*modes). Batch selected windows for long sliding records.
"""
struct TransientQuantumPlan{B,W}
    backend::B
    times::Vector{Float64}
    frequencies::Vector{Float64}
    ports::Vector{Int}
    coefficients::Matrix{ComplexF64}
    weights::W
    gram::Matrix{ComplexF64}
    vacuum::Matrix{Float64}
    commutator::Matrix{Float64}
    dt::Float64
end

function transientquantumgrid(times)
    ts = Float64.(collect(times))
    length(ts) >= 4 && all(isfinite, ts) ||
        throw(ArgumentError("At least four finite RF sample times are required."))
    # A late measurement window can lose several digits when subtracting
    # its first two RF times. Estimate the spacing over the full record.
    dt = (last(ts)-first(ts))/(length(ts)-1)
    spacingatol = 64eps(maximum(abs, ts))
    dt > 0 &&
    all(x -> x>0 && isapprox(x, dt; rtol = 1e-10, atol = spacingatol), diff(ts)) ||
        throw(ArgumentError("Quantum mode samples must be uniformly increasing."))
    return ts, dt, collect(1:fld(length(ts)-1, 2)) ./ (length(ts)*dt)
end

"""
    transientquantumplan(times, coefficients::AbstractMatrix; ports, backend=CPU())
    transientquantumplan(times, frequencies::AbstractVector; ports=fill(1,...),
        envelopes=nothing, backend=CPU())

Define photon-normalized temporal modes on a uniformly sampled half-open record
`[times[1], times[1]+length(times)*dt)`. Do not include the repeated right endpoint.
For N samples the coefficient rows are Fourier bins `k=1:fld(N-1,2)` at
`f=k/(N*dt)`; DC and a self-conjugate Nyquist bin are excluded. Each coefficient
column must have unit Euclidean norm. Different columns may overlap.

The frequency convenience form creates bin-aligned monochromatic modes when
`envelopes=nothing`. Otherwise each column of `envelopes[sample,mode]` defines
`g(t)=envelope(t)*exp(-2pi*im*f*(t-times[1]))`; project it onto positive Fourier
bins and normalize the resulting coefficients. Supply enough RF bandwidth to
resolve the carrier and envelope. Frequencies must lie strictly below Nyquist.

A canonical bin has real physical power wave
`w(t)=sqrt(h*f/(N*dt))*(X*cos(2pi*f*t)+P*sin(2pi*f*t))`.
The readout includes the frequency-dependent `1/sqrt(h*f)` weighting before
combining bins. A resolved cosine of peak amplitude A in one full-record bin
has X=A*sqrt(N*dt/(h*f)), P=0 and mean photon number X^2/2.
The P convention is the negative of `imag(transientiq(...))` for that cosine's
classical phasor. Existing peak-amplitude I/Q conventions remain unchanged.

Windows are finite-record mode definitions, not independent white-noise samples.
Use `gram`, `vacuum` and `commutator` when modes/windows overlap. The plan is
read-only after construction and has no mutable shared FFT workspace.
"""
function transientquantumplan(times, coefficients::AbstractMatrix;
        ports = fill(1, size(coefficients, 2)), backend::Backend = CPU())
    ts, dt, fs = transientquantumgrid(times)
    c = ComplexF64.(Array(coefficients))
    m = size(c, 2)
    size(c, 1) == length(fs) && m > 0 ||
        throw(DimensionMismatch("Incorrect positive-frequency coefficient shape."))
    length(ports) == m && all(p->p isa Integer && p>0, ports) ||
        throw(ArgumentError("Supply one positive compiled port index per mode."))
    all(isfinite, c) && all(j->isapprox(sum(abs2, view(c, :, j)), 1; rtol = 1e-10), 1:m) ||
        throw(ArgumentError("Each temporal-mode column must be finite and have unit norm."))
    n = length(ts)
    weights = zeros(n, 2m)
    spectrum = zeros(ComplexF64, n)
    for j in 1:m
        fill!(spectrum, 0)
        spectrum[2:(length(fs)+1)] .= conj.(c[:, j]) ./ sqrt.(planck_constant .* fs)
        # Unnormalized inverse FFT has the positive exponent extracting a_k.
        w = FFTW.bfft(spectrum) .* (2dt/sqrt(n*dt))
        weights[:, 2j-1] .= real.(w)
        weights[:, 2j] .= imag.(w)
    end
    gram = c'c
    vacuum, comm = zeros(2m, 2m), zeros(2m, 2m)
    for j in 1:m, k in 1:m

        ports[j] == ports[k] || (gram[j, k] = 0)
        r, s = real(gram[j, k]), imag(gram[j, k])
        vacuum[2j-1, 2k-1] = vacuum[2j, 2k] = r/2
        vacuum[2j-1, 2k] = -s/2
        vacuum[2j, 2k-1] = s/2
        comm[2j-1, 2k-1] = comm[2j, 2k] = s
        comm[2j-1, 2k] = r
        comm[2j, 2k-1] = -r
    end
    return TransientQuantumPlan(backend, ts, fs, Int.(ports), c,
        tobackend(backend, weights), gram, vacuum, comm, dt)
end

function transientquantumplan(times, frequencies::AbstractVector;
        ports = fill(1, length(frequencies)), envelopes = nothing, backend::Backend = CPU())
    ts, dt, fs = transientquantumgrid(times)
    f = Float64.(frequencies)
    all(x->isfinite(x) && 0<x<0.5/dt, f) ||
        throw(ArgumentError("Mode frequencies must be positive and below Nyquist."))
    n, m = length(ts), length(f)
    c = zeros(ComplexF64, length(fs), m)
    if isnothing(envelopes)
        for j in 1:m
            k = round(Int, f[j]*n*dt)
            1 <= k <= length(fs) && isapprox(f[j], fs[k]; rtol = 1e-10) ||
                throw(ArgumentError("Full-record tones must lie on Fourier bins; supply envelopes for other temporal modes."))
            c[k, j] = 1
        end
    else
        size(envelopes) == (n, m) && all(isfinite, envelopes) ||
            throw(DimensionMismatch("Envelopes must be a finite samples-by-modes matrix."))
        for j in 1:m
            g = envelopes[:, j] .* cispi.(-2f[j] .* (ts .- first(ts)))
            c[:, j] .= FFTW.bfft(g)[2:(length(fs)+1)]
            norm(view(c, :, j)) > 0 ||
                throw(ArgumentError("Temporal mode has no positive-frequency support."))
            c[:, j] ./= norm(view(c, :, j))
        end
    end
    return transientquantumplan(ts, c; ports, backend)
end

function transientquantumcheck(plan, traces, out)
    traces isa AbstractMatrix{<:Real} && out isa AbstractVector{<:Real} ||
        throw(ArgumentError("Quantum readouts require real RF traces and real quadrature vectors."))
    size(traces, 2) == length(plan.times) && size(traces, 1) >= maximum(plan.ports) &&
    length(out) == 2length(plan.ports) ||
        throw(DimensionMismatch("Quantum measurement dimensions do not match the plan."))
    b = KernelAbstractions.get_backend(plan.weights)
    KernelAbstractions.get_backend(traces) == b &&
    KernelAbstractions.get_backend(out) == b ||
        throw(ArgumentError("Measurement weights, RF traces and quadratures must share a backend."))
    return nothing
end

"""
    transientquantum!(out, plan, traces)

The photon normalized quadratures `(X1, P1, X2, P2, ...)` of the temporal
modes of `plan` measured on the real power wave traces `traces` in
sqrt(W), a `(port, time)` matrix on the plan's backend covering the
plan's record, written into the real vector `out`: one product of the
plan's weights with each mode's port trace, so a mode's quadratures
are exact linear functionals of the trace. The traces, the weights and
the output must share a backend.
"""
function transientquantum!(out, plan::TransientQuantumPlan, traces)
    transientquantumcheck(plan, traces, out)
    for j in eachindex(plan.ports)
        q = (2j-1):2j
        mul!(view(out, q), transpose(view(plan.weights, :, q)), view(traces, plan.ports[j], :))
    end
    return out
end
"""
    transientquantum(plan, traces)

The quadratures of the temporal modes of `plan` measured on `traces`,
allocated on the plan's backend; see [`transientquantum!`](@ref).
"""
function transientquantum(plan::TransientQuantumPlan, traces)
    out = KernelAbstractions.zeros(plan.backend, Float64, 2length(plan.ports))
    return transientquantum!(out, plan, traces)
end

"""
    transientquantumvjp!(out, plan, weights)

The exact transpose of [`transientquantum!`](@ref): the real trace
weights, a `(port, time)` matrix on the plan's backend, of the
functional `dot(weights, transientquantum(plan, traces))` for the
quadrature `weights` given, written into `out`. Placed on the port rows
of the outgoing waves over the plan's window, these are the `weights` of
[`transientadjoint`](@ref), which propagates a measured quadrature back
to the currents and the initial state; the noise builds its objectives
this way.
"""
function transientquantumvjp!(out, plan::TransientQuantumPlan, weights)
    transientquantumcheck(plan, out, weights)
    fill!(out, 0)
    for j in eachindex(plan.ports)
        q = (2j-1):2j
        mul!(view(out, plan.ports[j], :),
            view(plan.weights, :, q), view(weights, q), 1.0, 1.0)
    end
    return out
end

"""
    transientquantumdiagnostics(covariance, commutator, expected; rtol=1e-3)

Check commutator closure and the Gaussian uncertainty inequality
`covariance + im*expected/2 >= 0`. Small output matrices are inspected on the
host. Passing is necessary, not sufficient: refine RF sampling, bath cutoff,
frequency quadrature and solver tolerances independently.
"""
function transientquantumdiagnostics(covariance, commutator, expected; rtol = 1e-3)
    v, k, e = Matrix{Float64}(Array(covariance)), Matrix{Float64}(Array(commutator)),
    Matrix{Float64}(Array(expected))
    size(v) == size(k) == size(e) && size(v, 1) == size(v, 2) && !isempty(v) ||
        throw(DimensionMismatch("Quantum covariance/commutator shapes must agree."))
    isfinite(rtol) && rtol > 0 && all(isfinite, v) && all(isfinite, k) &&
    all(isfinite, e) ||
        throw(ArgumentError("Invalid quantum diagnostics data or tolerance."))
    isapprox(v, v'; rtol = 1e-10, atol = 1e-12) &&
    isapprox(e, -e'; rtol = 1e-10, atol = 1e-12) ||
        throw(ArgumentError("Expected a symmetric covariance and antisymmetric target commutator."))
    error = norm(k-e)/max(norm(e), eps())
    uncertainty = eigmin(Hermitian((v+v')/2+im*e/2))
    return (; commutationerror = error, minimumuncertaintyeigenvalue = uncertainty,
        passed = error<=rtol && uncertainty>=-rtol*max(opnorm(v), 1.0))
end

"""
    transientquantumefficiency(gain, covariance; rtol=1e-3)

Phase-preserving single-mode metrics for a real 2x2 quadrature gain and an
isotropic canonical output covariance (vacuum = I/2). `gain` is the incremental
response to a coherent displacement, not a ratio of large loaded amplitudes.
Returns photon gain, input-referred added noise, QE, QEideal and QE/QEideal.
Uses the existing HB ideal-efficiency convention. Phase-sensitive gain or
anisotropic noise is rejected: retain the full gain/covariance for those cases.
The caller must also verify commutator closure and bath-basis convergence.
"""
function transientquantumefficiency(gain, covariance; rtol = 1e-3)
    a, v = Matrix{Float64}(Array(gain)), Matrix{Float64}(Array(covariance))
    size(a) == size(v) == (2, 2) ||
        throw(DimensionMismatch("Single-mode gain and covariance must be 2x2."))
    all(isfinite, a) && all(isfinite, v) && isfinite(rtol) && rtol>0 ||
        throw(ArgumentError("Invalid quantum efficiency data."))
    g, s = sum(abs2, a)/2, tr(v)/2
    g>0 && s>0 && det(a)>0 && isapprox(a*a', g*I; rtol) && isapprox(v, s*I; rtol) ||
        throw(ArgumentError("Scalar QE requires phase-preserving gain and isotropic noise."))
    ideal = only(calcqeideal(reshape([sqrt(g)], 1, 1)))
    qe = g/(2s)
    return (;
        gain = g, addednoise = s/g-0.5, QE = qe, QEideal = ideal, normalizedQE = qe/ideal)
end
