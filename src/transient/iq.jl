# A windowed I/Q measurement of a real port trace: a causal sliding window
# demodulated at each carrier by zero padded FFT convolution, on the
# package's FFTW and cuFFT plans and its backend convention, and the exact
# transpose of that measurement for an adjoint objective.

"""
    TransientIQPlan

Reusable causal I/Q measurement plan. `times` are the right edges of complete
windows; `centertimes` subtract the filter's `groupdelay`. Frequencies and
`bandwidth3db` are in Hz, times in seconds. `noisebandwidth` is
one-sided, `sum(abs2, taps)/(2dt)`, for the unity-DC-gain low-pass filter.

The FFT workspace is mutable. Use separate plans for concurrent measurements.
"""
struct TransientIQPlan{B,V,F,R}
    backend::B
    frequencies::Vector{Float64}
    ports::Vector{Int}
    times::Vector{Float64}
    centertimes::Vector{Float64}
    groupdelay::Float64
    bandwidth3db::Float64
    noisebandwidth::Float64
    dt::Float64
    start::Float64
    phasereference::Float64
    nsamples::Int
    ntaps::Int
    stride::Int
    taps::Vector{Float64}
    work::V
    spectrum::V
    forward::F
    backward::R
end

function transientiqfftplans(work, ::CPU)
    return FFTW.plan_fft!(work; flags = FFTW.ESTIMATE),
    FFTW.plan_bfft!(work; flags = FFTW.ESTIMATE)
end

transientiqfftplans(work, backend) = throw(ArgumentError(
    "I/Q FFT plans are not implemented for $(typeof(backend)). Load CUDA for the CUDA backend."))

function transientiqbandwidth(taps, dt)
    # Symmetric taps have a real centered frequency response. The first
    # half-power crossing is in this interval for both supported windows.
    center = (length(taps)-1)/2
    magnitude(f) = abs(sum(taps[k]*cospi(2f*dt*(k-1-center)) for k in eachindex(taps)))
    lo, hi = 0.0, min(0.5/dt, 2/(length(taps)*dt))
    for _ in 1:60
        mid = (lo+hi)/2
        if magnitude(mid) > inv(sqrt(2))
            lo = mid
        else
            hi = mid
        end
    end
    return (lo+hi)/2
end

"""
    transientiqplan(times, frequencies; duration, window=:hann,
        ports=fill(1,length(frequencies)), stride=1, phasereference=first(times),
        backend=CPU())

Plan causal sliding I/Q measurements of uniformly sampled **real** port traces
with shape `(port, time)`. `ports` selects a trace row for each carrier in Hz.
For normalized window taps `h[k]`, the output is

`z[c,n] = 2 sum(h[k] x[ports[c],n-k] exp(-2pi*im*f[c]*(t[n-k]-phasereference)))`.

Thus a resolved cosine of peak amplitude `A` and phase `phi` gives approximately
`A*exp(im*phi)` when the doubled-carrier image is rejected by the window. The
units match the input (volts or sqrt(watts)); these are peak, not RMS or photon
amplitudes. I is `real(z)` and Q is `imag(z)`.

`window` is `:hann` or `:rectangular`; `duration` is rounded to an integer
number of sample intervals (at least three). Only complete windows are kept,
every `stride` samples, without assumed prehistory. Use `plan.times` for causal
availability or `plan.centertimes` for delay-corrected plotting. The reported
3 dB bandwidth is the positive-frequency half-width of the low-pass filter.
Neither window is a brick-wall filter: choose sampling, duration and stride
to control image leakage and aliasing. The input must already resolve the RF.

CPU uses FFTW and CUDA uses cuFFT with KernelAbstractions kernels. One FFT
workspace is reused across channels; this avoids a samples-by-window array.
"""
function transientiqplan(times, frequencies; duration, window = :hann,
        ports = fill(1, length(frequencies)), stride = 1,
        phasereference = first(times), backend = CPU())
    ts, fs = Float64.(collect(times)), Float64.(collect(frequencies))
    length(ts) >= 4 && all(isfinite, ts) ||
        throw(ArgumentError("At least four finite sample times are required."))
    dt = ts[2]-ts[1]
    dt > 0 && all(d -> isapprox(d, dt; rtol = 1e-10,
            atol = 64eps(maximum(abs, ts))), diff(ts)) ||
        throw(ArgumentError("Sample times must be uniformly increasing."))
    !isempty(fs) && all(f -> isfinite(f) && 0 < f < 0.5/dt, fs) ||
        throw(ArgumentError("Carrier frequencies must lie strictly between zero and Nyquist."))
    length(ports) == length(fs) && all(p -> p isa Integer && p >= 1, ports) ||
        throw(ArgumentError("Provide one positive integer trace row per carrier."))
    stride isa Integer && stride >= 1 ||
        throw(ArgumentError("stride must be a positive integer."))
    isfinite(duration) && 3dt <= duration <= (length(ts)-1)*dt*(1+1e-10) ||
        throw(ArgumentError("duration must span at least three intervals and fit in the trace."))
    isfinite(phasereference) || throw(ArgumentError("phasereference must be finite."))
    ntaps = round(Int, duration/dt)+1
    4 <= ntaps <= length(ts) ||
        throw(ArgumentError("Rounded window does not fit in the trace."))
    taps = if window === :hann
        [0.5-0.5cospi(2k/(ntaps-1)) for k in 0:(ntaps-1)]
    elseif window === :rectangular
        ones(ntaps)
    else
        throw(ArgumentError("window must be :hann or :rectangular."))
    end
    taps ./= sum(taps)
    nfft = nextpow(2, length(ts)+ntaps-1)
    work = KernelAbstractions.zeros(backend, ComplexF64, nfft)
    padded = zeros(ComplexF64, nfft)
    padded[1:ntaps] .= taps
    spectrum = tobackend(backend, padded)
    forward, backward = transientiqfftplans(work, backend)
    mul!(spectrum, forward, spectrum)
    KernelAbstractions.synchronize(backend)
    outputtimes = ts[ntaps:stride:end]
    delay = (ntaps-1)*dt/2
    return TransientIQPlan(
        backend, fs, Int.(collect(ports)), outputtimes, outputtimes .- delay,
        delay, transientiqbandwidth(taps, dt), sum(abs2, taps)/(2dt), dt, first(ts),
        Float64(phasereference), length(ts), ntaps, Int(stride), taps, work, spectrum, forward, backward)
end

@kernel function iqmixkernel!(work, @Const(traces), row, n, t0, dt, f)
    i = @index(Global)
    @inbounds work[i] = i <= n ? traces[row, i]*cispi(-2f*(t0+(i-1)*dt)) :
                        zero(eltype(work))
end

@kernel function iqextractkernel!(out, @Const(work), row, ntaps, stride, scale)
    j = @index(Global)
    @inbounds out[row, j] = scale*work[ntaps+(j-1)*stride]
end

@kernel function iqseedkernel!(work, @Const(weights), row, n, ntaps, stride)
    i = @index(Global)
    k = i-ntaps
    @inbounds work[i] = 0 <= k && i <= n && k % stride == 0 ?
                        weights[row, 1+k÷stride] : zero(eltype(work))
end

@kernel function iqpullbackkernel!(out, @Const(work), row, t0, dt, f, scale)
    i = @index(Global)
    @inbounds out[row, i] += scale*real(work[i]*cispi(2f*(t0+(i-1)*dt)))
end

function transientiqcheck(plan, rf, iq)
    rf isa AbstractMatrix{<:Real} && iq isa AbstractMatrix{<:Complex} ||
        throw(ArgumentError("RF traces must be a real matrix and I/Q a complex matrix."))
    size(rf, 2) == plan.nsamples && size(rf, 1) >= maximum(plan.ports) &&
    size(iq) == (length(plan.frequencies), length(plan.times)) ||
        throw(DimensionMismatch("Trace or I/Q dimensions do not match the measurement plan."))
    # Compare storage backends, not launch options such as CPU(static=true).
    storagebackend = KernelAbstractions.get_backend(plan.work)
    KernelAbstractions.get_backend(rf) == storagebackend &&
    KernelAbstractions.get_backend(iq) == storagebackend ||
        throw(ArgumentError("Trace, I/Q and plan must use the same backend."))
    return nothing
end

"""
    transientiq!(out, plan, traces)

The complex I/Q measurements of `plan` on the real port traces `traces`,
a `(port, time)` matrix on the plan's backend, written into `out`, a
complex `(channel, window)` matrix: one row per carrier of the plan and
one column per complete window, in the peak amplitude convention of
[`transientiqplan`](@ref). Each channel is one mixing kernel, one
forward and one backward FFT on the plan's workspace and one extraction
kernel, so no samples by windows array is formed.
"""
function transientiq!(out, plan::TransientIQPlan, traces)
    transientiqcheck(plan, traces, out)
    nfft = length(plan.work)
    for c in eachindex(plan.frequencies)
        iqmixkernel!(plan.backend, 256)(plan.work, traces, plan.ports[c],
            plan.nsamples, plan.start-plan.phasereference, plan.dt, plan.frequencies[c]; ndrange = nfft)
        KernelAbstractions.synchronize(plan.backend)
        mul!(plan.work, plan.forward, plan.work)
        plan.work .*= plan.spectrum
        mul!(plan.work, plan.backward, plan.work)
        iqextractkernel!(plan.backend, 256)(out, plan.work, c, plan.ntaps,
            plan.stride, 2/nfft; ndrange = length(plan.times))
        KernelAbstractions.synchronize(plan.backend)
    end
    return out
end

"""
    transientiq(plan, traces)

The complex I/Q measurements of `plan` on the real port traces `traces`,
allocated on the plan's backend; see [`transientiq!`](@ref).
"""
function transientiq(plan::TransientIQPlan, traces)
    out = KernelAbstractions.zeros(plan.backend, ComplexF64,
        length(plan.frequencies), length(plan.times))
    return transientiq!(out, plan, traces)
end

"""
    transientiqvjp!(rfweights, plan, iqweights)

Overwrite the real RF gradient of `real(dot(iqweights, transientiq(plan,rf)))`.
This is the exact transpose of the finite-window measurement, including
normalization, overlap, decimation and multiple channels reading the same port.
These RF weights, placed on the port rows of the outgoing waves, are the
`weights` of [`transientadjoint`](@ref), which propagates them to the
currents and the initial state. The routine does not itself integrate an
adjoint.
"""
function transientiqvjp!(out, plan::TransientIQPlan, weights)
    transientiqcheck(plan, out, weights)
    fill!(out, 0)
    nfft = length(plan.work)
    for c in eachindex(plan.frequencies)
        iqseedkernel!(plan.backend, 256)(plan.work, weights, c, plan.nsamples,
            plan.ntaps, plan.stride; ndrange = nfft)
        KernelAbstractions.synchronize(plan.backend)
        mul!(plan.work, plan.forward, plan.work)
        plan.work .*= conj.(plan.spectrum)
        mul!(plan.work, plan.backward, plan.work)
        iqpullbackkernel!(plan.backend, 256)(out, plan.work, plan.ports[c],
            plan.start-plan.phasereference, plan.dt, plan.frequencies[c], 2/nfft; ndrange = plan.nsamples)
        KernelAbstractions.synchronize(plan.backend)
    end
    return out
end
