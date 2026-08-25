"""
    JosephsonCircuitsCUDAExt

Package extension loaded with `using CUDA`. It supplies the one thing the
residual and the matrix-free Jacobian-vector and Hessian-vector products of
[`JosephsonCircuits.HBSystem`](@ref) need on a CUDA device that the core
package cannot provide without taking on the dependency: the real transform
plans, through CUFFT via the AbstractFFTs interface.

Everything else on that path is already device generic. The two linear maps
around the pointwise time domain nonlinearity run as KernelAbstractions
kernels of `NonlinearTermPlan` on whichever backend the `HBSystem` was built
for, one work item per output slot with no atomics; the pointwise and
normalization steps are broadcasts; and the transform is a batched real
transform over all but the last dimension, whose layout (the Josephson
junction index last, and so slowest varying in column major) is already the
one CUFFT batches over.

The Jacobian assembly is a separate matter and is not addressed here. Its
scatters were converted to segmented gathers, so it has the right shape, but
it is still written as host loops and its only consumer factorizes the
result, which is a CPU factorization. Supply a device `Factorization` triple
to `hbnlsolve` before expecting that path to run on a device.

!!! warning
    This extension has NOT been executed on a GPU by its author. The `CPU()`
    backend of the same kernels is the correctness oracle, and it is
    validated against the assembled Jacobian and against central finite
    differences of the residual. Test on hardware before relying on results.
"""
module JosephsonCircuitsCUDAExt

using CUDA
using CUDA.CUFFT
using KernelAbstractions
import JosephsonCircuits: fftplans

# Real transform plans on the device, with the same dimensions, direction and
# normalization convention as the FFTW plans of the CPU backend: the
# transform runs over all but the last dimension, and the caller applies the
# prod(size(td)[1:end-1]) normalization.
function fftplans(fd::AbstractArray{Complex{T}}, td::AbstractArray{T},
    stepsperperiod::Int, backend::CUDABackend) where T
    dims = 1:length(size(fd))-1
    irfftplan = CUFFT.plan_irfft(fd, stepsperperiod, dims)
    rfftplan = CUFFT.plan_rfft(td, dims)
    return irfftplan, rfftplan
end


end # module
