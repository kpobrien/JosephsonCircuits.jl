# The linearized harmonic balance system: the operating point's Josephson
# term and the linear matrices assembled into one sparse matrix per signal
# frequency, and the maps which refill it.

"""
    HBLinearizedSystem

Everything needed to assemble the linearized harmonic balance system matrix
on the signal frequency grid at any signal frequency, built with the same
plan machinery ([`plancomplexjacobian`](@ref) and
[`addjosephsonterm!`](@ref)) used for the Jacobians of the nonlinear system
by [`HBSystem`](@ref). [`hblinsolve`](@ref) constructs one of these and its
per-frequency solves are expressed through it, so any external use of this
object is exercising exactly the production linearized solver path.

Use [`assemblesystemmatrix!`](@ref) to assemble

    A(ws) = AoLjnm + invLnm + im*Gnm*W - Cnm*W^2

(with the negative frequency mode conjugations and the frequency
dependent values resolved, as the linearized solver does) into a matrix
sharing the sparsity
structure of the `Asparse` field, at the signal frequency `ws`, for either
the pump modulation `AoLjnm` or its complex conjugate. Operator products and
adjoint products are then sparse matrix-vector products with the assembled
matrix, `mul!(y, A, v)` and `mul!(y, A', v)`, suitable for iterative solves
and sensitivity adjoints.

Unlike [`HBSystem`](@ref), no time domain matrix-free product is provided,
deliberately: a Fourier transform based product computes a cyclic
convolution on the pump grid, whereas the assembled matrix uses the explicit
truncation of [`hbmatind`](@ref) (the zeros of `Amatrixindices`), dropping
the couplings which fall outside the grid instead of wrapping them. The
assembled matrix defines the linearized solver, so products must match it
exactly; since the pump modulation contribution is precomputed, assembling
at a frequency costs about the same as one matrix-vector product would, and
every product thereafter is a plain sparse matrix-vector product.
"""
struct HBLinearizedSystem{TinvL,TG,TC,TF}
    # sparsity structure of the system matrix, with the pump modulation
    # contribution in its values. per-frequency assembly operates on copies
    # sharing this structure, so this object can be shared across threads.
    Asparse
    # the Josephson map over Asparse ([`StructureComplexJosephsonPlan`](@ref)),
    # which writes the pump modulation contribution and its conjugate
    complexjacobianplan
    # the pump modulation contribution AoLjnm = Rbnm'*AoLjbm*Rbnm, assembled
    # once as a nonzero value vector aligned with the sparsity structure of
    # Asparse
    AoLjnmnzval::Vector{Complex{Float64}}
    # the linear term matrices and their index maps into Asparse; the
    # entries which depend on the frequency are held apart, with their own
    # index maps, in `frequencydependent`, so that the rest are added as
    # the numbers they are. `frequencydependent` is `nothing` when no entry
    # depends on the frequency, and `symbolicvalues` says whether one does,
    # in which case the stored values change with the frequency
    invLnm::TinvL
    Gnm::TG
    Cnm::TC
    invLnmindexmap::Vector{Int}
    Gnmindexmap::Vector{Int}
    Cnmindexmap::Vector{Int}
    frequencydependent::TF
    symbolicvalues::Bool
    # the frequency independent augmentation: the constitutive equations
    # and Kirchhoff current law couplings of the coupled inductor currents
    # and the scattering block port currents, and its index map
    Amna0::SparseMatrixCSC{Complex{Float64},Int}
    Amna0indexmap::Vector{Int}
    # the pump mode frequency offsets of the signal modes, and the numbers
    # of modes and nodes, for computing the mode frequency matrices
    wpumpmodes::Vector{Float64}
    Nmodes::Int
    Nnodes::Int
    # the hybrid stamps of the scattering block components (see
    # ScatteringStampSystem), or nothing
    scattering
end

"""
    HBLinearizedSystem(Amatrixindices::Matrix, Ljb::SparseVector,
        Rbnm::SparseMatrixCSC, Nmodes::Integer, Nbranches::Integer,
        phimatrix::Array, invLnmcopy::SparseMatrixCSC,
        Gnmcopy::SparseMatrixCSC, Cnmcopy::SparseMatrixCSC, invLnm, Gnm,
        Cnm, symbolicvalues::Bool, Amna0::SparseMatrixCSC,
        wpumpmodes, Nnodes::Integer; scattering = nothing)

Construct an [`HBLinearizedSystem`](@ref) from the signal frequency grid
index matrix `Amatrixindices` (see [`hbmatind`](@ref)), the Josephson
junction data, the Fourier coefficients of `cos(phi(t))` of the pump in
`phimatrix`, the numeric copies of the linear term matrices (which define
the sparsity structure) and the possibly symbolic originals, with
`symbolicvalues` saying whether any of their entries is symbolic (see
[`symbolicindices`](@ref)), and the frequency independent augmentation
`Amna0` (the coupled inductor and scattering block port current rows; see
[`calcAmnaind`](@ref)), whose entries are merged into the sparsity
structure. `scattering`
is the [`ScatteringStampSystem`](@ref) of the circuit's scattering blocks,
whose pattern is merged as well, or `nothing`. Builds the sparsity
structure and Josephson map with [`plancomplexjacobian`](@ref) and assembles
the pump modulation contribution (the adjoint takes its conjugate) with
[`addjosephsonterm!`](@ref).
"""
function HBLinearizedSystem(Amatrixindices::Matrix, Ljb::SparseVector,
    Rbnm::SparseMatrixCSC, Nmodes::Integer, Nbranches::Integer,
    phimatrix::Array, invLnmcopy::SparseMatrixCSC,
    Gnmcopy::SparseMatrixCSC, Cnmcopy::SparseMatrixCSC, invLnm, Gnm, Cnm,
    symbolicvalues::Bool, Amna0::SparseMatrixCSC, wpumpmodes,
    Nnodes::Integer; scattering = nothing)

    # the sparsity structure must contain the modified nodal analysis
    # augmentation as well, so merge its entries into the numeric copies
    # used only for the structure. the source index stride is the size of
    # the pump frequency grid, the leading dimensions of phimatrix.
    invLnmpattern = spaddkeepzeros(invLnmcopy, Amna0)
    if !isnothing(scattering)
        # the scattering block contribution shares the structure
        invLnmpattern = spaddkeepzeros(invLnmpattern, scattering.pattern)
    end
    Asparse, complexjacobianplan = plancomplexjacobian(Amatrixindices, Ljb,
        1, Rbnm, Nmodes, Nbranches, prod(size(phimatrix)[1:end-1]),
        invLnmpattern, Gnmcopy, Cnmcopy)

    # the index maps of the matrices added per frequency. A matrix holding
    # a frequency dependent entry holds every entry as an expression, which
    # resolving at each frequency would do entry by entry through dynamic
    # dispatch, so its constant entries are split off as numbers and only
    # the rest are resolved per mode frequency.
    frequencydependent = nothing
    if symbolicvalues
        invLnm, invLnmfd = splitfrequencydependent(invLnm)
        Gnm, Gnmfd = splitfrequencydependent(Gnm)
        Cnm, Cnmfd = splitfrequencydependent(Cnm)
        frequencydependent = (invLnm = invLnmfd,
            invLnmindexmap = sparseaddmap(Asparse, invLnmfd), Gnm = Gnmfd,
            Gnmindexmap = sparseaddmap(Asparse, Gnmfd), Cnm = Cnmfd,
            Cnmindexmap = sparseaddmap(Asparse, Cnmfd))
    end
    invLnmindexmap = sparseaddmap(Asparse, invLnm)
    Gnmindexmap = sparseaddmap(Asparse, Gnm)
    Cnmindexmap = sparseaddmap(Asparse, Cnm)
    Amna0indexmap = sparseaddmap(Asparse, Amna0)
    if !isnothing(scattering)
        setscatteringindexmap!(scattering, Asparse)
    end

    # assemble the pump modulation contribution once, so resetting the
    # system matrix at each frequency is a single copy.
    AoLjnmnzval = zeros(Complex{Float64}, nnz(Asparse))
    addjosephsonterm!(AoLjnmnzval, complexjacobianplan, phimatrix)

    # also place the pump modulation contribution in Asparse so it holds
    # something reasonable to factorize.
    copyto!(Asparse.nzval, AoLjnmnzval)

    return HBLinearizedSystem(Asparse, complexjacobianplan, AoLjnmnzval,
        invLnm, Gnm, Cnm, invLnmindexmap, Gnmindexmap,
        Cnmindexmap, frequencydependent, symbolicvalues, Amna0,
        Amna0indexmap, wpumpmodes, Nmodes, Nnodes, scattering)
end

"""
    splitfrequencydependent(M::SparseMatrixCSC)

The entries of `M` which do not depend on the frequency, as a matrix of
complex numbers, and those which do, as a matrix of the expressions, each
with the stored entries of `M` it holds and no others.
"""
function splitfrequencydependent(M::SparseMatrixCSC)
    I, J, V = findnz(M)
    fd = [CircuitValues.hasprovider(v) for v in V]
    constant = [v isa CircuitValues.Constant ? v.val : Complex{Float64}(v)
        for v in V[.!fd]]
    return sparse(I[.!fd], J[.!fd], constant, size(M)...),
        sparse(I[fd], J[fd], Vector{CircuitValue}(V[fd]), size(M)...)
end

"""
    assemblesystemmatrix!(A::SparseMatrixCSC, lsys::HBLinearizedSystem,
        wmodes::AbstractVector; conjugatepump::Bool = false,
        scatteringwork = ScatteringWorkspace())
    assemblesystemmatrix!(A::SparseMatrixCSC, lsys::HBLinearizedSystem,
        ws::Number; conjugatepump::Bool = false)

Assemble the linearized harmonic balance system matrix into `A`, which must
share the sparsity structure of `lsys.Asparse`, either from the mode
frequency vector `wmodes` or at the signal frequency `ws` (from which the
mode frequencies are computed as `wmodes = ws .+ lsys.wpumpmodes`). The
frequency scaling, the negative frequency conjugation, and any symbolic
frequency substitution are applied per column from the mode index, without
materializing system sized diagonals. With `conjugatepump = true` the complex
conjugate of the pump modulation contribution is used, which for a circuit
without scattering blocks is a similarity transformation of the transposed
system, as below. The
negative frequency mode entries of the linear term matrices are conjugated
and any frequency dependent ones resolved, exactly as in the
per-frequency loop of [`hblinsolve`](@ref), which calls this function.
Returns `A`.

The conjugated pump system is a diagonal similarity transformation of the
transposed forward system,

    A(conjugate pump) = D*transpose(A(pump))*inv(D),

with `D` diagonal, equal to one on every node flux row.
Nothing else contributes, so long as the circuit has no scattering blocks: the
auxiliary rows of the promoted coupled inductors are already symmetric (see
[`calcAmnaind`](@ref)); the linear term matrices are symmetric and mode
diagonal, so the column indexed frequency scaling and conjugation of
`sparseaddconjsubst!` are symmetric under transposition; and the transpose of
the pump modulation contribution exchanges the mode pair, mapping each
difference harmonic to its complex conjugate, which is the same as
conjugating the pump. The hybrid rows of a [`ScatteringParameters`](@ref) do break
it: their constant Kirchhoff couplings and frequency dependent constitutive
entries exchange under transposition, and no diagonal `D` undoes that, so with
blocks the two are different matrices.

The adjoint solutions the noise, quantum efficiency, commutation relation and
adjoint node output calculations read are in every case the solutions of the
*transposed* system: by the adjoint identity, the response at an output port
to a source anywhere in the circuit is that source contracted against the
transposed solution driven at the port. [`hblinsolve`](@ref) obtains them with
[`trysolvetranspose!`](@ref) on the factorization of the forward system, which
costs a pair of triangular solves rather than an assembly and a factorization
at every signal frequency. Where the similarity holds, the conjugated pump
assembly is the independent construction that equivalence is tested against.
`scatteringwork` is the scratch of the scattering block evaluation, one
per worker in a sweep, so that a frequency allocates nothing for it.
"""
function assemblesystemmatrix!(A::SparseMatrixCSC,
    lsys::HBLinearizedSystem, wmodes::AbstractVector;
    conjugatepump::Bool = false,
    scatteringwork::ScatteringWorkspace = ScatteringWorkspace())

    # the pump modulation contribution, precomputed, or its complex
    # conjugate
    if conjugatepump
        A.nzval .= conj.(lsys.AoLjnmnzval)
    else
        copyto!(A.nzval, lsys.AoLjnmnzval)
    end

    # take the complex conjugate of the negative frequency terms in
    # the capacitance and conductance matrices, and resolve a frequency
    # dependent entry at the mode frequency of its column. the frequency
    # scaling of each term is the per column mode frequency raised to the
    # given power.
    sparseaddconjsubst!(A, -1, lsys.Cnm, lsys.Cnmindexmap, wmodes, 2)
    sparseaddconjsubst!(A, im, lsys.Gnm, lsys.Gnmindexmap, wmodes, 1)
    sparseaddconjsubst!(A, 1, lsys.invLnm, lsys.invLnmindexmap, wmodes, 0)
    fd = lsys.frequencydependent
    if !isnothing(fd)
        sparseaddconjsubst!(A, -1, fd.Cnm, fd.Cnmindexmap, wmodes, 2)
        sparseaddconjsubst!(A, im, fd.Gnm, fd.Gnmindexmap, wmodes, 1)
        sparseaddconjsubst!(A, 1, fd.invLnm, fd.invLnmindexmap, wmodes, 0)
    end

    # the frequency independent augmentation: the coupled inductor and
    # scattering block port current rows. A complex (lossy) coupled
    # inductance is a value like any other, conjugated in the columns of a
    # negative frequency mode; the incidence entries are real.
    sparseaddconjsubst!(A, 1, lsys.Amna0, lsys.Amna0indexmap, wmodes, 0)

    # the hybrid stamps of the scattering block components (see
    # `ScatteringStampSystem`): `sign*im*w_m*scale*B[p,q](w_m)` on the node
    # flux columns and `-iscale*C[p,q](w_m)` on the auxiliary port current
    # columns of their constitutive rows, per contribution
    if !isnothing(lsys.scattering)
        assemblescattering!(A, lsys.scattering, wmodes, scatteringwork)
    end

    return A
end

function assemblesystemmatrix!(A::SparseMatrixCSC,
    lsys::HBLinearizedSystem, ws::Number; conjugatepump::Bool = false)
    return assemblesystemmatrix!(A, lsys, ws .+ lsys.wpumpmodes;
        conjugatepump = conjugatepump)
end
