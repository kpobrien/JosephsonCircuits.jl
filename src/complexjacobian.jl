

"""
    branchnodesandsigns(Rbnm::SparseMatrixCSC, Nmodes::Integer,
        Nbranches::Integer)

Recover, for each branch, the list of (node, sign) pairs from the incidence
matrix `Rbnm = diagrepeat(Rbn, Nmodes)` which converts node fluxes to branch
fluxes. Verifies that `Rbnm` has the expected mode-diagonal structure (each
entry connects a branch-mode to a node-mode of the same mode index, with the
same value for every mode) and throws an error otherwise.
"""
function branchnodesandsigns(Rbnm::SparseMatrixCSC, Nmodes::Integer,
    Nbranches::Integer)

    size(Rbnm, 1) == Nbranches * Nmodes || throw(DimensionMismatch(
        lazy"Rbnm has $(size(Rbnm,1)) rows, expected Nbranches*Nmodes = $(Nbranches*Nmodes)."))
    size(Rbnm, 2) % Nmodes == 0 || throw(DimensionMismatch(
        lazy"the number of columns of Rbnm is not a multiple of Nmodes."))

    nodesandsigns = [Tuple{Int,eltype(Rbnm)}[] for _ in 1:Nbranches]
    rows = rowvals(Rbnm)
    vals = nonzeros(Rbnm)
    @inbounds for j in axes(Rbnm, 2)
        n = (j - 1) ÷ Nmodes + 1
        mcol = (j - 1) % Nmodes + 1
        for k in nzrange(Rbnm, j)
            i = rows[k]
            b = (i - 1) ÷ Nmodes + 1
            mrow = (i - 1) % Nmodes + 1
            mrow == mcol || throw(ArgumentError(
                "Rbnm does not have the expected mode-diagonal (diagrepeat) structure."))
            if mcol == 1
                push!(nodesandsigns[b], (n, vals[k]))
            else
                # verify every mode carries the same incidence entry
                any(ns -> ns[1] == n && ns[2] == vals[k], nodesandsigns[b]) ||
                    throw(ArgumentError(
                        "Rbnm does not have the expected mode-diagonal (diagrepeat) structure."))
            end
        end
    end
    return nodesandsigns
end


"""
    jjnodeadjacency(Ljb::SparseVector, nodesandsigns, Nnodes::Integer)

For each node, return the sorted, deduplicated list of nodes which share at
least one Josephson junction branch with it (including itself), from the
per-branch (node, sign) lists in `nodesandsigns` (see
[`branchnodesandsigns`](@ref)) restricted to the Josephson branches in `Ljb`.
Used to enumerate the columns of the Jacobian sparsity structure directly in
compressed sparse column form.
"""
function jjnodeadjacency(Ljb::SparseVector, nodesandsigns, Nnodes::Integer)
    adjacency = [Int[] for _ in 1:Nnodes]
    for i in eachindex(Ljb.nzind)
        ns = nodesandsigns[Ljb.nzind[i]]
        for (n2, _) in ns, (n1, _) in ns
            push!(adjacency[n2], n1)
        end
    end
    for a in adjacency
        sort!(a)
        unique!(a)
    end
    return adjacency
end

"""
    activemoderows(Nmodes::Integer, Amatrixindices::Matrix,
        Bmatrixindices::Union{Matrix,Nothing} = nothing)

For each column mode, return the sorted list of row modes for which the
frequency domain index matrix `Amatrixindices` (or, if provided,
`Bmatrixindices`) has a nonzero entry, ie. the row modes which contribute to
the Jacobian sparsity structure in that column mode.
"""
function activemoderows(Nmodes::Integer, Amatrixindices::Matrix,
    Bmatrixindices::Union{Matrix,Nothing} = nothing)
    return [[m1 for m1 in 1:Nmodes if !iszero(Amatrixindices[m1, m2]) ||
        (!isnothing(Bmatrixindices) && !iszero(Bmatrixindices[m1, m2]))]
        for m2 in 1:Nmodes]
end



"""
    ComplexJacobianPlan{TJ}

Precomputed plan for assembling the complex (holomorphic) Jacobian of the
harmonic balance system, from [`plancomplexjacobian`](@ref).

# Fields
- `josephson`: the Josephson contribution, as a linear map from the Fourier
    coefficients of `cos(phi(t))` to the stored entries. This was two segmented
    gathers holding one entry per contribution; it is now the circuit's
    incidence triple product, which is four entries per two terminal junction
    and is read backwards at assembly time. See
    [`StructureComplexJosephsonPlan`](@ref).
- `invLnmindexmap`, `Gnmindexmap`, `Cnmindexmap`: index maps into
    `nonzeros(Jx)` for the frequency dependent linear terms, which are
    scattered at assembly time so the mode frequencies may change between
    assemblies.
"""
struct ComplexJacobianPlan{TJ}
    josephson::TJ
    invLnmindexmap::Vector{Int}
    Gnmindexmap::Vector{Int}
    Cnmindexmap::Vector{Int}
end

"""
    plancomplexjacobian(Amatrixindices::Matrix, Ljb::SparseVector, Lmean,
        Rbnm::SparseMatrixCSC, Nmodes::Integer, Nbranches::Integer,
        Nfreq::Integer, invLnm::SparseMatrixCSC, Gnm::SparseMatrixCSC,
        Cnm::SparseMatrixCSC)

Build the complex Jacobian sparse matrix `Jx` (with the same sparsity
structure `spaddkeepzeros` applied to `Rbnm'*AoLjbm*Rbnm` and the linear term
matrices would produce, including stored numerical zeros) and a
[`ComplexJacobianPlan`](@ref) for assembling it directly with
[`assemblecomplexjacobian!`](@ref).

The plan folds together, at build time, the map from the Fourier coefficients
of `cos(phi(t))` to the Josephson branch matrix `AoLjbm` (`Amatrixindices`,
with negative entries denoting complex conjugation and zeros denoting
dropped couplings) and the incidence matrix triple product
`Rbnm'*AoLjbm*Rbnm`. The frequency dependent linear terms `invLnm`, `Gnm` and
`Cnm` are stored as [`sparseaddmap`](@ref) index maps and scattered at
assembly time so the mode frequencies may change between assemblies.

Returns the tuple `(Jx, plan)`.
"""
function plancomplexjacobian(Amatrixindices::Matrix, Ljb::SparseVector,
    Lmean, Rbnm::SparseMatrixCSC, Nmodes::Integer, Nbranches::Integer,
    Nfreq::Integer, invLnm::SparseMatrixCSC, Gnm::SparseMatrixCSC,
    Cnm::SparseMatrixCSC)

    size(Amatrixindices) == (Nmodes, Nmodes) || throw(DimensionMismatch(
        lazy"Amatrixindices must be Nmodes x Nmodes."))
    isreal(Lmean) || throw(ArgumentError(
        "plancomplexjacobian requires a real Lmean."))
    isempty(Ljb.nzval) || all(isreal, Ljb.nzval) || throw(ArgumentError(
        "plancomplexjacobian requires real Josephson inductances."))

    # a circuit with no Josephson junctions has an Ljb with element type
    # Nothing and an empty nzval, in which case only the linear terms
    # contribute to the Jacobian.
    T = if isempty(Ljb.nzval) || eltype(Ljb) === Nothing
        real(float(typeof(Lmean)))
    else
        real(promote_type(typeof(Lmean), real(eltype(Ljb))))
    end
    Tc = Complex{float(T)}

    nodesandsigns = branchnodesandsigns(Rbnm, Nmodes, Nbranches)

    # Build the sparsity structure (the union of the Josephson contributions
    # and the linear terms, keeping numerical zeros as stored entries so the
    # structure does not change between assemblies) directly in compressed
    # sparse column form. This reproduces the structure produced by
    # spaddkeepzeros on Rbnm'*AoLjbm*Rbnm and the linear term matrices.
    n = size(Rbnm, 2)
    adjacency = jjnodeadjacency(Ljb, nodesandsigns, n ÷ Nmodes)
    activem1 = activemoderows(Nmodes, Amatrixindices)
    C = devicecomplexjacobianpattern(n, Nmodes, adjacency, activem1,
        (invLnm, Gnm, Cnm), CPU())
    colptr, rowval = Array(C.colptr), Array(C.rowval)
    Jx = SparseMatrixCSC(n, n, colptr, rowval, zeros(Tc, length(rowval)))

    # Second pass: emit the runtime scatter entries for the nonlinear part,
    # with destinations resolved to indices into nonzeros(Jx). use 32 bit
    # indices for the scatter lists when the destination and source ranges
    # permit, halving the memory traffic of the plan during both construction
    # and assembly. the branch on the index type is resolved before calling
    # the function barrier fillcomplexscatterlists so the inner loops are
    # type stable.
    josephson = planstructurecomplexjosephson(Jx, T, Amatrixindices, Ljb,
        Lmean, nodesandsigns, Nmodes, Nfreq, CPU())

    # index maps for the frequency dependent linear terms
    invLnmindexmap = sparseaddmap(Jx, invLnm)
    Gnmindexmap = sparseaddmap(Jx, Gnm)
    Cnmindexmap = sparseaddmap(Jx, Cnm)

    plan = ComplexJacobianPlan(josephson, invLnmindexmap, Gnmindexmap,
        Cnmindexmap)

    return Jx, plan
end

"""
    assemblecomplexjacobian!(Jx::SparseMatrixCSC, plan::ComplexJacobianPlan,
        phimatrix::Array, invLnm::SparseMatrixCSC, Gnm::SparseMatrixCSC,
        Cnm::SparseMatrixCSC, wmodesm::Diagonal, wmodes2m::Diagonal)

Assemble the complex (holomorphic) Jacobian `Jx` in place from the Fourier
coefficients of `cos(phi(t))` in `phimatrix` and the linear term matrices,
using a [`ComplexJacobianPlan`](@ref) from [`plancomplexjacobian`](@ref).
Computes

    Jx = Rbnm'*AoLjbm*Rbnm + invLnm + im*Gnm*wmodesm - Cnm*wmodes2m

with no branch matrix update and no sparse matrix multiplications, just a
flat scatter loop for the nonlinear part and index mapped additions
([`sparseadd!`](@ref)) for the linear part.
"""
function assemblecomplexjacobian!(Jx::SparseMatrixCSC,
    plan::ComplexJacobianPlan, phimatrix::Array, invLnm::SparseMatrixCSC,
    Gnm::SparseMatrixCSC, Cnm::SparseMatrixCSC, wmodesm::Diagonal,
    wmodes2m::Diagonal)

    Jxnz = nonzeros(Jx)

    # nonlinear (Josephson) part: a segmented gather from phimatrix, which
    # writes every nonzero, so no zeroing pass is needed
    addjosephsonterm!(Jxnz, plan, phimatrix)

    # frequency dependent linear part
    sparseadd!(Jx, invLnm, plan.invLnmindexmap)
    sparseadd!(Jx, im, Gnm, wmodesm, plan.Gnmindexmap)
    sparseadd!(Jx, -1, Cnm, wmodes2m, plan.Cnmindexmap)

    return Jx
end

"""
    addjosephsonterm!(nzval::AbstractVector, plan::ComplexJacobianPlan,
        phimatrix::Array, conjugate::Bool = false)

Write the Josephson (nonlinear) contribution `Rbnm'*AoLjbm*Rbnm` into the
nonzero value vector `nzval`, which must be aligned with the sparsity
structure the [`ComplexJacobianPlan`](@ref) was built for, by gathering the
Fourier coefficients of `cos(phi(t))` in `phimatrix` through the plan. Every
entry of `nzval` is written, so the destination need not be zeroed first.
With `conjugate = true` the elementwise complex conjugate of the contribution
is written instead, as needed for the adjoint of the linearized harmonic
balance system. Used by [`assemblecomplexjacobian!`](@ref) for the Jacobian of the
nonlinear system and by [`hblinsolve`](@ref) for the pump modulation term of
the linearized system, so the two constructions share the same machinery.
"""
function addjosephsonterm!(nzval::AbstractVector, plan::ComplexJacobianPlan,
    phimatrix::AbstractArray, conjugate::Bool = false)
    return addjosephsonterm!(nzval, plan.josephson, phimatrix, conjugate)
end
