

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
    plancomplexjacobian(Amatrixindices::Matrix, Ljb::SparseVector, Lscale,
        Rbnm::SparseMatrixCSC, Nmodes::Integer, Nbranches::Integer,
        Nfreq::Integer, invLnm::SparseMatrixCSC, Gnm::SparseMatrixCSC,
        Cnm::SparseMatrixCSC)

Build the complex Jacobian sparse matrix `Jx` (with the same sparsity
structure `spaddkeepzeros` applied to `Rbnm'*AoLjbm*Rbnm` and the linear term
matrices would produce, including stored numerical zeros) and the
[`StructureComplexJosephsonPlan`](@ref) which writes the Josephson term into
it: the map from the Fourier coefficients of `cos(phi(t))` to the Josephson
branch matrix `AoLjbm` (`Amatrixindices`, with negative entries denoting
complex conjugation and zeros denoting dropped couplings) and the circuit's
incidence triple product `Rbnm'*AoLjbm*Rbnm` as a per node pair table
([`junctionpairtable`](@ref)), read backwards at assembly time.

The plan reads a [`JunctionStructure`](@ref) built here on the host in the
precision the scale and the inductances give. The linear term is not part
of the plan. The nonlinear solve gathers it once for its mode frequencies
with [`planstructurecomplexjacobian`](@ref), which takes this plan as its
`josephson`; the linearized solve builds its own index maps for the
matrices it substitutes per frequency.

Returns the tuple `(Jx, josephson)`.
"""
function plancomplexjacobian(Amatrixindices::Matrix, Ljb::SparseVector,
    Lscale, Rbnm::SparseMatrixCSC, Nmodes::Integer, Nbranches::Integer,
    Nfreq::Integer, invLnm::SparseMatrixCSC, Gnm::SparseMatrixCSC,
    Cnm::SparseMatrixCSC)

    size(Amatrixindices) == (Nmodes, Nmodes) || throw(DimensionMismatch(
        lazy"Amatrixindices must be Nmodes x Nmodes."))
    isreal(Lscale) || throw(ArgumentError(
        "plancomplexjacobian requires a real Lscale."))
    isempty(Ljb.nzval) || all(isreal, Ljb.nzval) || throw(ArgumentError(
        "plancomplexjacobian requires real Josephson inductances."))

    # a circuit with no Josephson junctions has an empty Ljb, in which case
    # only the linear terms contribute to the Jacobian
    T = real(promote_type(float(typeof(Lscale)), real(eltype(Ljb))))
    Tc = Complex{float(T)}

    junctions = junctionstructure(T, Amatrixindices,
        zeros(Int, Nmodes, Nmodes), Ljb, Lscale, Rbnm, Nmodes, Nbranches,
        Nfreq, CPU())
    nodesandsigns = junctions.nodesandsigns

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

    # Second pass: the Josephson term as a structure plan over the pattern
    # just built, which reads the incidence triple product backwards at
    # assembly time
    josephson = planstructurecomplexjosephson(Jx, junctions, CPU())

    return Jx, josephson
end
