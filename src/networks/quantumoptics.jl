
"""
    symplectic_form_block(n::Integer)

Return the `2n x 2n` matrix representing the symplectic form `Ω` for `n` modes
in the real quadrature operator basis with block order
`r = [x_1,...,x_n,p_1,...,p_n]` where `Ω = [0_n 1_n;-1_n 0_n]`. `0_n` is an
`n` by `n` matrix of zeros and `1_n` is an `n` by `n` identity matrix.

# Examples
```jldoctest
julia> JosephsonCircuits.symplectic_form_block(2)
4×4 SparseArrays.SparseMatrixCSC{Int64, Int64} with 4 stored entries:
  ⋅   ⋅  1  ⋅
  ⋅   ⋅  ⋅  1
 -1   ⋅  ⋅  ⋅
  ⋅  -1  ⋅  ⋅
```
"""
function symplectic_form_block(n::Integer)
    # column pointer has length number of columns + 1
    colptr = [i for i in 1:2*n+1]

    #
    rowval = Vector{Int}(undef, 2 * n)
    nzval = Vector{Int}(undef, 2 * n)

    for i in 1:n
        rowval[i] = n + i
        nzval[i] = -1
    end
    for i in n+1:2*n
        rowval[i] = i - n
        nzval[i] = 1
    end

    return SparseMatrixCSC(2 * n, 2 * n, colptr, rowval, nzval)
end

"""
    direct_sum(A; n::Integer=1)

Return the direct sum of `n` copies of the matrix `A`, the block diagonal
matrix `kron(I(n), A)`.
"""
function direct_sum(A; n::Integer=1)
    return kron(I(n), A)
end

"""
    symplectic_form_pair(n::Integer)

Return the `2n x 2n` matrix representing the symplectic form `Ω` for `n` modes
in the real quadrature operator basis with pair order
`r = [x_1,p_1,...,x_n,p_n]` where `Ω` = direct sum of `n` of `Ω1` where
`Ω1 = [0 1; -1 0]`.

# Examples
```jldoctest
julia> JosephsonCircuits.symplectic_form_pair(2)
4×4 SparseArrays.SparseMatrixCSC{Int64, Int64} with 4 stored entries:
  ⋅  1   ⋅  ⋅
 -1  ⋅   ⋅  ⋅
  ⋅  ⋅   ⋅  1
  ⋅  ⋅  -1  ⋅
```
"""
function symplectic_form_pair(n::Integer)
    Omega = sparse([0 1; -1 0])
    return direct_sum(Omega; n=n)
end

"""
    indefinite_hermitian_form_pair(n)

Return the `2n x 2n` matrix representing the indefinite Hermitian form `Σ` for
`n` modes in the annihilation and creation operator basis with pair
order `ξ = [a_1,adag_1,...,a_n,adag_n]` where `Σ` = direct sum of `n` of `σ3`
where `σ3 = [1 0; 0 -1]`.

# Examples
```jldoctest
julia> JosephsonCircuits.indefinite_hermitian_form_pair(2)
4×4 SparseArrays.SparseMatrixCSC{Int64, Int64} with 4 stored entries:
 1   ⋅  ⋅   ⋅
 ⋅  -1  ⋅   ⋅
 ⋅   ⋅  1   ⋅
 ⋅   ⋅  ⋅  -1
```
"""
function indefinite_hermitian_form_pair(n::Integer)
    Sigma = sparse([1 0; 0 -1])
    return direct_sum(Sigma; n=n)
end

"""
    indefinite_hermitian_form_block(n::Int)

Return the `2n x 2n` matrix representing the indefinite Hermitian form `Σ` for
`n` modes in the annihilation and creation operator basis with block order
`ξ = [a_1,...,a_n,adag_1,...,adag_n]` where `Σ = [1_n 0_n;0_n -1_n]`. `0_n` is
an `n` by `n` matrix of zeros and `1_n` is an `n` by `n` identity matrix.

# Examples
```jldoctest
julia> using LinearAlgebra

julia> JosephsonCircuits.indefinite_hermitian_form_block(2)
4×4 Diagonal{Int64, Vector{Int64}}:
 1  ⋅   ⋅   ⋅
 ⋅  1   ⋅   ⋅
 ⋅  ⋅  -1   ⋅
 ⋅  ⋅   ⋅  -1
```
"""
function indefinite_hermitian_form_block(n::Integer)
    d = Vector{Int}(undef, 2 * n)
    for i in 1:n
        d[i] = 1
    end
    for i in n+1:2*n
        d[i] = -1
    end
    return Diagonal(d)
end

# the default relative tolerance of the approximate checks, that of
# `isapprox`: the square root of the machine epsilon of the element type
# `T`, and zero when an absolute tolerance is given
approxrtol(T, atol) = sqrt(eps(real(float(oneunit(T))))) * iszero(atol)

"""
    is_positive_semi_definite(M) -> Bool

Return `true` if the matrix `M` is positive semi-definite and `false`
otherwise.

"""
function is_positive_semi_definite(M)
    # a pivoted Cholesky with error checking off, so that a rank deficient
    # matrix does not throw; positive semidefiniteness is then judged by
    # whether the factorization of the rank it reports reproduces the
    # matrix, in its pivoted order, to numerical error
    if !isapprox(M, M')
        return false
    end
    C = cholesky(Hermitian(M), RowMaximum(); check=false)
    issuccess(C) && return true
    L = C.L[:, 1:C.rank]
    return isapprox(M[C.p, C.p], L * L')
end

"""
    is_positive_definite(M) -> Bool

Return `true` if the matrix `M` is positive definite and `false` otherwise.

"""
function is_positive_definite(M)
    if isapprox(M, M')
        return isposdef(Hermitian(M))
    else
        return false
    end
end


"""
    is_unitary(M) -> Bool

Return `true` if the matrix `M` is unitary, `M ∈ U(n)`, and `false` otherwise.

Tests if `M` is square and satisfies the condition `M*M'==I` where `I` is
the identity matrix. A rectangular matrix with orthonormal rows, such as
`[1 0]`, is not unitary.

"""
function is_unitary(M)
    size(M, 1) == size(M, 2) || return false
    return isapprox(M * adjoint(M), I(size(M, 1)))
end

"""
    is_orthogonal(M) -> Bool

Return `true` if the matrix `M` is orthogonal, `M ∈ O(n)`, and `false`
otherwise.

Tests if `M` is square and satisfies the condition `M*transpose(M)==I`
where `I` is the identity matrix. A rectangular matrix with orthonormal
rows, such as `[1 0]`, is not orthogonal.

"""
function is_orthogonal(M)
    size(M, 1) == size(M, 2) || return false
    return isapprox(M * transpose(M), I(size(M, 1)))
end

"""
    is_symplectic(Ω, S) -> Bool

Return `true` if the matrix `S` is symplectic, `S ∈ Sp(2n, ℝ)` or
`S ∈ Sp(2n, ℂ)`, and `false` otherwise.

Tests if `S` satisfies the symplectic condition `S*Ω*transpose(S)==Ω` where
`Ω` is a user supplied symplectic form.

See also [`symplectic_form_block`](@ref),
[`symplectic_form_pair`](@ref), [`is_symplectic_block`](@ref), and
[`is_symplectic_pair`](@ref).

"""
function is_symplectic(Ω, S)
    if isodd(size(S, 1)) || isodd(size(S, 2))
        error(lazy"The dimensions of the input matrix must be even.")
    end
    # Ω is not checked to be a symplectic form, so swapped arguments are
    # not detected
    return isapprox(S * Ω * transpose(S), Ω)
end

"""
    is_symplectic_block(S) -> Bool

Return `true` if the matrix `S` is symplectic, `S ∈ Sp(2n, ℝ)` or
`S ∈ Sp(2n, ℂ)`, with block operator order and `false` otherwise.

Tests if `S` satisfies the symplectic condition `S*Ω*transpose(S)==Ω` where
`Ω` is the matrix representing the symplectic form with block operator order.

See also [`symplectic_form_block`](@ref).

"""
function is_symplectic_block(S)
    Omega = symplectic_form_block(size(S, 1) ÷ 2)
    return is_symplectic(Omega, S)
end

"""
    is_symplectic_pair(S) -> Bool

Return `true` if the matrix `S` is symplectic, `S ∈ Sp(2n, ℝ)` or
`S ∈ Sp(2n, ℂ)`, with pair operator order and `false` otherwise.

Tests if `S` satisfies the symplectic condition `S*Ω*transpose(S)==Ω` where
`Ω` is the matrix representing the symplectic form with pair operator
order.

See also [`symplectic_form_pair`](@ref).

"""
function is_symplectic_pair(S)
    Omega = symplectic_form_pair(size(S, 1) ÷ 2)
    return is_symplectic(Omega, S)
end

"""
    is_orthogonal_symplectic_block(M) -> Bool

Return `true` if the matrix `M` is orthogonal symplectic,
`M ∈ Sp(2n, ℝ) ∩ O(2n) ≅ U(n)`, with block operator order and `false`
otherwise.

"""
function is_orthogonal_symplectic_block(M)
    return is_symplectic_block(M) && is_orthogonal(M)
end

"""
    is_orthogonal_symplectic_pair(M) -> Bool

Return `true` if the matrix `M` is orthogonal symplectic,
`M ∈ Sp(2n, ℝ) ∩ O(2n) ≅ U(n)`, with pair operator order and `false`
otherwise.

"""
function is_orthogonal_symplectic_pair(M)
    return is_symplectic_pair(M) && is_orthogonal(M)
end

"""
    is_conjugate_symplectic_pair(M) -> Bool

Return `true` if the matrix `M` is conjugate symplectic, `M*Ω*M' == Ω`
with `Ω` the symplectic form of pair operator order, and `false`
otherwise.

"""
function is_conjugate_symplectic_pair(M)
    Omega = symplectic_form_pair(size(M, 1) ÷ 2)
    return is_pseudo_unitary(Omega, M)
end

"""
    is_conjugate_symplectic_block(M) -> Bool

Return `true` if the matrix `M` is conjugate symplectic, `M*Ω*M' == Ω`
with `Ω` the symplectic form of block operator order, and `false`
otherwise.

"""
function is_conjugate_symplectic_block(M)
    Omega = symplectic_form_block(size(M, 1) ÷ 2)
    return is_pseudo_unitary(Omega, M)
end


"""
    is_pseudo_unitary(Σ,M) -> Bool

Return `true` if the matrix `M` is pseudo-unitary, `M ∈ U(n, n)`, and `false`
otherwise.

Tests if `M` satisfies the pseudo-unitary condition `M*Σ*M'==Σ` where `Σ` is a
user supplied matrix representing the indefinite Hermitian form.

See also [`indefinite_hermitian_form_pair`](@ref),
[`indefinite_hermitian_form_block`](@ref), [`is_pseudo_unitary_block`](@ref),
and [`is_pseudo_unitary_pair`](@ref).

"""
function is_pseudo_unitary(Sigma, M)
    return isapprox(M * Sigma * M', Sigma)
end

"""
    is_pseudo_unitary_block(M) -> Bool

Return `true` if the matrix `M` is pseudo-unitary, `M ∈ U(n, n)`, with block
operator order and `false` otherwise.

Tests if `M*Σ*M'==Σ` where `Σ` is a matrix representing the indefinite
Hermitian form with block operator order.

See also [`indefinite_hermitian_form_block`](@ref).

"""
function is_pseudo_unitary_block(M)
    Sigma = indefinite_hermitian_form_block(size(M, 1) ÷ 2)
    return is_pseudo_unitary(Sigma, M)
end

"""
    is_pseudo_unitary_pair(M) -> Bool

Return `true` if the matrix `M` is pseudo-unitary, `M ∈ U(n, n)`, with
pair operator order and `false` otherwise.

Tests if `M*Σ*M'==Σ` where `Σ` is a matrix representing the indefinite
Hermitian form with pair operator order.

See also [`indefinite_hermitian_form_pair`](@ref).

"""
function is_pseudo_unitary_pair(M)
    Sigma = indefinite_hermitian_form_pair(size(M, 1) ÷ 2)
    return is_pseudo_unitary(Sigma, M)
end


"""
    is_positive_definite_symplectic_block(S) -> Bool

Return `true` if the matrix `S` is positive definite and symplectic,
`S ∈ Sp(2n, ℝ)` or `S ∈ Sp(2n, ℂ)`, with block operator order and `false`
otherwise.

"""
function is_positive_definite_symplectic_block(S)
    # `ishermitian` tests exact equality; test up to numerical error instead
    return is_symplectic_block(S) && isapprox(S, S') && is_positive_definite(Hermitian(S))
end

"""
    is_positive_definite_symplectic_pair(S) -> Bool

Return `true` if the matrix `S` is positive definite and symplectic,
`S ∈ Sp(2n, ℝ)` or `S ∈ Sp(2n, ℂ)`, with pair operator order and
`false` otherwise.

"""
function is_positive_definite_symplectic_pair(S)
    # `ishermitian` tests exact equality; test up to numerical error instead
    return is_symplectic_pair(S) && isapprox(S, S') && is_positive_definite(Hermitian(S))
end

"""
    is_bogoliubov_pair(M) -> Bool

Return `true` if the matrix `M` is Bogoliubov, `M ∈ Sp(2n, ℂ) ∩ U(n, n)`, with
pair operator order and `false` otherwise.

"""
function is_bogoliubov_pair(M)
    return is_pseudo_unitary_pair(M) && is_symplectic_pair(M)
end

"""
    is_bogoliubov_block(M) -> Bool

Return `true` if the matrix `M` is Bogoliubov, `M ∈ Sp(2n, ℂ) ∩ U(n, n)`, with
block operator order and `false` otherwise.

"""
function is_bogoliubov_block(M)
    return is_pseudo_unitary_block(M) && is_symplectic_block(M)
end

"""
    is_orthogonal_bogoliubov_block(M) -> Bool

Return `true` if the matrix `M` is orthogonal Bogoliubov,
`M ∈ Sp(2n, ℂ) ∩ U(n, n) ∩ U(2n) ≅ U(n)`, with block operator order and
`false` otherwise.

"""
function is_orthogonal_bogoliubov_block(M)
    return is_bogoliubov_block(M) && is_unitary(M)
end

"""
    is_orthogonal_bogoliubov_pair(M) -> Bool

Return `true` if the matrix `M` is orthogonal Bogoliubov,
`M ∈ Sp(2n, ℂ) ∩ U(n, n) ∩ U(2n) ≅ U(n)`, with pair operator order and
`false` otherwise.

"""
function is_orthogonal_bogoliubov_pair(M)
    return is_bogoliubov_pair(M) && is_unitary(M)
end


# the scale `hbar` of the quadrature functions, finite and positive; see
# `is_cptp`
function checkhbar(hbar::Real)
    if !(isfinite(hbar) && hbar > 0)
        throw(ArgumentError(lazy"`hbar` must be finite and positive; got $(hbar)."))
    end
    return hbar
end

"""
    is_cptp(Omega, X, Y; hbar = 1, atol = 0, rtol = ...) -> Bool

Return `true` if the Gaussian map `V -> X*V*X' + Y` of the transformation
`X` and the noise `Y` is completely positive and trace preserving for the
symplectic form `Omega`, that is if
`K = Y + im*(hbar/2)*(Omega - X*Omega*X')` is Hermitian and positive
semi-definite, and `false` otherwise.

The quadrature functions share one convention, with the scale `hbar`. The
quadratures of a mode with the ladder operators `[a, adag] = 1` are
`x = sqrt(hbar/2)*(a + adag)` and `p = -im*sqrt(hbar/2)*(a - adag)`, with
`[x, p] = im*hbar`; the covariance of the quadratures `r` is
`V = ⟨{Δr, Δrᵀ}⟩/2`, and that of the vacuum `(hbar/2)*I`. The default
`hbar = 1` counts the vacuum as half a photon, as circuit quantum
electrodynamics does; `hbar = 2` gives the covariance `σ` of Serafini,
whose vacuum is `I`, and the condition of his Eq. 5.37. `hbar` must be
finite and positive.

`K` is a difference of terms which cancel, exactly so for a noiseless map
(`Y = 0` with `X` symplectic), so it is judged against their scale
`s = norm(Y) + (hbar/2)*norm(Omega)*(1 + opnorm(X)^2)` rather than against
itself: `K` is Hermitian when `norm(K - K') <= tol` and positive
semi-definite when the smallest eigenvalue of its Hermitian part is at least
`-tol`, with `tol = max(atol, rtol*s)`. The default `rtol` is the square root
of the machine epsilon of the element type, and zero when `atol` is given.

# References
A. Serafini, "Quantum Continuous Variables: A Primer of Theoretical
Methods," CRC Press (2017).
"""
function is_cptp(Omega, X, Y; hbar::Real = 1, atol::Real = 0,
        rtol::Real = approxrtol(promote_type(eltype(X), eltype(Y)), atol))
    K, tol = cptpcondition(Omega, X, Y, hbar, atol, rtol)
    if norm(K - K') > tol
        return false
    end
    return isempty(K) || eigmin(Hermitian((K + K') / 2)) >= -tol
end

# the matrix `K` of the condition of `is_cptp` and the tolerance `tol` it
# is judged to, which the dilation of the map shares
function cptpcondition(Omega, X, Y, hbar, atol, rtol)
    checkhbar(hbar)
    if size(X) != size(Omega) || size(Y) != size(Omega) || size(Omega, 1) != size(Omega, 2)
        throw(DimensionMismatch(lazy"`Omega`, `X` and `Y` must be square matrices of one size; got $(size(Omega)), $(size(X)) and $(size(Y))."))
    end
    Omegam = Matrix(Omega)
    K = Y + im * (hbar / 2) * (Omegam - X * Omegam * X')
    tol = max(atol, rtol * (norm(Y) + (hbar / 2) * norm(Omegam) * (1 + opnorm(X)^2)))
    return K, tol
end

"""
    is_cptp_quadrature_block(X, Y; hbar = 1, atol = 0, rtol = ...) -> Bool

[`is_cptp`](@ref) for the quadrature transformation `X` and noise `Y` in
block operator order, in the units of `hbar` stated there.
"""
function is_cptp_quadrature_block(X, Y; kwargs...)
    n = size(X, 1) ÷ 2
    Omega = symplectic_form_block(n)
    return is_cptp(Omega, X, Y; kwargs...)
end

"""
    is_cptp_quadrature_pair(X, Y; hbar = 1, atol = 0, rtol = ...) -> Bool

[`is_cptp`](@ref) for the quadrature transformation `X` and noise `Y` in
pair operator order, in the units of `hbar` stated there.
"""
function is_cptp_quadrature_pair(X, Y; kwargs...)
    n = size(X, 1) ÷ 2
    Omega = symplectic_form_pair(n)
    return is_cptp(Omega, X, Y; kwargs...)
end

"""
    is_cptp_ladder_pair(X, Y; atol = 0, rtol = ...) -> Bool

[`is_cptp`](@ref) for the ladder (Bogoliubov) transformation `X` and noise
`Y` in pair operator order. The ladder operators are dimensionless,
`[a, adag] = 1`, and `Y` is a symmetrized covariance `⟨{Δξ, Δξ†}⟩/2` of
`ξ = [a_1, adag_1, ...]`, whose vacuum is `I/2`;
[`ladder_to_quadrature_pair`](@ref) takes it to the quadrature covariance
at `hbar = 1`.
"""
function is_cptp_ladder_pair(X, Y; atol::Real = 0,
        rtol::Real = approxrtol(promote_type(eltype(X), eltype(Y)), atol))
    n = size(X, 1) ÷ 2
    # the ladder image of the quadrature symplectic form is -im*Σ, with Σ
    # the indefinite Hermitian form; +im*Σ, used here, gives the complex
    # conjugate of the condition in the quadrature basis, which has the same
    # eigenvalues, so the verdict is the same
    Omega = im * indefinite_hermitian_form_pair(n)
    return is_cptp(Omega, X, Y; hbar = 1, atol = atol, rtol = rtol)
end

"""
    is_cptp_ladder_block(X, Y; atol = 0, rtol = ...) -> Bool

[`is_cptp_ladder_pair`](@ref) in block operator order.
"""
function is_cptp_ladder_block(X, Y; atol::Real = 0,
        rtol::Real = approxrtol(promote_type(eltype(X), eltype(Y)), atol))
    n = size(X, 1) ÷ 2
    # as in is_cptp_ladder_pair
    Omega = im * indefinite_hermitian_form_block(n)
    return is_cptp(Omega, X, Y; hbar = 1, atol = atol, rtol = rtol)
end

"""
    rand_positive_definite(T, n::Integer)
    rand_positive_definite(n::Integer)

A random `n` by `n` positive definite matrix of element type `T` (`Float64`
by default), `A*A'` for a random matrix `A`; see
[`rand_positive_semi_definite`](@ref).
"""
function rand_positive_definite(T, n::Integer)
    A = rand(T, n, n)
    return A * A'
end

function rand_positive_definite(n::Integer)
    return rand_positive_definite(Float64, n)
end

"""
    rand_unitary(T, n::Integer)
    rand_unitary(n::Integer)

A random `n` by `n` unitary matrix, the exponential of the skew-Hermitian
part of a random matrix of element type `T` (`ComplexF64` by default); see
[`is_unitary`](@ref).
"""
function rand_unitary(T, n::Integer)
    A = rand(T, n, n)
    # make a skew-Hermitian matrix
    H = (A - A') / 2
    return exp(H)
end

function rand_unitary(n::Integer)
    return rand_unitary(Complex{Float64}, n)
end

"""
    rand_orthogonal(n::Integer)

A random `n` by `n` orthogonal matrix, [`rand_unitary`](@ref) of a real
matrix.
"""
function rand_orthogonal(n::Integer)
    return rand_unitary(Float64, n)
end


"""
    rand_positive_semi_definite(T, n, m)
    rand_positive_semi_definite(n, m)

A random `(n + m)` by `(n + m)` positive semidefinite matrix of rank at
most `n` and element type `T` (`Float64` by default), as `A*A'` for a
random `(n + m)` by `n` matrix `A` with entries uniform in `[0, 1)`.

"""
function rand_positive_semi_definite(T, n, m)
    A = rand(T, n + m, n)
    return A * A'
end

function rand_positive_semi_definite(n, m)
    return rand_positive_semi_definite(Float64, n, m)
end


"""
    rand_symplectic_block(T, n::Integer)

Return a random `2n x 2n` symplectic matrix `S`, `S ∈ Sp(2n, ℝ)` or
`S ∈ Sp(2n, ℂ)`, depending on the type `T` with block operator order.

"""
function rand_symplectic_block(T::DataType, n::Integer)
    A = randn(T, 2 * n, 2 * n)
    Omega = symplectic_form_block(n)
    # generate a random symmetric matrix
    M = (transpose(A) + A) / 2
    # the matrix exponential is numerically more accurate than the Cayley
    # transform here
    return exp(Omega * M)
end

"""
    rand_symplectic_block(n::Integer)

Return a random `2n x 2n` symplectic matrix `S`, `S ∈ Sp(2n, ℝ)`, with block
operator order.

"""
function rand_symplectic_block(n::Integer)
    return rand_symplectic_block(Float64, n)
end

function rand_symplectic_pair(T::DataType, n::Integer)
    A = randn(T, 2 * n, 2 * n)
    Omega = symplectic_form_pair(n)
    # generate a random symmetric matrix
    M = (transpose(A) + A) / 2
    # the matrix exponential is numerically more accurate than the Cayley
    # transform here
    return exp(Omega * M)
end

"""
    rand_symplectic_pair(n::Integer)

Return a random `2n x 2n` symplectic matrix `S`, `S ∈ Sp(2n, ℝ)`, with
pair operator order.

"""
function rand_symplectic_pair(n::Integer)
    return rand_symplectic_pair(Float64, n)
end

function rand_orthogonal_symplectic_block(T, n::Integer)
    P = randn(T, n, n)
    # make P skew-Symmetric
    P = (P - transpose(P)) / 2

    Q = randn(T, n, n)
    # make Q complex symmetric
    Q = (Q + transpose(Q)) / 2

    # assemble M
    M = [-Q -P; P -Q]

    # compute the symplectic matrix using the Cayley transform
    Omega = symplectic_form_block(n)
    return cayley_transform(Omega, M)
end

"""
    rand_orthogonal_symplectic_block(n::Integer)

Return a random `2n x 2n` orthogonal symplectic matrix `S`,
`S ∈ Sp(2n, ℝ) ∩ O(2n) ≅ U(n)`, with block operator order.

"""
function rand_orthogonal_symplectic_block(n::Integer)
    return rand_orthogonal_symplectic_block(Float64, n)
end


"""
    rand_orthogonal_symplectic_pair(n::Integer)

Return a random `2n x 2n` orthogonal symplectic matrix `S`,
`S ∈ Sp(2n, ℝ) ∩ O(2n) ≅ U(n)` with pair operator order.

"""
function rand_orthogonal_symplectic_pair(n::Integer)
    return rand_orthogonal_symplectic_pair(Float64, n)
end

function rand_orthogonal_symplectic_pair(T, n::Integer)
    return block_to_pair(rand_orthogonal_symplectic_block(T, n))
end


"""
    rand_positive_definite_symplectic_block(T, n::Integer)

Return a random `2n x 2n` positive definite symplectic matrix `S`,
`S ∈ Sp(2n, ℝ)` or `S ∈ Sp(2n, ℂ)`, depending on the type `T` with block
operator order.

"""
function rand_positive_definite_symplectic_block(T, n::Integer)

    # the positive definite factor of the polar decomposition of a random
    # symplectic matrix, which is more robust numerically than an explicit
    # construction from a random positive definite matrix
    A = rand_symplectic_block(T, n)
    S, _ = polar(A)
    return S


end

"""
    rand_positive_definite_symplectic_block(n::Integer)

Return a random `2n x 2n` positive definite symplectic matrix `S`,
`S ∈ Sp(2n, ℝ)`, with block operator order.

"""
function rand_positive_definite_symplectic_block(n::Integer)
    return rand_positive_definite_symplectic_block(Float64, n)
end

"""
    rand_positive_definite_symplectic_pair(T, n::Integer)

Return a random `2n x 2n` positive definite symplectic matrix `S`,
`S ∈ Sp(2n, ℝ)` or `S ∈ Sp(2n, ℂ)`, depending on the type `T` with pair
operator order.

"""
function rand_positive_definite_symplectic_pair(T, n::Integer)
    return block_to_pair(rand_positive_definite_symplectic_block(T, n))
end

"""
    rand_positive_definite_symplectic_pair(n::Integer)

Return a random `2n x 2n` positive definite symplectic matrix `S`,
`S ∈ Sp(2n, ℝ)`, with pair operator order.

"""
function rand_positive_definite_symplectic_pair(n::Integer)
    return rand_positive_definite_symplectic_pair(Float64, n)
end


"""
    rand_conjugate_symplectic_block(T, n::Integer)

Return a random `2n x 2n` conjugate symplectic matrix `S`, `S*Ω*S' == Ω`,
of element type `T`, with block operator order.
"""
function rand_conjugate_symplectic_block(T, n::Integer)
    A = randn(T, 2 * n, 2 * n)
    Omega = symplectic_form_block(n)
    # generate a random Hermitian matrix
    M = (A + A') / 2
    # the matrix exponential is numerically more accurate than the Cayley
    # transform here
    return exp(Omega * M)

end

function rand_conjugate_symplectic_block(n::Integer)
    return rand_conjugate_symplectic_block(Complex{Float64}, n)
end

"""
    rand_conjugate_symplectic_pair(T, n::Integer)
    rand_conjugate_symplectic_pair(n::Integer)

The pair ordered form of [`rand_conjugate_symplectic_block`](@ref).
"""
function rand_conjugate_symplectic_pair(T, n::Integer)
    return block_to_pair(rand_conjugate_symplectic_block(T, n))
end

function rand_conjugate_symplectic_pair(n::Integer)
    return rand_conjugate_symplectic_pair(Complex{Float64}, n)
end

"""
    rand_bogoliubov_block(n::Integer)

Return a random `2n x 2n` Bogoliubov matrix `S`,
`S ∈ Sp(2n, ℂ) ∩ U(n, n) =: Bog(n)`, with block operator order.

"""
function rand_bogoliubov_block(n::Integer)
    return rand_bogoliubov_block(Complex{Float64}, n)
end

function rand_bogoliubov_block(T, n::Integer)
    P = randn(T, n, n)
    # make P symmetric
    P = (P + transpose(P)) / 2

    Q = randn(T, n, n)
    # make Q skew-Hermitian
    Q = (Q - Q') / 2

    # assemble M
    M = [P Q; transpose(Q) -conj(P)]

    # the matrix exponential is numerically more accurate than the Cayley
    # transform here
    Omega = symplectic_form_block(n)
    return exp(Omega * M)

end

"""
    rand_bogoliubov_pair(n::Integer)

Return a random `2n x 2n` Bogoliubov matrix `S`,
`S ∈ Sp(2n, ℂ) ∩ U(n, n) =: Bog(n)`, with pair operator order.

"""
function rand_bogoliubov_pair(n::Integer)
    return rand_bogoliubov_pair(Complex{Float64}, n)
end

function rand_bogoliubov_pair(T, n::Integer)
    return block_to_pair(rand_bogoliubov_block(T, n))
end

"""
    rand_orthogonal_bogoliubov_block(n::Integer)

Return a random `2n x 2n` orthogonal Bogoliubov matrix `S`,
`S ∈ Sp(2n, ℂ) ∩ U(n, n) ∩ U(2n) ≅ U(n)`, with block operator order.

"""
function rand_orthogonal_bogoliubov_block(n::Integer)
    return rand_orthogonal_bogoliubov_block(Complex{Float64}, n)
end

function rand_orthogonal_bogoliubov_block(T, n::Integer)

    Q = randn(T, n, n)
    # make Q skew-Hermitian
    Q = (Q - Q') / 2

    # assemble M
    M = [0*I(n) Q; transpose(Q) 0*I(n)]

    # compute the symplectic matrix using the Cayley transform
    Omega = symplectic_form_block(n)
    return cayley_transform(Omega, M)
end

"""
    rand_orthogonal_bogoliubov_pair(n::Integer)

Return a random `2n x 2n` orthogonal Bogoliubov matrix `S`,
`S ∈ Sp(2n, ℂ) ∩ U(n, n) ∩ U(2n) ≅ U(n)`, with pair operator order.

"""
function rand_orthogonal_bogoliubov_pair(n::Integer)
    return rand_orthogonal_bogoliubov_pair(Complex{Float64}, n)
end

function rand_orthogonal_bogoliubov_pair(T, n::Integer)
    return block_to_pair(rand_orthogonal_bogoliubov_block(T, n))
end


"""
    rand_pseudo_unitary_block(n::Integer)

Return a random `2n x 2n` pseudo-unitary matrix `S`, `S ∈ U(n, n)`, with block
operator order.

"""
function rand_pseudo_unitary_block(n::Integer)
    return rand_pseudo_unitary_block(Complex{Float64}, n)
end

function rand_pseudo_unitary_block(T, n::Integer)
    A = randn(T, 2 * n, 2 * n)
    K = indefinite_hermitian_form_block(n)
    # generate a random skew-Hermitian
    M = (A - A') / 2
    return cayley_transform(K, M)
end

"""
    rand_pseudo_unitary_pair(n::Integer)

Return a random `2n x 2n` pseudo-unitary matrix `S`, `S ∈ U(n, n)`, with
pair operator order.

"""
function rand_pseudo_unitary_pair(n::Integer)
    return rand_pseudo_unitary_pair(Complex{Float64}, n)
end

function rand_pseudo_unitary_pair(T, n::Integer)
    A = randn(T, 2 * n, 2 * n)
    K = indefinite_hermitian_form_pair(n)
    # generate a random skew-Hermitian
    M = (A - A') / 2
    return cayley_transform(K, M)
end

"""
    cayley_transform(Omega, M)

The Cayley transform `(I + Omega*M)*inv(I - Omega*M)`, which takes `Omega*M`
of the Lie algebra of the group preserving the form `Omega` to a member of
that group: a symplectic matrix for the symplectic form and a symmetric `M`,
a pseudo-unitary one for an indefinite Hermitian form and a skew-Hermitian
`M`, as the random matrices of these groups are made.
"""
function cayley_transform(Omega, M)
    n = size(M, 1)
    S = (I(n) + Omega * M) * inv(I(n) - Omega * M)
    return S
end

"""
    rand_cptp_quadrature_block(T, nsys; nenv = nsys, hbar = 1,
        sigma_env = hbar*I(2*nenv))
    rand_cptp_quadrature_block(nsys; nenv = nsys, hbar = 1)

Return a random completely positive trace preserving (CPTP) map `(X, Y)`
of `nsys` modes in the quadrature basis with block operator order: the
system block `X` of a random symplectic matrix of `nsys + nenv` modes and
the noise `Y = B*sigma_env*transpose(B)` its system-environment block `B`
adds from an environment of `nenv` modes with covariance `sigma_env`. The
covariances are in the units of `hbar` of [`is_cptp`](@ref), whose vacuum
is `(hbar/2)*I`; the default `sigma_env = hbar*I` is a thermal state of
half a photon in each mode, and a given `sigma_env` is taken in these units
as it is. `T` is the element type, `Float64` by default.
"""
function rand_cptp_quadrature_block(T, nsys::Integer; nenv::Integer=nsys,
    hbar::Real=1, sigma_env=hbar * I(2 * nenv))
    checkhbar(hbar)
    # start from a pair ordered matrix
    S = rand_symplectic_pair(T, nsys + nenv)
    # and convert each block to the block ordering
    A = pair_to_block(S[1:2*nsys, 1:2*nsys])
    B = pair_to_block(S[1:2*nsys, 2*nsys+1:end])

    X = A
    Y = B * sigma_env * transpose(B)
    return (X=X, Y=Y)
end

function rand_cptp_quadrature_block(nsys::Integer; nenv::Integer=nsys, hbar::Real=1)
    return rand_cptp_quadrature_block(Float64, nsys; nenv=nenv, hbar=hbar)
end

"""
    rand_cptp_quadrature_pair(T, nsys; nenv = nsys, hbar = 1,
        sigma_env = hbar*I(2*nenv))
    rand_cptp_quadrature_pair(nsys; nenv = nsys, hbar = 1)

Return a random CPTP map `(X, Y)` of `nsys` modes in the quadrature basis
with pair operator order, with an environment of `nenv` modes, in the units
of `hbar` of [`is_cptp`](@ref); see [`rand_cptp_quadrature_block`](@ref).
"""
function rand_cptp_quadrature_pair(T, nsys::Integer; nenv::Integer=nsys,
    hbar::Real=1, sigma_env=hbar * I(2 * nenv))
    checkhbar(hbar)
    S = rand_symplectic_pair(T, nsys + nenv)
    A = S[1:2*nsys, 1:2*nsys]
    B = S[1:2*nsys, 2*nsys+1:end]

    X = A
    Y = B * sigma_env * transpose(B)
    return (X=X, Y=Y)
end

function rand_cptp_quadrature_pair(nsys::Integer; nenv::Integer=nsys, hbar::Real=1)
    return rand_cptp_quadrature_pair(Float64, nsys; nenv=nenv, hbar=hbar)
end

"""
    rand_cptp_ladder_pair(T, nsys; nenv = nsys, sigma_env = I(2*nenv))
    rand_cptp_ladder_pair(nsys; nenv = nsys)

Return a random CPTP map `(X, Y)` of `nsys` modes in the ladder basis with
pair operator order: the system block `X` of a random Bogoliubov matrix of
`nsys + nenv` modes and the noise `Y = B*sigma_env*B'` its
system-environment block `B` adds from an environment of `nenv` modes with
covariance `sigma_env`. The covariances are symmetrized, with the vacuum
`I/2` of [`is_cptp_ladder_pair`](@ref); the default `sigma_env = I` is a
thermal state of half a photon in each mode. `T` is the element type,
`Complex{Float64}` by default.
"""
function rand_cptp_ladder_pair(T, nsys::Integer; nenv::Integer=nsys,
    sigma_env=I(2 * nenv))
    S = rand_bogoliubov_pair(T, nsys + nenv)
    A = S[1:2*nsys, 1:2*nsys]
    B = S[1:2*nsys, 2*nsys+1:end]

    X = A
    Y = B * sigma_env * B'
    return (X=X, Y=Y)
end

function rand_cptp_ladder_pair(nsys::Integer; nenv::Integer=nsys)
    return rand_cptp_ladder_pair(Complex{Float64}, nsys; nenv=nenv)
end

"""
    rand_cptp_ladder_block(T, nsys; nenv = nsys, sigma_env = I(2*nenv))
    rand_cptp_ladder_block(nsys; nenv = nsys)

Return a random CPTP map `(X, Y)` of `nsys` modes in the ladder basis with
block operator order, with an environment of `nenv` modes and the
symmetrized covariances of [`rand_cptp_ladder_pair`](@ref).
"""
function rand_cptp_ladder_block(T, nsys::Integer; nenv::Integer=nsys,
    sigma_env=I(2 * nenv))
    # start from a pair ordered matrix
    S = rand_bogoliubov_pair(T, nsys + nenv)
    # and convert each block to the block ordering
    A = pair_to_block(S[1:2*nsys, 1:2*nsys])
    B = pair_to_block(S[1:2*nsys, 2*nsys+1:end])

    X = A
    Y = B * sigma_env * B'
    return (X=X, Y=Y)
end

function rand_cptp_ladder_block(nsys::Integer; nenv::Integer=nsys)
    return rand_cptp_ladder_block(Complex{Float64}, nsys; nenv=nenv)
end




"""
    block_to_pair_perm(n::Int)

Return a vector `p` which permutes the block operator ordering into the
pair operator ordering.

Return a `2n` length vector `p` which permutes the block operator ordering
`r = [x1,...,xn,p1,...,pn]` to the pair operator ordering
`r[p] = (x_1,p_1,...,x_n,p_n)`.

# Examples
```jldoctest
r = [:x1, :x2, :x3, :x4, :p1, :p2, :p3, :p4]
p = JosephsonCircuits.block_to_pair_perm(4)
r[p]

# output
8-element Vector{Symbol}:
 :x1
 :p1
 :x2
 :p2
 :x3
 :p3
 :x4
 :p4
```
"""
function block_to_pair_perm(n::Integer)
    p = Vector{Int}(undef, 2 * n)
    for k in 1:n
        # positions, the odd terms, are from k
        p[2*k-1] = k
        # momenta, the even terms, are from n+k
        p[2*k] = n + k
    end
    return p
end

"""
    pair_to_block_perm(n::Int)

Return a `2n` length vector `p` which permutes the pair operator
ordering `r = [x_1,p_1,...,x_n,p_n]` to the block operator ordering
`r[p] = (x1,...,xn,p1,...,pn)`.

# Examples
```jldoctest
r = [:x1, :p1, :x2, :p2, :x3, :p3, :x4, :p4]
p = JosephsonCircuits.pair_to_block_perm(4)
r[p]

# output
8-element Vector{Symbol}:
 :x1
 :x2
 :x3
 :x4
 :p1
 :p2
 :p3
 :p4
```
"""
function pair_to_block_perm(n::Integer)
    p = Vector{Int}(undef, 2 * n)
    for k in 1:n
        # the first block, position terms, are from the odd pair terms
        p[k] = 2 * k - 1
        # the second block, momentum terms, are from the even pair terms
        p[n+k] = 2 * k
    end
    return p
end


"""
    R_block_to_pair(n::Integer)

Return the `2n x 2n` permutation matrix which takes a vector of `n` modes
in block order to pair order.

# Examples
```jldoctest
julia> JosephsonCircuits.R_block_to_pair(2)
4×4 SparseArrays.SparseMatrixCSC{Int64, Int64} with 4 stored entries:
 1  ⋅  ⋅  ⋅
 ⋅  ⋅  1  ⋅
 ⋅  1  ⋅  ⋅
 ⋅  ⋅  ⋅  1
```
"""
function R_block_to_pair(n::Integer)
    # column pointer has length number of columns + 1
    colptr = [i for i in 1:2*n+1]
    rowval = Vector{Int}(undef, 2 * n)
    for i in 1:n
        rowval[i] = 2 * i - 1
    end
    for i in n+1:2*n
        rowval[i] = 2 * i - 2 * n
    end
    nzval = ones(Int, 2 * n)

    return SparseMatrixCSC(2 * n, 2 * n, colptr, rowval, nzval)
end

"""
    block_to_pair(r::AbstractVector)

# Examples
```jldoctest
r = [:x1, :x2, :x3, :x4, :p1, :p2, :p3, :p4]
JosephsonCircuits.block_to_pair(r)

# output
8-element Vector{Symbol}:
 :x1
 :p1
 :x2
 :p2
 :x3
 :p3
 :x4
 :p4
```
"""
function block_to_pair(r::AbstractVector)
    p = block_to_pair_perm(length(r) ÷ 2)
    return r[p]
end

"""
    block_to_pair2(r::AbstractVector)

Return the vector `r` in block order reordered to pair order, as
[`block_to_pair`](@ref), by a product with the permutation matrix
[`R_block_to_pair`](@ref).
"""
function block_to_pair2(r::AbstractVector)
    R = R_block_to_pair(length(r) ÷ 2)
    return R * r
end

"""
    block_to_pair(S::AbstractMatrix)

Return the matrix `S` with its rows and columns in block order reordered
to pair order.
"""
function block_to_pair(S::AbstractMatrix)
    p1 = block_to_pair_perm(size(S, 1) ÷ 2)
    p2 = block_to_pair_perm(size(S, 2) ÷ 2)
    return S[p1, p2]
end

"""
    block_to_pair2(S::AbstractMatrix)

Return the matrix `S` with its rows and columns in block order reordered
to pair order, as [`block_to_pair`](@ref), by products with the
permutation matrices [`R_block_to_pair`](@ref).
"""
function block_to_pair2(S::AbstractMatrix)
    R1 = R_block_to_pair(size(S, 1) ÷ 2)
    R2 = R_block_to_pair(size(S, 2) ÷ 2)
    return R1 * S * R2'
end

"""
    R_pair_to_block(n::Integer)

Return the `2n x 2n` permutation matrix which takes a vector of `n` modes
in pair order to block order.

# Examples
```jldoctest
julia> JosephsonCircuits.R_pair_to_block(2)
4×4 SparseArrays.SparseMatrixCSC{Int64, Int64} with 4 stored entries:
 1  ⋅  ⋅  ⋅
 ⋅  ⋅  1  ⋅
 ⋅  1  ⋅  ⋅
 ⋅  ⋅  ⋅  1
```
"""
function R_pair_to_block(n::Integer)
    # column pointer has length number of columns + 1
    colptr = [i for i in 1:2*n+1]
    rowval = Vector{Int}(undef, 2 * n)
    j = 1
    for i in 1:2:2*n
        rowval[i] = j
        j += 1
    end
    for i in 2:2:2*n
        rowval[i] = j
        j += 1
    end
    nzval = ones(Int, 2 * n)
    return SparseMatrixCSC(2 * n, 2 * n, colptr, rowval, nzval)
end

"""
    pair_to_block(r::AbstractVector)

# Examples
```jldoctest
r = [:x1, :p1, :x2, :p2, :x3, :p3, :x4, :p4]
JosephsonCircuits.pair_to_block(r)

# output
8-element Vector{Symbol}:
 :x1
 :x2
 :x3
 :x4
 :p1
 :p2
 :p3
 :p4
```
"""
function pair_to_block(r::AbstractVector)
    p = pair_to_block_perm(length(r) ÷ 2)
    return r[p]
end

"""
    pair_to_block2(r::AbstractVector)

Return the vector `r` in pair order reordered to block order, as
[`pair_to_block`](@ref), by a product with the permutation matrix
[`R_pair_to_block`](@ref).
"""
function pair_to_block2(r::AbstractVector)
    R = R_pair_to_block(length(r) ÷ 2)
    return R * r
end

"""
    pair_to_block(S::AbstractMatrix)

Return the matrix `S` with its rows and columns in pair order reordered
to block order.
"""
function pair_to_block(S::AbstractMatrix)
    p1 = pair_to_block_perm(size(S, 1) ÷ 2)
    p2 = pair_to_block_perm(size(S, 2) ÷ 2)
    return S[p1, p2]
end

"""
    pair_to_block2(S::AbstractMatrix)

Return the matrix `S` with its rows and columns in pair order reordered
to block order, as [`pair_to_block`](@ref), by products with the
permutation matrices [`R_pair_to_block`](@ref).
"""
function pair_to_block2(S::AbstractMatrix)
    R1 = R_pair_to_block(size(S, 1) ÷ 2)
    R2 = R_pair_to_block(size(S, 2) ÷ 2)
    return R1 * S * R2'
end

"""
    R_ladder_to_quadrature_pair(n::Integer)

Return the `2n x 2n` unitary matrix which takes the ladder operators of `n`
modes in pair order to their quadratures, `x = (a + adag)/sqrt(2)` and
`p = -im*(a - adag)/sqrt(2)`, in pair order.

# Examples
```
julia> JosephsonCircuits.R_ladder_to_quadrature_pair(1)
2×2 SparseArrays.SparseMatrixCSC{ComplexF64, Int64} with 4 stored entries:
 0.707107+0.0im       0.707107+0.0im
      0.0-0.707107im       0.0+0.707107im
```
"""
function R_ladder_to_quadrature_pair(n::Integer)
    return direct_sum(sparse([1 1; -im im] / sqrt(2)); n=n)
end

"""
    ladder_to_quadrature_pair(r::AbstractVector)

Return the vector `r` of ladder operators in pair order `[a_1, adag_1, ...]`
in the basis of the quadratures in pair order `[x_1, p_1, ...]`, `R*r` with
`R` the matrix of [`R_ladder_to_quadrature_pair`](@ref).
"""
function ladder_to_quadrature_pair(r::AbstractVector)
    R = R_ladder_to_quadrature_pair(length(r) ÷ 2)
    return R * r
end

"""
    ladder_to_quadrature_pair(S::AbstractMatrix)

Return the matrix `S`, which maps ladder operators in pair order
`[a_1, adag_1, ...]` to others, as the map between the quadratures in pair
order `[x_1, p_1, ...]`: `R1*S*R2'` with `R1` and `R2` the matrices of
[`R_ladder_to_quadrature_pair`](@ref) for the rows and the columns of `S`.
"""
function ladder_to_quadrature_pair(S::AbstractMatrix)
    R1 = R_ladder_to_quadrature_pair(size(S, 1) ÷ 2)
    R2 = R_ladder_to_quadrature_pair(size(S, 2) ÷ 2)
    return R1 * S * R2'
end

"""
    R_quadrature_to_ladder_pair(n::Integer)

Return the `2n x 2n` unitary matrix which takes the quadratures of `n`
modes in pair order to their ladder operators, `a = (x + im*p)/sqrt(2)` and
`adag = (x - im*p)/sqrt(2)`, in pair order; the inverse of
[`R_ladder_to_quadrature_pair`](@ref).

# Examples
```jldoctest
julia> JosephsonCircuits.R_quadrature_to_ladder_pair(2)
4×4 SparseArrays.SparseMatrixCSC{ComplexF64, Int64} with 8 stored entries:
 0.707107+0.0im  0.0+0.707107im           ⋅          ⋅
 0.707107+0.0im  0.0-0.707107im           ⋅          ⋅
          ⋅          ⋅           0.707107+0.0im  0.0+0.707107im
          ⋅          ⋅           0.707107+0.0im  0.0-0.707107im
```
"""
function R_quadrature_to_ladder_pair(n::Integer)
    return direct_sum(sparse([1 im; 1 -im] / sqrt(2)); n=n)
end


"""
    quadrature_to_ladder_pair(r::AbstractVector)

Return the vector `r` of quadratures in pair order `[x_1, p_1, ...]` in the
basis of the ladder operators in pair order `[a_1, adag_1, ...]`, `R*r` with
`R` the matrix of [`R_quadrature_to_ladder_pair`](@ref).
"""
function quadrature_to_ladder_pair(r::AbstractVector)
    R = R_quadrature_to_ladder_pair(length(r) ÷ 2)
    return R * r
end

"""
    quadrature_to_ladder_pair(S::AbstractMatrix)

Return the matrix `S`, which maps quadratures in pair order
`[x_1, p_1, ...]` to others, as the map between the ladder operators in pair
order `[a_1, adag_1, ...]`: `R1*S*R2'` with `R1` and `R2` the matrices of
[`R_quadrature_to_ladder_pair`](@ref) for the rows and the columns of `S`.
"""
function quadrature_to_ladder_pair(S::AbstractMatrix)
    R1 = R_quadrature_to_ladder_pair(size(S, 1) ÷ 2)
    R2 = R_quadrature_to_ladder_pair(size(S, 2) ÷ 2)
    return R1 * S * R2'
end

"""
    R_ladder_to_quadrature_block(n::Integer)

The block ordered form of [`R_ladder_to_quadrature_pair`](@ref).

# Examples
```jldoctest
julia> JosephsonCircuits.R_ladder_to_quadrature_block(1)
2×2 Matrix{ComplexF64}:
 0.707107+0.0im       0.707107+0.0im
      0.0-0.707107im       0.0+0.707107im
```
"""
function R_ladder_to_quadrature_block(n::Integer)
    return [I(n)/sqrt(2) I(n)/sqrt(2); -im*I(n)/sqrt(2) im*I(n)/sqrt(2)]
end

"""
    ladder_to_quadrature_block(r::AbstractVector)

Return the vector `r` of ladder operators in block order
`[a_1, ..., adag_1, ...]` in the basis of the quadratures in block order
`[x_1, ..., p_1, ...]`, `R*r` with `R` the matrix of
[`R_ladder_to_quadrature_block`](@ref).
"""
function ladder_to_quadrature_block(r::AbstractVector)
    R = R_ladder_to_quadrature_block(length(r) ÷ 2)
    return R * r
end

"""
    ladder_to_quadrature_block(S::AbstractMatrix)

Return the matrix `S`, which maps ladder operators in block order
`[a_1, ..., adag_1, ...]` to others, as the map between the quadratures in
block order `[x_1, ..., p_1, ...]`: `R1*S*R2'` with `R1` and `R2` the
matrices of [`R_ladder_to_quadrature_block`](@ref) for the rows and the
columns of `S`.
"""
function ladder_to_quadrature_block(S::AbstractMatrix)
    R1 = R_ladder_to_quadrature_block(size(S, 1) ÷ 2)
    R2 = R_ladder_to_quadrature_block(size(S, 2) ÷ 2)
    return R1 * S * R2'
end

"""
    R_quadrature_to_ladder_block(n::Integer)

The block ordered form of [`R_quadrature_to_ladder_pair`](@ref).

# Examples
```jldoctest
julia> JosephsonCircuits.R_quadrature_to_ladder_block(2)
4×4 Matrix{ComplexF64}:
 0.707107+0.0im       0.0+0.0im  0.0+0.707107im  0.0+0.0im
      0.0+0.0im  0.707107+0.0im  0.0+0.0im       0.0+0.707107im
 0.707107+0.0im       0.0+0.0im  0.0-0.707107im  0.0+0.0im
      0.0+0.0im  0.707107+0.0im  0.0+0.0im       0.0-0.707107im
```
"""
function R_quadrature_to_ladder_block(n::Integer)
    return [I(n)/sqrt(2) im*I(n)/sqrt(2); I(n)/sqrt(2) -im*I(n)/sqrt(2)]
end

"""
    quadrature_to_ladder_block(r::AbstractVector)

Return the vector `r` of quadratures in block order `[x_1, ..., p_1, ...]`
in the basis of the ladder operators in block order
`[a_1, ..., adag_1, ...]`, `R*r` with `R` the matrix of
[`R_quadrature_to_ladder_block`](@ref).
"""
function quadrature_to_ladder_block(r::AbstractVector)
    R = R_quadrature_to_ladder_block(length(r) ÷ 2)
    return R * r
end

"""
    quadrature_to_ladder_block(S::AbstractMatrix)

Return the matrix `S`, which maps quadratures in block order
`[x_1, ..., p_1, ...]` to others, as the map between the ladder operators in
block order `[a_1, ..., adag_1, ...]`: `R1*S*R2'` with `R1` and `R2` the
matrices of [`R_quadrature_to_ladder_block`](@ref) for the rows and the
columns of `S`.
"""
function quadrature_to_ladder_block(S::AbstractMatrix)
    R1 = R_quadrature_to_ladder_block(size(S, 1) ÷ 2)
    R2 = R_quadrature_to_ladder_block(size(S, 2) ÷ 2)
    return R1 * S * R2'
end



# The ladder and quadrature forms of a scattering matrix. A scattering
# matrix relates the amplitudes of modes of signed frequencies `w`: a mode
# of positive frequency is an annihilation operator `a` and any other a
# creation operator `a'`, the mode of scattering index `i` having the
# frequency `w[mod(i-1, length(w))+1]`, so that `w` may hold the frequency
# of each mode of a port. The ladder (Bogoliubov) form relates the pairs of
# operators `(a, a')` of the scattering indices, and the quadrature
# (symplectic) form the pairs `(x, p)` with `a = (x + im*p)/sqrt(2)`. The two
# operators of scattering index `i` of an axis of `n` indices are at `2i-1`
# and `2i` in pair order and at `i` and `i+n` in block order; `pairops` and
# `blockops` give their positions, and one kernel per conversion serves
# both orders. The kernels step the index into `w` along each axis,
# `nextmode`, rather than divide for it at every entry.

# the positions of the two operators of scattering index `i` of an axis of
# `n` scattering indices, in pair and in block order
pairops(i, n) = (2 * i - 1, 2 * i)
blockops(i, n) = (i, i + n)

# true when a mode of frequency `wi` is an annihilation operator
isannihilation(wi) = wi > zero(wi)

# the index into `w` of the scattering index after one at index `k`
nextmode(k, Nmodes) = k == Nmodes ? 1 : k + 1

# the positions of the operator an amplitude belongs to and of the one its
# conjugate belongs to: the first and the second of the pair for an
# annihilation operator, the second and the first for a creation operator
operatorpair((i1, i2), annihilation) = annihilation ? (i1, i2) : (i2, i1)

# the sign an imaginary part takes in the quadrature form: that of an
# annihilation operator, and its opposite for a creation operator, whose
# amplitude is a conjugate
operatorsign(annihilation) = annihilation ? 1 : -1

# the element types of the forms of a scattering matrix of element type
# `T`, and of a scattering matrix read back from a form of element type `T`
quadraturetype(T) = real(T)
scatteringtype(T) = Complex{real(typeof(zero(T) / 2))}

# the default relative tolerance of the conversions back from a form: the
# smaller dimension of `S` times the machine epsilon of its element type,
# and zero when an absolute tolerance is given
defaultrtol(S, atol) =
    (min(size(S, 1), size(S, 2)) * eps(real(float(oneunit(eltype(S)))))) * iszero(atol)

# the checks of the kernels: a vector form holds two operators per
# amplitude and whole ports of modes, and a matrix form is twice the size
# of its scattering matrix
function checkvectorform(S_form, S_scattering, w, form)
    checkmodes(w)
    if length(S_form) != 2 * length(S_scattering)
        throw(DimensionMismatch(lazy"The length of the $(form) vector must be double that of the scattering parameter vector."))
    end
    if mod(length(S_scattering), length(w)) != 0
        throw(DimensionMismatch(lazy"Length of scattering vector must be integer multiples of the number of modes."))
    end
    return nothing
end
function checkmatrixform(S_form, S_scattering, w)
    checkmodes(w)
    if size(S_form) != 2 .* size(S_scattering)
        throw(DimensionMismatch(lazy"The size $(size(S_form)) of the ladder or quadrature form must be twice the size $(size(S_scattering)) of the scattering matrix."))
    end
    return nothing
end

function checkmodes(w)
    if isempty(w)
        throw(ArgumentError("The vector of mode frequencies `w` must not be empty."))
    end
    return nothing
end

# refuse a departure `err` of an entry of a form from the structure of the
# form beyond the tolerances; an exact zero is not compared, so that
# symbolic input converts
function checkformerror(err, atol, rtol, normS, form)
    if !iszero(err) && abs(err) > max(atol, rtol * normS)
        error(lazy"Error in $(form) to scattering parameter conversion larger than `atol` and `rtol`.")
    end
    return nothing
end

# --- the kernels, one per conversion, in either order ---------------------

function scattering_to_quadrature!(S_symplectic::AbstractVector,
        S_scattering::AbstractVector, w::AbstractVector, ops)
    checkvectorform(S_symplectic, S_scattering, w, "symplectic")
    n = length(S_scattering)
    k = 1
    @inbounds for i in 1:n
        i1, i2 = ops(i, n)
        s = operatorsign(isannihilation(w[k]))
        k = nextmode(k, length(w))
        S_symplectic[i1] = sqrt(2) * real(S_scattering[i])
        S_symplectic[i2] = s * sqrt(2) * imag(S_scattering[i])
    end
    return S_symplectic
end

function scattering_to_ladder!(S_bogoliubov::AbstractVector,
        S_scattering::AbstractVector, w::AbstractVector, ops)
    checkvectorform(S_bogoliubov, S_scattering, w, "bogoliubov")
    n = length(S_scattering)
    k = 1
    @inbounds for i in 1:n
        p, q = operatorpair(ops(i, n), isannihilation(w[k]))
        k = nextmode(k, length(w))
        S_bogoliubov[p] = S_scattering[i]
        S_bogoliubov[q] = conj(S_scattering[i])
    end
    return S_bogoliubov
end

function scattering_to_ladder!(S_bogoliubov::AbstractMatrix,
        S_scattering::AbstractMatrix, w::AbstractVector, ops)
    checkmatrixform(S_bogoliubov, S_scattering, w)
    n, m = size(S_scattering)
    z = zero(eltype(S_bogoliubov))
    # an entry from the column mode to the row mode takes the positions of
    # the operators of their amplitudes, and its conjugate those of their
    # conjugates; the other two entries of the block are zero
    kj = 1
    @inbounds for j in 1:m
        r, t = operatorpair(ops(j, m), isannihilation(w[kj]))
        kj = nextmode(kj, length(w))
        ki = 1
        for i in 1:n
            p, q = operatorpair(ops(i, n), isannihilation(w[ki]))
            ki = nextmode(ki, length(w))
            Sij = S_scattering[i, j]
            S_bogoliubov[p, r] = Sij
            S_bogoliubov[q, t] = conj(Sij)
            S_bogoliubov[p, t] = z
            S_bogoliubov[q, r] = z
        end
    end
    return S_bogoliubov
end

function ladder_to_scattering!(S_scattering::AbstractMatrix,
        S_bogoliubov::AbstractMatrix, w::AbstractVector, ops; atol::Real = 0,
        rtol::Real = defaultrtol(S_bogoliubov, atol))
    checkmatrixform(S_bogoliubov, S_scattering, w)
    n, m = size(S_scattering)
    normS = norm(S_bogoliubov)
    kj = 1
    @inbounds for j in 1:m
        r, t = operatorpair(ops(j, m), isannihilation(w[kj]))
        kj = nextmode(kj, length(w))
        ki = 1
        for i in 1:n
            p, q = operatorpair(ops(i, n), isannihilation(w[ki]))
            ki = nextmode(ki, length(w))
            # the entry is the average of the amplitude and the conjugate of
            # the conjugate's entry; their difference and the two entries
            # which should be zero are the departure from a ladder form
            direct = S_bogoliubov[p, r]
            conjugate = conj(S_bogoliubov[q, t])
            checkformerror((direct - conjugate) / 2, atol, rtol, normS, "Bogoliubov")
            checkformerror(S_bogoliubov[p, t], atol, rtol, normS, "Bogoliubov")
            checkformerror(S_bogoliubov[q, r], atol, rtol, normS, "Bogoliubov")
            S_scattering[i, j] = (direct + conjugate) / 2
        end
    end
    return S_scattering
end

function scattering_to_quadrature!(S_symplectic::AbstractMatrix,
        S_scattering::AbstractMatrix, w::AbstractVector, ops)
    checkmatrixform(S_symplectic, S_scattering, w)
    n, m = size(S_scattering)
    # an entry takes the real two by two block of multiplication by it, with
    # the imaginary part, and the second row or column, turned for a
    # creation operator
    kj = 1
    @inbounds for j in 1:m
        j1, j2 = ops(j, m)
        sj = operatorsign(isannihilation(w[kj]))
        kj = nextmode(kj, length(w))
        ki = 1
        for i in 1:n
            i1, i2 = ops(i, n)
            si = operatorsign(isannihilation(w[ki]))
            ki = nextmode(ki, length(w))
            Sij = S_scattering[i, j]
            S_symplectic[i1, j1] = real(Sij)
            S_symplectic[i1, j2] = -sj * imag(Sij)
            S_symplectic[i2, j1] = si * imag(Sij)
            S_symplectic[i2, j2] = si * sj * real(Sij)
        end
    end
    return S_symplectic
end

function quadrature_to_scattering!(S_scattering::AbstractMatrix,
        S_symplectic::AbstractMatrix, w::AbstractVector, ops; atol::Real = 0,
        rtol::Real = defaultrtol(S_symplectic, atol))
    checkmatrixform(S_symplectic, S_scattering, w)
    n, m = size(S_scattering)
    normS = norm(S_symplectic)
    kj = 1
    @inbounds for j in 1:m
        j1, j2 = ops(j, m)
        sj = operatorsign(isannihilation(w[kj]))
        kj = nextmode(kj, length(w))
        ki = 1
        for i in 1:n
            i1, i2 = ops(i, n)
            si = operatorsign(isannihilation(w[ki]))
            ki = nextmode(ki, length(w))
            q11 = S_symplectic[i1, j1]
            q12 = S_symplectic[i1, j2]
            q21 = S_symplectic[i2, j1]
            q22 = S_symplectic[i2, j2]
            # the real and the imaginary part each held twice in the block,
            # averaged; the differences are the departure from the block of
            # an entry
            checkformerror((q11 - si * sj * q22 + im * (si * q21 + sj * q12)) / 2,
                atol, rtol, normS, "symplectic")
            S_scattering[i, j] = (q11 + si * sj * q22 + im * (si * q21 - sj * q12)) / 2
        end
    end
    return S_scattering
end

# --- the conversions ---------------------------------------------------

"""
    scattering_to_quadrature_pair(S_scattering::AbstractVector, w)

Return the quadrature vector `[x_1, p_1, ..., x_n, p_n]` of the vector of
amplitudes `S_scattering` of modes of signed frequencies `w`, with
`x = sqrt(2)*real(a)` and `p = sqrt(2)*imag(a)` for the amplitude `a` of an
annihilation operator (a mode of positive frequency) and `p` of the
opposite sign for that of a creation operator. The length of
`S_scattering` must be a multiple of that of `w`.
"""
scattering_to_quadrature_pair(S_scattering::AbstractVector, w) =
    scattering_to_quadrature_pair!(
        zeros(quadraturetype(eltype(S_scattering)), 2 * length(S_scattering)),
        S_scattering, w)

"""
    scattering_to_quadrature_block(S_scattering::AbstractVector, w)

Return the quadrature vector `[x_1, ..., x_n, p_1, ..., p_n]` of the vector
of amplitudes `S_scattering` of modes of signed frequencies `w`; see
[`scattering_to_quadrature_pair`](@ref).
"""
scattering_to_quadrature_block(S_scattering::AbstractVector, w) =
    scattering_to_quadrature_block!(
        zeros(quadraturetype(eltype(S_scattering)), 2 * length(S_scattering)),
        S_scattering, w)

"""
    scattering_to_ladder_pair(S_scattering::AbstractVector, w)

Return the ladder vector `[a_1, a_1', ..., a_n, a_n']` of the vector of
amplitudes `S_scattering` of modes of signed frequencies `w`: each
amplitude and its conjugate, in the order of an annihilation operator
(a mode of positive frequency) and its conjugate, or of a creation
operator's conjugate and itself. The length of `S_scattering` must be a
multiple of that of `w`.
"""
scattering_to_ladder_pair(S_scattering::AbstractVector, w) =
    scattering_to_ladder_pair!(
        zeros(eltype(S_scattering), 2 * length(S_scattering)), S_scattering, w)

"""
    scattering_to_ladder_block(S_scattering::AbstractVector, w)

Return the ladder vector `[a_1, ..., a_n, a_1', ..., a_n']` of the vector of
amplitudes `S_scattering` of modes of signed frequencies `w`; see
[`scattering_to_ladder_pair`](@ref).
"""
scattering_to_ladder_block(S_scattering::AbstractVector, w) =
    scattering_to_ladder_block!(
        zeros(eltype(S_scattering), 2 * length(S_scattering)), S_scattering, w)

"""
    scattering_to_ladder_pair(S_scattering::AbstractMatrix, w)

Return the ladder (Bogoliubov) form, in pair operator order
`ξ = [a_1, a_1', ..., a_n, a_n']`, of the scattering matrix `S_scattering`
between modes of signed frequencies `w`. A mode of positive frequency is
an annihilation operator `a` and any other a creation operator `a'`, the
mode of scattering index `i` having the frequency
`w[mod(i-1, length(w))+1]`. Each entry and its conjugate take two of the
four entries of the two by two block of its row and column, according to
the operators of the two modes, and the other two are zero.
`S_scattering` may be rectangular.

See also [`ladder_to_scattering_pair`](@ref),
[`scattering_to_ladder_block`](@ref) and
[`scattering_to_quadrature_pair`](@ref).
"""
function scattering_to_ladder_pair(S_scattering::AbstractMatrix, w)
    n, m = size(S_scattering)
    return scattering_to_ladder_pair!(zeros(eltype(S_scattering), 2 * n, 2 * m),
        S_scattering, w)
end

"""
    scattering_to_ladder_block(S_scattering::AbstractMatrix, w)

Return the ladder (Bogoliubov) form, in block operator order
`ξ = [a_1, ..., a_n, a_1', ..., a_n']`, of the scattering matrix
`S_scattering` between modes of signed frequencies `w`; see
[`scattering_to_ladder_pair`](@ref).
"""
function scattering_to_ladder_block(S_scattering::AbstractMatrix, w)
    n, m = size(S_scattering)
    return scattering_to_ladder_block!(zeros(eltype(S_scattering), 2 * n, 2 * m),
        S_scattering, w)
end

"""
    scattering_to_quadrature_pair(S_scattering::AbstractMatrix, w)

Return the quadrature (symplectic) form, in pair operator order
`r = [x_1, p_1, ..., x_n, p_n]`, of the scattering matrix `S_scattering`
between modes of signed frequencies `w`, with `a = (x + im*p)/sqrt(2)` for
an annihilation operator `a`; see [`scattering_to_ladder_pair`](@ref) for
the modes. Each entry takes the real two by two block of multiplication by
it, turned for creation operators. The form is real even when
`S_scattering` is complex.
"""
function scattering_to_quadrature_pair(S_scattering::AbstractMatrix, w)
    n, m = size(S_scattering)
    return scattering_to_quadrature_pair!(
        zeros(quadraturetype(eltype(S_scattering)), 2 * n, 2 * m), S_scattering, w)
end

"""
    scattering_to_quadrature_block(S_scattering::AbstractMatrix, w)

Return the quadrature (symplectic) form, in block operator order
`r = [x_1, ..., x_n, p_1, ..., p_n]`, of the scattering matrix
`S_scattering` between modes of signed frequencies `w`; see
[`scattering_to_quadrature_pair`](@ref).
"""
function scattering_to_quadrature_block(S_scattering::AbstractMatrix, w)
    n, m = size(S_scattering)
    return scattering_to_quadrature_block!(
        zeros(quadraturetype(eltype(S_scattering)), 2 * n, 2 * m), S_scattering, w)
end

"""
    scattering_to_quadrature_pair!(S_symplectic, S_scattering, w)

In place version of [`scattering_to_quadrature_pair`](@ref), writing into
`S_symplectic`.
"""
scattering_to_quadrature_pair!(S_symplectic::AbstractArray,
    S_scattering::AbstractArray, w::AbstractVector) =
    scattering_to_quadrature!(S_symplectic, S_scattering, w, pairops)
"""
    scattering_to_quadrature_block!(S_symplectic, S_scattering, w)

In place version of [`scattering_to_quadrature_block`](@ref), writing into
`S_symplectic`.
"""
scattering_to_quadrature_block!(S_symplectic::AbstractArray,
    S_scattering::AbstractArray, w::AbstractVector) =
    scattering_to_quadrature!(S_symplectic, S_scattering, w, blockops)
"""
    scattering_to_ladder_pair!(S_bogoliubov, S_scattering, w)

In place version of [`scattering_to_ladder_pair`](@ref), writing into
`S_bogoliubov`.
"""
scattering_to_ladder_pair!(S_bogoliubov::AbstractArray,
    S_scattering::AbstractArray, w::AbstractVector) =
    scattering_to_ladder!(S_bogoliubov, S_scattering, w, pairops)
"""
    scattering_to_ladder_block!(S_bogoliubov, S_scattering, w)

In place version of [`scattering_to_ladder_block`](@ref), writing into
`S_bogoliubov`.
"""
scattering_to_ladder_block!(S_bogoliubov::AbstractArray,
    S_scattering::AbstractArray, w::AbstractVector) =
    scattering_to_ladder!(S_bogoliubov, S_scattering, w, blockops)

# the scattering matrix read back from the form `S`
scatteringfromform(S::AbstractMatrix) =
    zeros(scatteringtype(eltype(S)), size(S, 1) ÷ 2, size(S, 2) ÷ 2)

"""
    ladder_to_scattering_pair(S_bogoliubov, w; atol = 0, rtol = ...)

Return the scattering matrix whose ladder (Bogoliubov) form in pair
operator order is `S_bogoliubov`, between modes of signed frequencies `w`;
the inverse of [`scattering_to_ladder_pair`](@ref). Each entry is the
average of the two entries of the form holding it and its conjugate; if
their difference, or an entry of the form which should be zero, exceeds
`max(atol, rtol*norm(S_bogoliubov))`, it is an error. The default `rtol`
is the smaller dimension of `S_bogoliubov` times the machine epsilon of
its element type, and zero when `atol` is given. Input which is not of
floating point type, such as symbolic input, is converted with both
tolerances zero.
"""
function ladder_to_scattering_pair(S_bogoliubov::AbstractMatrix{T}, w; atol::Real=0,
    rtol::Real=defaultrtol(S_bogoliubov, atol)) where T<:Union{AbstractFloat,Complex{<:AbstractFloat}}
    return ladder_to_scattering_pair!(scatteringfromform(S_bogoliubov),
        S_bogoliubov, w; atol = atol, rtol = rtol)
end

function ladder_to_scattering_pair(S_bogoliubov::AbstractArray, w)
    return ladder_to_scattering_pair!(scatteringfromform(S_bogoliubov),
        S_bogoliubov, w; atol = 0, rtol = 0)
end

"""
    ladder_to_scattering_block(S_bogoliubov, w; atol = 0, rtol = ...)

Return the scattering matrix whose ladder (Bogoliubov) form in block
operator order is `S_bogoliubov`, between modes of signed frequencies `w`;
the inverse of [`scattering_to_ladder_block`](@ref). See
[`ladder_to_scattering_pair`](@ref) for the tolerances.
"""
function ladder_to_scattering_block(S_bogoliubov::AbstractMatrix{T}, w; atol::Real=0,
    rtol::Real=defaultrtol(S_bogoliubov, atol)) where T<:Union{AbstractFloat,Complex{<:AbstractFloat}}
    return ladder_to_scattering_block!(scatteringfromform(S_bogoliubov),
        S_bogoliubov, w; atol = atol, rtol = rtol)
end

function ladder_to_scattering_block(S_bogoliubov::AbstractArray, w)
    return ladder_to_scattering_block!(scatteringfromform(S_bogoliubov),
        S_bogoliubov, w; atol = 0, rtol = 0)
end

"""
    quadrature_to_scattering_pair(S_symplectic, w; atol = 0, rtol = ...)

Return the scattering matrix whose quadrature (symplectic) form in pair
operator order is `S_symplectic`, between modes of signed frequencies `w`;
the inverse of [`scattering_to_quadrature_pair`](@ref). The real and the
imaginary part of each entry are the averages of the two entries of its
block holding each; if the differences exceed
`max(atol, rtol*norm(S_symplectic))`, it is an error. See
[`ladder_to_scattering_pair`](@ref) for the default tolerances.
"""
function quadrature_to_scattering_pair(S_symplectic::AbstractMatrix{T}, w; atol::Real=0,
    rtol::Real=defaultrtol(S_symplectic, atol)) where T<:Union{AbstractFloat,Complex{<:AbstractFloat}}
    return quadrature_to_scattering_pair!(scatteringfromform(S_symplectic),
        S_symplectic, w; atol = atol, rtol = rtol)
end

function quadrature_to_scattering_pair(S_symplectic::AbstractArray, w)
    return quadrature_to_scattering_pair!(scatteringfromform(S_symplectic),
        S_symplectic, w; atol = 0, rtol = 0)
end

"""
    quadrature_to_scattering_block(S_symplectic, w; atol = 0, rtol = ...)

Return the scattering matrix whose quadrature (symplectic) form in block
operator order is `S_symplectic`, between modes of signed frequencies `w`;
the inverse of [`scattering_to_quadrature_block`](@ref). See
[`quadrature_to_scattering_pair`](@ref) for the tolerances.
"""
function quadrature_to_scattering_block(S_symplectic::AbstractMatrix{T}, w; atol::Real=0,
    rtol::Real=defaultrtol(S_symplectic, atol)) where T<:Union{AbstractFloat,Complex{<:AbstractFloat}}
    return quadrature_to_scattering_block!(scatteringfromform(S_symplectic),
        S_symplectic, w; atol = atol, rtol = rtol)
end

function quadrature_to_scattering_block(S_symplectic::AbstractArray, w)
    return quadrature_to_scattering_block!(scatteringfromform(S_symplectic),
        S_symplectic, w; atol = 0, rtol = 0)
end

"""
    ladder_to_scattering_pair!(S_scattering, S_bogoliubov, w; atol = 0,
        rtol = ...)

In place version of [`ladder_to_scattering_pair`](@ref), writing into
`S_scattering`.
"""
ladder_to_scattering_pair!(S_scattering, S_bogoliubov, w; atol::Real=0,
    rtol::Real=defaultrtol(S_bogoliubov, atol)) =
    ladder_to_scattering!(S_scattering, S_bogoliubov, w, pairops;
        atol = atol, rtol = rtol)
"""
    ladder_to_scattering_block!(S_scattering, S_bogoliubov, w; atol = 0,
        rtol = ...)

In place version of [`ladder_to_scattering_block`](@ref), writing into
`S_scattering`.
"""
ladder_to_scattering_block!(S_scattering, S_bogoliubov, w; atol::Real=0,
    rtol::Real=defaultrtol(S_bogoliubov, atol)) =
    ladder_to_scattering!(S_scattering, S_bogoliubov, w, blockops;
        atol = atol, rtol = rtol)
"""
    quadrature_to_scattering_pair!(S_scattering, S_symplectic, w; atol = 0,
        rtol = ...)

In place version of [`quadrature_to_scattering_pair`](@ref), writing into
`S_scattering`.
"""
quadrature_to_scattering_pair!(S_scattering, S_symplectic, w; atol::Real=0,
    rtol::Real=defaultrtol(S_symplectic, atol)) =
    quadrature_to_scattering!(S_scattering, S_symplectic, w, pairops;
        atol = atol, rtol = rtol)
"""
    quadrature_to_scattering_block!(S_scattering, S_symplectic, w; atol = 0,
        rtol = ...)

In place version of [`quadrature_to_scattering_block`](@ref), writing into
`S_scattering`.
"""
quadrature_to_scattering_block!(S_scattering, S_symplectic, w; atol::Real=0,
    rtol::Real=defaultrtol(S_symplectic, atol)) =
    quadrature_to_scattering!(S_scattering, S_symplectic, w, blockops;
        atol = atol, rtol = rtol)


"""
    ports_modes_to_modes_ports_perm(Nports,Nmodes)

Return a permutation vector that converts one axis of a scattering matrix with
`Nports` ports and `Nmodes` modes from (ports, modes) ordering to
(modes, ports) ordering. In (ports, modes) ordering the port index runs
fastest, as in an array of size `(Nports, Nmodes)`: port `p` of mode `m` is
at `p + (m-1)*Nports`. In (modes, ports) ordering the mode index runs
fastest, which is the order of the solvers' scattering matrices: mode `m`
of port `p` is at `m + (p-1)*Nmodes`. Entry `k` of the permutation is the
index in (ports, modes) ordering of the entry that goes to index `k`.

# Examples
For 2 ports and 4 modes the (ports, modes) ordering is
`[(p1,m1), (p2,m1), (p1,m2), (p2,m2), (p1,m3), (p2,m3), (p1,m4), (p2,m4)]`
and the (modes, ports) ordering is
`[(p1,m1), (p1,m2), (p1,m3), (p1,m4), (p2,m1), (p2,m2), (p2,m3), (p2,m4)]`:
```jldoctest
julia> p = JosephsonCircuits.ports_modes_to_modes_ports_perm(2,4)
8-element Vector{Int64}:
 1
 3
 5
 7
 2
 4
 6
 8
```
"""
function ports_modes_to_modes_ports_perm(Nports, Nmodes)
    # the index in (ports, modes) ordering of each port and mode, read in
    # (modes, ports) ordering
    return vec(permutedims(reshape(1:Nports*Nmodes, Nports, Nmodes)))
end

"""
    modes_ports_to_ports_modes_perm(Nports,Nmodes)

Return a permutation vector that converts one axis of a scattering matrix with
`Nports` ports and `Nmodes` modes from (modes, ports) ordering to
(ports, modes) ordering, the inverse of
[`ports_modes_to_modes_ports_perm`](@ref), where the orderings are
described.

# Examples
For 2 ports and 4 modes the (modes, ports) ordering is
`[(p1,m1), (p1,m2), (p1,m3), (p1,m4), (p2,m1), (p2,m2), (p2,m3), (p2,m4)]`
and the (ports, modes) ordering is
`[(p1,m1), (p2,m1), (p1,m2), (p2,m2), (p1,m3), (p2,m3), (p1,m4), (p2,m4)]`:
```jldoctest
julia> p = JosephsonCircuits.modes_ports_to_ports_modes_perm(2,4)
8-element Vector{Int64}:
 1
 5
 2
 6
 3
 7
 4
 8
```
"""
function modes_ports_to_ports_modes_perm(Nports, Nmodes)
    # the index in (modes, ports) ordering of each port and mode, read in
    # (ports, modes) ordering
    return vec(permutedims(reshape(1:Nports*Nmodes, Nmodes, Nports)))
end

"""
    scattering_to_pair_perm(p0::Vector{Int})

Return the permutation of an axis of the pair ordered form of a scattering
matrix, in which scattering index `i` has its two operators at `2i-1` and
`2i`, that moves the scattering indices as the permutation `p0` of the
scattering matrix does, each pair together. For example the permutation
from (ports, modes) to (modes, ports) ordering of 2 ports and 4 modes,
[`ports_modes_to_modes_ports_perm`](@ref), becomes:

# Examples
```jldoctest
julia> JosephsonCircuits.scattering_to_pair_perm(JosephsonCircuits.ports_modes_to_modes_ports_perm(2,4))
16-element Vector{Int64}:
  1
  2
  5
  6
  9
 10
 13
 14
  3
  4
  7
  8
 11
 12
 15
 16
```
"""
function scattering_to_pair_perm(p0::Vector{Int})

    # the two operators of scattering index p0[i] go to 2i-1 and 2i
    p = Vector{Int}(undef, 2 * length(p0))
    for i in eachindex(p0)
        p[2*i-1] = 2 * p0[i] - 1
        p[2*i] = 2 * p0[i]
    end

    return p
end

"""
    scattering_to_block_perm(p0::Vector{Int})

Return the permutation of an axis of the block ordered form of a scattering
matrix, in which scattering index `i` has its two operators at `i` and
`i+n` for `n` scattering indices, that moves the scattering indices as the
permutation `p0` of the scattering matrix does, in each block alike. For
example the permutation from (ports, modes) to (modes, ports) ordering of 2
ports and 4 modes, [`ports_modes_to_modes_ports_perm`](@ref), becomes:

# Examples
```jldoctest
julia> JosephsonCircuits.scattering_to_block_perm(JosephsonCircuits.ports_modes_to_modes_ports_perm(2,4))
16-element Vector{Int64}:
  1
  3
  5
  7
  2
  4
  6
  8
  9
 11
 13
 15
 10
 12
 14
 16
```
"""
function scattering_to_block_perm(p0::Vector{Int})

    # the two operators of scattering index p0[i] go to i and i+n
    p = Vector{Int}(undef, 2 * length(p0))
    for i in eachindex(p0)
        p[i] = p0[i]
        p[i+length(p0)] = p0[i] + length(p0)
    end

    return p
end


# The reorderings between (ports, modes) and (modes, ports) index order,
# for scattering matrices and for their pair and block symplectic forms.
# The pair and block forms permute within the symplectic structure: the
# scattering permutation is built first and then adapted to each.

# the number of ports of an axis of `n` scattering indices with `Nmodes`
# modes each
function axisports(n::Integer, Nmodes::Integer)
    if Nmodes < 1 || mod(n, Nmodes) != 0
        throw(DimensionMismatch(lazy"The number of scattering indices $(n) of an axis must be a multiple of the number of modes $(Nmodes)."))
    end
    return n ÷ Nmodes
end

# the number of scattering indices of an axis of length `n` of a pair or
# block form, two operators per scattering index
function axisindices(n::Integer)
    if isodd(n)
        throw(DimensionMismatch(lazy"The length $(n) of an axis of a pair or block form must be even."))
    end
    return n ÷ 2
end

# the permutation `perm(Nports, Nmodes)` of an axis of length `n` of a
# scattering matrix, of its pair form, and of its block form
scatteringaxisperm(perm, n, Nmodes) = perm(axisports(n, Nmodes), Nmodes)
pairaxisperm(perm, n, Nmodes) =
    scattering_to_pair_perm(perm(axisports(axisindices(n), Nmodes), Nmodes))
blockaxisperm(perm, n, Nmodes) =
    scattering_to_block_perm(perm(axisports(axisindices(n), Nmodes), Nmodes))

"""
    modes_ports_to_ports_modes_scattering(S::AbstractMatrix,Nmodes::Int)

Reorder both axes of the scattering matrix `S` of `Nmodes` modes per port
from (modes, ports) to (ports, modes) ordering; see
[`modes_ports_to_ports_modes_perm`](@ref). The length of each axis must be
a multiple of `Nmodes`.

# Examples
```jldoctest
S = [:S11 :S12 :S13 :S14 :S15 :S16 :S17 :S18; :S21 :S22 :S23 :S24 :S25 :S26 :S27 :S28; :S31 :S32 :S33 :S34 :S35 :S36 :S37 :S38; :S41 :S42 :S43 :S44 :S45 :S46 :S47 :S48; :S51 :S52 :S53 :S54 :S55 :S56 :S57 :S58; :S61 :S62 :S63 :S64 :S65 :S66 :S67 :S68; :S71 :S72 :S73 :S74 :S75 :S76 :S77 :S78; :S81 :S82 :S83 :S84 :S85 :S86 :S87 :S88]
JosephsonCircuits.modes_ports_to_ports_modes_scattering(S,2)

# output
8×8 Matrix{Symbol}:
 :S11  :S13  :S15  :S17  :S12  :S14  :S16  :S18
 :S31  :S33  :S35  :S37  :S32  :S34  :S36  :S38
 :S51  :S53  :S55  :S57  :S52  :S54  :S56  :S58
 :S71  :S73  :S75  :S77  :S72  :S74  :S76  :S78
 :S21  :S23  :S25  :S27  :S22  :S24  :S26  :S28
 :S41  :S43  :S45  :S47  :S42  :S44  :S46  :S48
 :S61  :S63  :S65  :S67  :S62  :S64  :S66  :S68
 :S81  :S83  :S85  :S87  :S82  :S84  :S86  :S88
```
"""
function modes_ports_to_ports_modes_scattering(S::AbstractMatrix, Nmodes::Int)
    return S[scatteringaxisperm(modes_ports_to_ports_modes_perm, size(S, 1), Nmodes),
        scatteringaxisperm(modes_ports_to_ports_modes_perm, size(S, 2), Nmodes)]
end

"""
    ports_modes_to_modes_ports_scattering(S::AbstractMatrix,Nmodes::Int)

Reorder both axes of the scattering matrix `S` of `Nmodes` modes per port
from (ports, modes) to (modes, ports) ordering; see
[`ports_modes_to_modes_ports_perm`](@ref). The length of each axis must be
a multiple of `Nmodes`.

# Examples
```jldoctest
S = [:S11 :S12 :S13 :S14 :S15 :S16 :S17 :S18; :S21 :S22 :S23 :S24 :S25 :S26 :S27 :S28; :S31 :S32 :S33 :S34 :S35 :S36 :S37 :S38; :S41 :S42 :S43 :S44 :S45 :S46 :S47 :S48; :S51 :S52 :S53 :S54 :S55 :S56 :S57 :S58; :S61 :S62 :S63 :S64 :S65 :S66 :S67 :S68; :S71 :S72 :S73 :S74 :S75 :S76 :S77 :S78; :S81 :S82 :S83 :S84 :S85 :S86 :S87 :S88]
JosephsonCircuits.ports_modes_to_modes_ports_scattering(S,2)

# output
8×8 Matrix{Symbol}:
 :S11  :S15  :S12  :S16  :S13  :S17  :S14  :S18
 :S51  :S55  :S52  :S56  :S53  :S57  :S54  :S58
 :S21  :S25  :S22  :S26  :S23  :S27  :S24  :S28
 :S61  :S65  :S62  :S66  :S63  :S67  :S64  :S68
 :S31  :S35  :S32  :S36  :S33  :S37  :S34  :S38
 :S71  :S75  :S72  :S76  :S73  :S77  :S74  :S78
 :S41  :S45  :S42  :S46  :S43  :S47  :S44  :S48
 :S81  :S85  :S82  :S86  :S83  :S87  :S84  :S88
```
"""
function ports_modes_to_modes_ports_scattering(S::AbstractMatrix, Nmodes::Int)
    return S[scatteringaxisperm(ports_modes_to_modes_ports_perm, size(S, 1), Nmodes),
        scatteringaxisperm(ports_modes_to_modes_ports_perm, size(S, 2), Nmodes)]
end

"""
    ports_modes_to_modes_ports_pair(S::AbstractMatrix,Nmodes::Int)

Reorder both axes of the pair ordered ladder or quadrature form `S` of a
scattering matrix of `Nmodes` modes per port from (ports, modes) to
(modes, ports) ordering, moving the two operators of each scattering index
together; see [`ports_modes_to_modes_ports_perm`](@ref). Half the length
of each axis must be a multiple of `Nmodes`.

# Examples
```jldoctest
S = [:S11 :S12 :S13 :S14 :S15 :S16 :S17 :S18; :S21 :S22 :S23 :S24 :S25 :S26 :S27 :S28; :S31 :S32 :S33 :S34 :S35 :S36 :S37 :S38; :S41 :S42 :S43 :S44 :S45 :S46 :S47 :S48; :S51 :S52 :S53 :S54 :S55 :S56 :S57 :S58; :S61 :S62 :S63 :S64 :S65 :S66 :S67 :S68; :S71 :S72 :S73 :S74 :S75 :S76 :S77 :S78; :S81 :S82 :S83 :S84 :S85 :S86 :S87 :S88]
JosephsonCircuits.ports_modes_to_modes_ports_pair(S,2)

# output
8×8 Matrix{Symbol}:
 :S11  :S12  :S15  :S16  :S13  :S14  :S17  :S18
 :S21  :S22  :S25  :S26  :S23  :S24  :S27  :S28
 :S51  :S52  :S55  :S56  :S53  :S54  :S57  :S58
 :S61  :S62  :S65  :S66  :S63  :S64  :S67  :S68
 :S31  :S32  :S35  :S36  :S33  :S34  :S37  :S38
 :S41  :S42  :S45  :S46  :S43  :S44  :S47  :S48
 :S71  :S72  :S75  :S76  :S73  :S74  :S77  :S78
 :S81  :S82  :S85  :S86  :S83  :S84  :S87  :S88
```
"""
function ports_modes_to_modes_ports_pair(S::AbstractMatrix, Nmodes::Int)
    return S[pairaxisperm(ports_modes_to_modes_ports_perm, size(S, 1), Nmodes),
        pairaxisperm(ports_modes_to_modes_ports_perm, size(S, 2), Nmodes)]
end

"""
    modes_ports_to_ports_modes_pair(S::AbstractMatrix,Nmodes::Int)

Reorder both axes of the pair ordered ladder or quadrature form `S` of a
scattering matrix of `Nmodes` modes per port from (modes, ports) to
(ports, modes) ordering, moving the two operators of each scattering index
together; see [`modes_ports_to_ports_modes_perm`](@ref). Half the length
of each axis must be a multiple of `Nmodes`.

# Examples
```jldoctest
S = [:S11 :S12 :S13 :S14 :S15 :S16 :S17 :S18; :S21 :S22 :S23 :S24 :S25 :S26 :S27 :S28; :S31 :S32 :S33 :S34 :S35 :S36 :S37 :S38; :S41 :S42 :S43 :S44 :S45 :S46 :S47 :S48; :S51 :S52 :S53 :S54 :S55 :S56 :S57 :S58; :S61 :S62 :S63 :S64 :S65 :S66 :S67 :S68; :S71 :S72 :S73 :S74 :S75 :S76 :S77 :S78; :S81 :S82 :S83 :S84 :S85 :S86 :S87 :S88]
JosephsonCircuits.modes_ports_to_ports_modes_pair(S,2)

# output
8×8 Matrix{Symbol}:
 :S11  :S12  :S15  :S16  :S13  :S14  :S17  :S18
 :S21  :S22  :S25  :S26  :S23  :S24  :S27  :S28
 :S51  :S52  :S55  :S56  :S53  :S54  :S57  :S58
 :S61  :S62  :S65  :S66  :S63  :S64  :S67  :S68
 :S31  :S32  :S35  :S36  :S33  :S34  :S37  :S38
 :S41  :S42  :S45  :S46  :S43  :S44  :S47  :S48
 :S71  :S72  :S75  :S76  :S73  :S74  :S77  :S78
 :S81  :S82  :S85  :S86  :S83  :S84  :S87  :S88
```
"""
function modes_ports_to_ports_modes_pair(S::AbstractMatrix, Nmodes::Int)
    return S[pairaxisperm(modes_ports_to_ports_modes_perm, size(S, 1), Nmodes),
        pairaxisperm(modes_ports_to_ports_modes_perm, size(S, 2), Nmodes)]
end

"""
    ports_modes_to_modes_ports_block(S::AbstractMatrix,Nmodes::Int)

Reorder both axes of the block ordered ladder or quadrature form `S` of a
scattering matrix of `Nmodes` modes per port from (ports, modes) to
(modes, ports) ordering, each of the two blocks of an axis alike; see
[`ports_modes_to_modes_ports_perm`](@ref). Half the length of each axis
must be a multiple of `Nmodes`.

# Examples
```jldoctest
S = [:S11 :S12 :S13 :S14 :S15 :S16 :S17 :S18; :S21 :S22 :S23 :S24 :S25 :S26 :S27 :S28; :S31 :S32 :S33 :S34 :S35 :S36 :S37 :S38; :S41 :S42 :S43 :S44 :S45 :S46 :S47 :S48; :S51 :S52 :S53 :S54 :S55 :S56 :S57 :S58; :S61 :S62 :S63 :S64 :S65 :S66 :S67 :S68; :S71 :S72 :S73 :S74 :S75 :S76 :S77 :S78; :S81 :S82 :S83 :S84 :S85 :S86 :S87 :S88]
JosephsonCircuits.ports_modes_to_modes_ports_block(S,2)

# output
8×8 Matrix{Symbol}:
 :S11  :S13  :S12  :S14  :S15  :S17  :S16  :S18
 :S31  :S33  :S32  :S34  :S35  :S37  :S36  :S38
 :S21  :S23  :S22  :S24  :S25  :S27  :S26  :S28
 :S41  :S43  :S42  :S44  :S45  :S47  :S46  :S48
 :S51  :S53  :S52  :S54  :S55  :S57  :S56  :S58
 :S71  :S73  :S72  :S74  :S75  :S77  :S76  :S78
 :S61  :S63  :S62  :S64  :S65  :S67  :S66  :S68
 :S81  :S83  :S82  :S84  :S85  :S87  :S86  :S88
```
"""
function ports_modes_to_modes_ports_block(S::AbstractMatrix, Nmodes::Int)
    return S[blockaxisperm(ports_modes_to_modes_ports_perm, size(S, 1), Nmodes),
        blockaxisperm(ports_modes_to_modes_ports_perm, size(S, 2), Nmodes)]
end

"""
    modes_ports_to_ports_modes_block(S::AbstractMatrix,Nmodes::Int)

Reorder both axes of the block ordered ladder or quadrature form `S` of a
scattering matrix of `Nmodes` modes per port from (modes, ports) to
(ports, modes) ordering, each of the two blocks of an axis alike; see
[`modes_ports_to_ports_modes_perm`](@ref). Half the length of each axis
must be a multiple of `Nmodes`.

# Examples
```jldoctest
S = [:S11 :S12 :S13 :S14 :S15 :S16 :S17 :S18; :S21 :S22 :S23 :S24 :S25 :S26 :S27 :S28; :S31 :S32 :S33 :S34 :S35 :S36 :S37 :S38; :S41 :S42 :S43 :S44 :S45 :S46 :S47 :S48; :S51 :S52 :S53 :S54 :S55 :S56 :S57 :S58; :S61 :S62 :S63 :S64 :S65 :S66 :S67 :S68; :S71 :S72 :S73 :S74 :S75 :S76 :S77 :S78; :S81 :S82 :S83 :S84 :S85 :S86 :S87 :S88]
JosephsonCircuits.modes_ports_to_ports_modes_block(S,2)

# output
8×8 Matrix{Symbol}:
 :S11  :S13  :S12  :S14  :S15  :S17  :S16  :S18
 :S31  :S33  :S32  :S34  :S35  :S37  :S36  :S38
 :S21  :S23  :S22  :S24  :S25  :S27  :S26  :S28
 :S41  :S43  :S42  :S44  :S45  :S47  :S46  :S48
 :S51  :S53  :S52  :S54  :S55  :S57  :S56  :S58
 :S71  :S73  :S72  :S74  :S75  :S77  :S76  :S78
 :S61  :S63  :S62  :S64  :S65  :S67  :S66  :S68
 :S81  :S83  :S82  :S84  :S85  :S87  :S86  :S88
```
"""
function modes_ports_to_ports_modes_block(S::AbstractMatrix, Nmodes::Int)
    return S[blockaxisperm(modes_ports_to_ports_modes_perm, size(S, 1), Nmodes),
        blockaxisperm(modes_ports_to_ports_modes_perm, size(S, 2), Nmodes)]
end

"""
    optimum_eigenvalue_angle(values; target_angle = pi)

For eigenvalues `values` on or near the unit circle, the angle of the
midpoint of the widest empty arc between them, and the angle of the
rotation which moves that midpoint to `target_angle` (by default `pi`, so
that the rotated eigenvalues lie as far as they can from the branch cut
of the logarithm and the square root). The arcs are the counterclockwise
gaps between the sorted angles of the eigenvalues, the last one running
from the largest angle around to the smallest, so that eigenvalues
clustered on a short arc leave the rest of the circle as the widest gap.
After https://github.com/XanaduAI/thewalrus/pull/403.

# Examples
Eigenvalues clustered about `1.5` leave the widest gap about `1.5 - pi`,
and the rotation by `-1.5` moves it to `pi`:
```jldoctest
julia> optimum_angle, optimum_rotation = JosephsonCircuits.optimum_eigenvalue_angle(cis.([1.4, 1.5, 1.6]));

julia> round(optimum_angle, digits = 4), round(optimum_rotation, digits = 4)
(-1.6416, -1.5)
```
"""
function optimum_eigenvalue_angle(values; target_angle=pi)
    # the angles of the eigenvalues in increasing order
    angles = sort!(angle.(values))
    n = length(angles)

    # the widest counterclockwise gap between consecutive angles, starting
    # with the one from the largest angle around to the smallest
    max_arc_width = angles[1] + 2 * pi - angles[n]
    arc_midpoint = angles[n] + max_arc_width / 2
    for i in 1:n-1
        arc_width = angles[i+1] - angles[i]
        if arc_width > max_arc_width
            max_arc_width = arc_width
            arc_midpoint = angles[i] + arc_width / 2
        end
    end
    return angle(cis(arc_midpoint)), angle(cis(target_angle - arc_midpoint))
end

"""
    polar(A)

Return a positive semi-definite matrix `P` and unitary matrix `U` such that
`A = P U`. `A` is a square real or complex matrix.

If `A` is symplectic, both `P` and `U` are symplectic.

This definition is different from wikipedia
https://en.wikipedia.org/wiki/Polar_decomposition

# References
[1] M. Houde, W. McCutcheon, and N. Quesada, “Matrix decompositions in Quantum
Optics: Takagi/Autonne, Bloch-Messiah/Euler, Iwasawa, and Williamson,” Can.
J. Phys., vol. 102, no. 10, pp. 497–507, Oct. 2024, doi: 10.1139/cjp-2024-0070.
[2] https://en.wikipedia.org/wiki/Polar_decomposition#Relation_to_the_SVD
"""
function polar(A::AbstractMatrix)

    F = svd(A)
    P = F.U * Diagonal(F.S) * F.U'

    U = F.U * F.Vt

    return (P=P, U=U)
end

"""
    williamson_pair(M::AbstractMatrix{<:Real}; atol = 0, rtol = ...)

For a symmetric positive semi-definite matrix `M`, return a vector of values `d`
and a real symplectic matrix `S` such that `M = S Diagonal(d) S^T`. `S` is
symplectic with respect to the pair ordered symplectic form `Ω`.

Such a decomposition exists for a positive definite `M`, and for a
positive semi-definite one on whose range the symplectic form is
nondegenerate, which requires an even rank. Otherwise, as for
`Diagonal([1, 0])` of rank one or `Diagonal([1, 0, 1, 0])`, whose range
holds the positions alone, an `ArgumentError` is thrown.

The values `d` are unique but the matrix `S` is not. `S` is computed from
a Schur decomposition [2] rather than an eigendecomposition, which is
robust to degenerate values.

`M` must be symmetric to the tolerances `atol` and `rtol` of `isapprox`
(the default `rtol` is the square root of the machine epsilon of its
element type, and zero when `atol` is given), and its symmetric part
`(M + transpose(M))/2` is decomposed, so that a covariance which rounding
left slightly asymmetric, such as `S*Diagonal(d)*transpose(S)`, is
accepted.

# References
[1] M. Idel, S. Soto Gaona, and M. M. Wolf, “Perturbation bounds for
Williamson’s symplectic normal form,” Linear Algebra and its Applications,
vol. 525, pp. 45–58, Jul. 2017, doi: 10.1016/j.laa.2017.03.013.
[2] M. Houde, W. McCutcheon, and N. Quesada, “Matrix decompositions in Quantum
Optics: Takagi/Autonne, Bloch-Messiah/Euler, Iwasawa, and Williamson,” Can.
J. Phys., vol. 102, no. 10, pp. 497–507, Oct. 2024, doi: 10.1139/cjp-2024-0070.
"""
function williamson_pair(M::AbstractMatrix{<:Real}; atol::Real = 0,
        rtol::Real = approxrtol(eltype(M), atol))
    n = size(M, 1) ÷ 2
    Omega = symplectic_form_pair(n)
    d, S = _williamson(Omega, symmetricpart(M, atol, rtol))
    return d, S
end

"""
    williamson_block(M::AbstractMatrix{<:Real}; atol = 0, rtol = ...)

The block ordered form of [`williamson_pair`](@ref): `d` is in block order
and `S` is symplectic with respect to the block ordered symplectic form.
"""
function williamson_block(M::AbstractMatrix{<:Real}; atol::Real = 0,
        rtol::Real = approxrtol(eltype(M), atol))
    n = size(M, 1) ÷ 2
    # computed directly in block ordering rather than through
    # `williamson_pair(block_to_pair(M))`
    Omega = symplectic_form_block(n)
    d, S = _williamson(Omega, symmetricpart(M, atol, rtol))
    return pair_to_block(d), S * R_block_to_pair(n)
end

# the symmetric part (M + transpose(M))/2 of a matrix `M` which is
# symmetric to the tolerances `atol` and `rtol`, exactly symmetric, as the
# factorizations of a symmetric matrix require
function symmetricpart(M, atol, rtol)
    if !isapprox(M, transpose(M); atol = atol, rtol = rtol)
        error(lazy"M must be symmetric.")
    end
    return (M + transpose(M)) / 2
end

# the 2x2 blocks of the real Schur form `T` of a real skew-symmetric matrix,
# the first `nblocks` of them: the half difference `a` of the off-diagonal
# entries of each, and the scaling of its two columns, `1/sqrt(|a|)` with
# `inverse` and `sqrt(|a|)` without, the second column's times the sign of
# `a`. The scaling with `inverse` takes the block to the symplectic form of
# one mode in pair order; without, it is the factor which takes that form
# to the block.
function schurblockscaling(T, nblocks; inverse::Bool)
    a = [(T[i, i+1] - T[i+1, i]) / 2 for i in 1:2:2*nblocks]
    scaling = similar(a, 2 * nblocks)
    for (k, ak) in enumerate(a)
        s = inverse ? inv(sqrt(abs(ak))) : sqrt(abs(ak))
        scaling[2k-1] = s
        scaling[2k] = sign(ak) * s
    end
    return a, scaling
end

function _williamson(Omega, M::AbstractMatrix{<:Real})
    n = size(M, 1) ÷ 2

    # a pivoted Cholesky factorization, which is robust against small
    # negative eigenvalues from rounding; done in a function because the
    # sparse and dense factorizations have different syntax
    L, rankL = cholesky_williamson(M)

    K = transpose(L) * Omega * L

    r = rankL ÷ 2

    # K is skew symmetric, so its real Schur form has 2x2 blocks which are
    # either zero or have opposite values on the off diagonal; the
    # symplectic normal form follows from it
    F = schur(K)
    T = Matrix(F.T)
    Z = Matrix(F.Z)

    # the values and the scalings of the blocks, in pair ordering; the
    # caller permutes to block ordering
    a, phi = schurblockscaling(T, r; inverse = true)
    # a zero block: the symplectic form is degenerate on the range of M, as
    # on the span of x_1 and x_2 alone, and no symplectic matrix brings M to
    # a diagonal form
    if any(ak -> !(abs(ak) > 0), a)
        throw(ArgumentError("The symplectic form is degenerate on the range of M, so it has no symplectic normal form."))
    end
    d = repeat(abs.(a); inner = 2)

    S1 = L * Z * Diagonal(phi)

    # a full rank matrix is done; a rank deficient one is completed with
    # the symplectic complement
    if r == n
        d, S = d, S1
    else
        S2 = symplectic_complement(Omega, S1)
        d, S = vcat(d, zeros(eltype(d), 2 * (n - r))), hcat(S1, S2)
    end
    # a range on which the symplectic form is degenerate to rounding, such
    # as that of the positions alone rotated, gives blocks at rounding level
    # and can give a matrix S which is not symplectic; S takes the pair
    # ordered form of its columns to Omega
    if !isapprox(S * symplectic_form_pair(n) * transpose(S), Omega)
        throw(ArgumentError("The symplectic form is degenerate on the range of M to working precision, so it has no symplectic normal form."))
    end
    return (d=d, S=S)
end

function cholesky_williamson(M::AbstractArray)
    # a pivoted Cholesky factorization with error checking off, robust
    # against small negative eigenvalues from rounding, of the matrix taken
    # as symmetric, in its dense form: the pivoted factorization is
    # LAPACK's, which a structured matrix such as a `Diagonal` does not
    # reach on every Julia version
    C = cholesky(Symmetric(convert(Matrix, M)), RowMaximum(); check=false)

    # the pivoted factorization sometimes reports a rank one larger than
    # the true (even) rank; an odd rank is reduced by one and the reduced
    # factorization is then checked against the matrix below
    if isodd(C.rank)
        rankL = C.rank - 1
    else
        rankL = C.rank
    end

    L = C.L[invperm(C.p), 1:rankL]

    # if the factorization reports failure, or the rank was reduced, check
    # that L*L' still reproduces the matrix; if only the factor of the
    # reported, odd, rank does, the matrix is positive semi-definite of odd
    # rank, which no symplectic normal form has
    if !issuccess(C) || isodd(C.rank)
        if !isapprox(M, L * L')
            Lodd = C.L[invperm(C.p), 1:C.rank]
            if isodd(C.rank) && isapprox(M, Lodd * Lodd')
                throw(ArgumentError(lazy"The rank $(C.rank) of M is odd, so it has no symplectic normal form."))
            end
            error(lazy"Cholesky factorization has failed. Input matrix is not positive semi-definite.")
        end
    end

    return L, rankL
end

# the sparse factorization of a positive definite matrix, of full rank; a
# matrix it does not factorize, semidefinite or indefinite, takes the rank
# revealing dense path, so that both storages accept the same matrices
function cholesky_williamson(M::SparseMatrixCSC)
    C = cholesky(M; check = false)
    issuccess(C) || return cholesky_williamson(Matrix(M))
    return Matrix(sparse(C.L)[invperm(C.p), :]), size(M, 1)
end

function symplectic_complement(Omega, S1)

    n = size(S1, 1) ÷ 2
    r = size(S1, 2) ÷ 2

    # the symplectic complement of range(S1) = range(M) is ker(S1' Ω),
    # 2n by (2n-r) with orthonormal columns, for the symplectic form Ω of
    # the order of the input
    S2 = nullspace(Matrix(transpose(S1) * Omega))
    # (2n-r) × (2n-r), skew, full-rank
    K2 = transpose(S2) * Omega * S2
    F2 = schur(K2)
    T2 = Matrix(F2.T)
    Z2 = Matrix(F2.Z)
    _, phi2 = schurblockscaling(T2, n - r; inverse = true)
    # 2n × (2n-r)
    return S2 * Z2 * Diagonal(phi2)
end



"""
    symplectic_normal_form_pair(A::AbstractMatrix{<:Real})

For a real skew-symmetric matrix `A` return `Q` such that `A = Q Ω Q^T` where
`Q` is an invertible matrix and `Ω` is the pair symplectic form. If `A`
is singular, the decomposition works, but `Q` is no longer invertible.
"""
function symplectic_normal_form_pair(A::AbstractMatrix{<:Real})

    # check if the matrix A is skew-symmetric
    if !isapprox(A, -transpose(A))
        error(lazy"A must be skew-symmetric.")
    end

    n = size(A, 1)

    if isodd(n)
        error(lazy"A must have even dimensions for a symplectic normal form.")
    end

    # the symplectic normal form from the real Schur decomposition
    F = schur(A)
    return schurnormalform(Matrix(F.T), Matrix(F.Z))
end

# The factor `Q` of `Z T Z^T = Q Ω Q^T`, where `T` and `Z` are the real Schur
# form of a skew-symmetric matrix and `Ω` is the pair symplectic form. The
# 2x2 blocks of `T` hold the eigenvalue pairs ±ia; a singular matrix's also
# has 1x1 zero blocks, which can fall between them. `Q` takes the columns
# of `Z` with each 2x2 block's two together, then the 1x1 blocks' in pairs,
# each pair scaled as `schurblockscaling` scales its block.
function schurnormalform(T, Z)
    n = size(T, 1)
    order, singles = Int[], Int[]
    i = 1
    while i <= n
        if i < n && !iszero(T[i+1, i])
            push!(order, i, i + 1)
            i += 2
        else
            push!(singles, i)
            i += 1
        end
    end
    append!(order, singles)
    _, d = schurblockscaling(T[order, order], n ÷ 2; inverse = false)
    return Z[:, order] * Diagonal(d)
end

"""
    symplectic_normal_form_block(A::AbstractMatrix{<:Real})

The block ordered form of [`symplectic_normal_form_pair`](@ref): `A` is in
block order and `Ω` the block ordered symplectic form.
"""
function symplectic_normal_form_block(A::AbstractMatrix{<:Real})
    Q = symplectic_normal_form_pair(block_to_pair(A))
    return pair_to_block(Q)
end

"""
    inv_symplectic_pair(S)

Return the inverse of the symplectic matrix `S` computed from
`transpose(Ω)*transpose(S)*Ω` where `Ω` is the symplectic form for pair
operator ordering.

# Examples
```jldoctest
julia> S = JosephsonCircuits.rand_symplectic_pair(2);isapprox(inv(S),JosephsonCircuits.inv_symplectic_pair(S))
true
```
"""
function inv_symplectic_pair(S)
    Omega = symplectic_form_pair(size(S, 2) ÷ 2)
    return transpose(Omega) * transpose(S) * Omega
end

"""
    inv_symplectic_block(S)

Return the inverse of the symplectic matrix `S` computed from
`transpose(Ω)*transpose(S)*Ω` where `Ω` is the symplectic form for block
operator ordering.

# Examples
```jldoctest
julia> S = JosephsonCircuits.rand_symplectic_block(2);isapprox(inv(S),JosephsonCircuits.inv_symplectic_block(S))
true
```
"""
function inv_symplectic_block(S)
    Omega = symplectic_form_block(size(S, 2) ÷ 2)
    return transpose(Omega) * transpose(S) * Omega
end

"""
    inv_bogoliubov_pair(S)

Return the inverse of the Bogoliubov matrix `S` computed from
`transpose(Ω)*transpose(S)*Ω` where `Ω` is the symplectic form for pair
operator ordering.

# Examples
```jldoctest
julia> S = JosephsonCircuits.rand_bogoliubov_pair(2);isapprox(inv(S),JosephsonCircuits.inv_bogoliubov_pair(S))
true
```
"""
function inv_bogoliubov_pair(S)
    return inv_symplectic_pair(S)
end

"""
    inv_bogoliubov_block(S)

Return the inverse of the Bogoliubov matrix `S` computed from
`transpose(Ω)*transpose(S)*Ω` where `Ω` is the symplectic form for block
operator ordering.

# Examples
```jldoctest
julia> S = JosephsonCircuits.rand_bogoliubov_block(2);isapprox(inv(S),JosephsonCircuits.inv_bogoliubov_block(S))
true
```
"""
function inv_bogoliubov_block(S)
    return inv_symplectic_block(S)
end


"""
    autonne_takagi(M::AbstractMatrix; atol = 0, rtol = ...)

Return a vector `Λ` and a unitary matrix `W` for a symmetric complex input
matrix `M` such that `M == W*Diagonal(Λ)*transpose(W)` where `M` satisfies
`M = transpose(M)`. Note that if `M` complex this means `M` is not Hermitian.
`M` must be symmetric to the tolerances `atol` and `rtol` of `isapprox`
(the default `rtol` is the square root of the machine epsilon of its
element type, and zero when `atol` is given); its symmetric part
`(M + transpose(M))/2` is factorized.

They are returned as the named tuple `(Λ, W)`, with `Λ` the singular
values of `M` in decreasing order, as for a real `M`. `W` is `U*sqrt(Z)`
for the singular value
decomposition `M = U*Diagonal(Λ)*V'` and the unitary `Z = U'*conj(V)`,
whose square root is taken with its eigenvalues rotated away from the
branch cut by [`optimum_eigenvalue_angle`](@ref).

# References
[1] A. M. Chebotarev and A. E. Teretenkov, “Singular value decomposition for
the Takagi factorization of symmetric matrices,” Applied Mathematics and
Computation, vol. 234, pp. 380–384, May 2014, doi: 10.1016/j.amc.2014.01.170.
[2] M. Houde, W. McCutcheon, and N. Quesada, “Matrix decompositions in Quantum
Optics: Takagi/Autonne, Bloch-Messiah/Euler, Iwasawa, and Williamson,” Can.
J. Phys., vol. 102, no. 10, pp. 497–507, Oct. 2024, doi: 10.1139/cjp-2024-0070.
[3] https://github.com/XanaduAI/thewalrus/pull/403
"""
function autonne_takagi(M::AbstractMatrix; atol::Real = 0,
        rtol::Real = approxrtol(eltype(M), atol))

    # `svd` has no method for a complex `Symmetric` matrix in Julia 1.10
    # and 1.11, so the plain matrix is factorized
    F = svd(symmetricpart(Matrix(M), atol, rtol))
    # an empty matrix has nothing to rotate
    isempty(F.S) && return (Λ = F.S, W = F.U)
    # this is how Chebotarev and Teretenkov 2014 define Z
    Z = F.U' * transpose(F.Vt)

    # Z is unitary, so normal, and its complex Schur factor is diagonal up
    # to rounding: the square root is that of its eigenvalues, rotated
    # away from the branch cut and back
    E = schur(Z)
    optimum_angle, optimum_rotation = optimum_eigenvalue_angle(E.values)
    shift = exp(im * optimum_rotation)
    invshifto2 = exp(-im * optimum_rotation / 2)
    Zsqrt = invshifto2 * E.vectors * Diagonal(sqrt.(shift .* E.values)) * E.vectors'
    W = F.U * Zsqrt
    return (Λ = F.S, W = W)
end


"""
    autonne_takagi(M::AbstractMatrix{<:Real}; atol = 0, rtol = ...)

Return a vector `Λ` and a unitary matrix `W` for input matrix `M` such that
`M == W*Diagonal(Λ)*transpose(W)` where `M` is a symmetric real matrix
`M = transpose(M)`, to the tolerances `atol` and `rtol` as for a complex
`M`; its symmetric part is factorized.

They are returned as the named tuple `(Λ, W)`, as for a complex `M`, with
`Λ` the magnitudes of the eigenvalues of `M` in decreasing order and each
column of `W` the eigenvector times `1` or `im`, the square root of the
sign of its eigenvalue (`1` for a zero eigenvalue).

"""
function autonne_takagi(M::AbstractMatrix{<:Real}; atol::Real = 0,
        rtol::Real = approxrtol(eltype(M), atol))

    F = eigen(Symmetric(symmetricpart(M, atol, rtol)); sortby = λ -> -abs(λ))

    # M = V*Diagonal(λ)*transpose(V) with V real orthogonal, so each column
    # of V times the square root of the sign of its eigenvalue, 1 or im,
    # gives W*Diagonal(abs.(λ))*transpose(W) = M; a zero eigenvalue takes 1
    # like a positive one
    T = eltype(F.values)
    phases = [λ < zero(λ) ? Complex(zero(T), one(T)) : Complex(one(T), zero(T))
        for λ in F.values]
    return (Λ = abs.(F.values), W = F.vectors * Diagonal(phases))

end


"""
    bloch_messiah_block(S::AbstractMatrix{<:Real})

Return the Bloch-Messiah (Euler) decomposition  `O`, `D`, `Q` of the
symplectic matrix `S = O*Diagonal(D)*Q` where `O` and `Q` are
orthogonal-symplectic matrices and `Diagonal(D)` is a symplectic-diagonal and
positive definite matrix. This is also called the symplectic singular value
decomposition (SVD). The matrices are symplectic with respect to the block
symplectic form `Ω`.

This decomposition is unique up to permutations and or degeneracies of the
Takagi-Autonne singular values.

The singular values are the same as the regular SVD, but the ordering of the
singular values and order/signs of the factors are different, in order to make
them orthogonal-symplectic.

# References
[1] G. Cariolaro and G. Pierobon, “Reexamination of Bloch-Messiah reduction,”
Phys. Rev. A, vol. 93, no. 6, p. 062115, Jun. 2016,
doi: 10.1103/PhysRevA.93.062115.
[2] G. Cariolaro and G. Pierobon, “Bloch-Messiah reduction of Gaussian
unitaries by Takagi factorization,” Phys. Rev. A, vol. 94, no. 6, p. 062109,
Dec. 2016, doi: 10.1103/PhysRevA.94.062109.
[3] M. Houde, W. McCutcheon, and N. Quesada, “Matrix decompositions in Quantum
Optics: Takagi/Autonne, Bloch-Messiah/Euler, Iwasawa, and Williamson,” Can.
J. Phys., vol. 102, no. 10, pp. 497–507, Oct. 2024, doi: 10.1139/cjp-2024-0070.
"""
function bloch_messiah_block(S::AbstractMatrix{<:Real})

    # get the size
    n = size(S, 1) ÷ 2

    # test if S is symplectic
    if !is_symplectic_block(S)
        error(lazy"A must be symplectic.")
    end

    # the polar decomposition S = P*Y
    P, Y = polar(S)

    # partition the symplectic matrix P
    # P = [A B; B^T C]
    A = view(P, 1:n, 1:n)
    B = view(P, 1:n, n+1:2*n)
    Bt = view(P, n+1:2*n, 1:n)
    C = view(P, n+1:2*n, n+1:2*n)

    # M = 1/2*(A - C + im*(B + B^T))
    M = 1 / 2 * (A .- C .+ im * B .+ im * Bt)

    # perform Takagi-Autonne decomposition
    # M = W*Λ*W^T
    Λ, W = autonne_takagi(Symmetric(M))

    # Γ = Λ + sqrt(I + Λ^2) and D = Γ ⊕ 1/Γ
    Γ = Λ .+ sqrt.(1 .+ Λ .^ 2)
    D = vcat(Γ, inv.(Γ))

    # O = [Re(W) -Im(W); Im(W) Re(W)]
    O = [real(W) -imag(W); imag(W) real(W)]

    # Q = O^T * Y
    Q = O' * Y

    return (O=O, D=D, Q=Q)

end

"""
    bloch_messiah_pair(S::AbstractMatrix{<:Real})

The pair ordered form of [`bloch_messiah_block`](@ref): `S`, the factors
and `D` are in pair order, and the factors symplectic with respect to the
pair ordered symplectic form.
"""
function bloch_messiah_pair(S::AbstractMatrix{<:Real})
    F = bloch_messiah_block(pair_to_block(S))
    return (O=block_to_pair(F.O), D=block_to_pair(F.D),
        Q=block_to_pair(F.Q))
end

"""
    pre_iwasawa_block(S::AbstractMatrix)

Return the pre-Iwasawa decomposition `E`, `D`, `F` of the symplectic matrix
`S = [S11 S12; S21 S22]` in block operator order, `S = E*D*F`. With
`A0 = sqrt(S11*transpose(S11) + S12*transpose(S12))`, `D = [A0 0; 0 inv(A0)]`
is block diagonal, `E = [I 0; C0*inv(A0) I]` with
`C0 = (S21*transpose(S11) + S22*transpose(S12))*inv(A0)` is lower block
triangular with identity diagonal blocks, and `F = [X Y; -Y X]` with
`X = inv(A0)*S11` and `Y = inv(A0)*S12` is symplectic, and orthogonal for
a real `S`.

# References
[1] M. Houde, W. McCutcheon, and N. Quesada, “Matrix decompositions in Quantum
Optics: Takagi/Autonne, Bloch-Messiah/Euler, Iwasawa, and Williamson,” Can.
J. Phys., vol. 102, no. 10, pp. 497–507, Oct. 2024, doi: 10.1139/cjp-2024-0070.
[2] Arvind, B. Dutta, N. Mukunda, and R. Simon, “The real symplectic groups in
quantum mechanics and optics,” Pramana - J Phys, vol. 45, no. 6, pp. 471–497,
Dec. 1995, doi: 10.1007/BF02848172.
"""
function pre_iwasawa_block(S::AbstractMatrix)

    # get the size
    n = size(S, 1) ÷ 2

    # test if S is symplectic
    if !is_symplectic_block(S)
        error(lazy"A must be symplectic.")
    end

    # partition the symplectic matrix S = [S11 S12; S21 S22]
    S11 = view(S, 1:n, 1:n)
    S12 = view(S, 1:n, n+1:2*n)
    S21 = view(S, n+1:2*n, 1:n)
    S22 = view(S, n+1:2*n, n+1:2*n)

    # the first block row of S has full rank, so for a real S the matrix
    # under the root is symmetric positive definite, whose principal root
    # is that of its eigendecomposition; for a complex S it is complex
    # symmetric, and takes the general root
    G = S11 * transpose(S11) + S12 * transpose(S12)
    A0 = eltype(S) <: Real ? Matrix(sqrt(Symmetric(G))) : sqrt(G)
    invA0 = inv(A0)
    C0 = (S21 * transpose(S11) + S22 * transpose(S12)) * invA0
    X = invA0 * S11
    Y = invA0 * S12

    # the factors in the element type of the root: E the identity with
    # C0*inv(A0) below its diagonal, and D block diagonal
    T = eltype(A0)
    E = Matrix{T}(I, 2*n, 2*n)
    mul!(view(E, n+1:2*n, 1:n), C0, invA0)
    D = zeros(T, 2*n, 2*n)
    D[1:n, 1:n] .= A0
    D[n+1:2*n, n+1:2*n] .= invA0
    F = [X Y; -Y X]

    return (E=E, D=D, F=F)
end

"""
    pre_iwasawa_pair(S::AbstractMatrix)

The pair ordered form of [`pre_iwasawa_block`](@ref): `S` and the factors
are in pair order.
"""
function pre_iwasawa_pair(S::AbstractMatrix)
    F = pre_iwasawa_block(pair_to_block(S))
    return (E=block_to_pair(F.E), D=block_to_pair(F.D),
        F=block_to_pair(F.F))
end

"""
    iwasawa_block(S::AbstractMatrix)

Return the Iwasawa (KAN) decomposition `K`, `A`, `N` of the symplectic matrix
`S`. `K` is a unitary symplectic matrix (maximal compact), `A` is a diagonal
symplectic matrix (Abelian), and `N` is a block upper triangular symplectic
matrix whose first diagonal block is unit upper triangular (nilpotent). The
matrices are symplectic with respect to the block symplectic form `Ω`. This
decomposition is unique.

# References
[1] Arvind, B. Dutta, N. Mukunda, and R. Simon, “The real symplectic groups in
quantum mechanics and optics,” Pramana - J Phys, vol. 45, no. 6, pp. 471–497,
Dec. 1995, doi: 10.1007/BF02848172.
[2] M. Benzi and N. Razouk, “On the Iwasawa decomposition of a symplectic
matrix,” Applied Mathematics Letters, vol. 20, no. 3, pp. 260–265, Mar. 2007,
doi: 10.1016/j.aml.2006.04.004.
[3] M. Houde, W. McCutcheon, and N. Quesada, “Matrix decompositions in Quantum
Optics: Takagi/Autonne, Bloch-Messiah/Euler, Iwasawa, and Williamson,” Can.
J. Phys., vol. 102, no. 10, pp. 497–507, Oct. 2024, doi: 10.1139/cjp-2024-0070.
"""
function iwasawa_block(S::AbstractMatrix)
    # algorithm 2.3 from [2] based on QR factorization

    # get the size
    n = size(S, 1) ÷ 2

    # test if S is symplectic
    if !is_symplectic_block(S)
        error(lazy"A must be symplectic.")
    end

    # partition the symplectic matrix S
    # S = [S11 S12; S21 S22]
    S11 = view(S, 1:n, 1:n)
    S12 = view(S, 1:n, n+1:2*n)
    S21 = view(S, n+1:2*n, 1:n)
    S22 = view(S, n+1:2*n, n+1:2*n)

    S1 = [S11; S21]
    F = qr(S1)
    R11 = Matrix(F.R)
    Q = Matrix(F.Q)

    # define some views
    Q11 = view(Q, 1:n, 1:n)
    Q21 = view(Q, n+1:2*n, 1:n)

    # factor upper triangular matrix R11 = R
    # as R11 = H*U where H is diagonal and U is unit upper triangular
    H = Diagonal(diag(R11))
    U = inv(H) * R11

    # D as a vector; abs2 so that complex matrices work
    D = abs2.(diag(R11))
    Dsqrt = sqrt.(D)
    Dinvsqrt = inv.(Dsqrt)

    # define K, A, N
    A = Diagonal(vcat(Dsqrt, Dinvsqrt))

    K11 = Q11 * H * Diagonal(Dinvsqrt)
    K12 = -Q21 * H * Diagonal(Dinvsqrt)
    # with conjugates so that complex matrices work
    K = [K11 conj(K12); -K12 conj(K11)]


    # adjoint rather than transpose so that complex matrices work
    N = A \ (K' * [S12; S22])

    N12 = view(N, 1:n, 1:n)
    N22 = view(N, n+1:2*n, 1:n)

    N = [U N12; 0*I(n) N22]

    return (K=K, A=A, N=N)
end

"""
    iwasawa_pair(S::AbstractMatrix)

The pair ordered form of [`iwasawa_block`](@ref): `S` and the factors are
in pair order.
"""
function iwasawa_pair(S::AbstractMatrix)
    F = iwasawa_block(pair_to_block(S))
    return (K=block_to_pair(F.K), A=block_to_pair(F.A),
        N=block_to_pair(F.N))
end

"""
    iwasawa_bogoliubov_pair(S::AbstractMatrix)

The Iwasawa decomposition [`iwasawa_block`](@ref) of the Bogoliubov matrix
`S` of the ladder operators in pair order, through its quadrature form:
the factors are Bogoliubov matrices in pair order.
"""
function iwasawa_bogoliubov_pair(S::AbstractMatrix)
    F = iwasawa_block(ladder_to_quadrature_block(pair_to_block(S)))
    return (
        K=quadrature_to_ladder_pair(block_to_pair(F.K)),
        A=quadrature_to_ladder_pair(block_to_pair(F.A)),
        N=quadrature_to_ladder_pair(block_to_pair(F.N)),
    )
end

"""
    iwasawa_bogoliubov_block(S::AbstractMatrix)

The Iwasawa decomposition [`iwasawa_block`](@ref) of the Bogoliubov matrix
`S` of the ladder operators in block order, through its quadrature form:
the factors are Bogoliubov matrices in block order.
"""
function iwasawa_bogoliubov_block(S::AbstractMatrix)
    F = iwasawa_block(ladder_to_quadrature_block(S))
    return (
        K=quadrature_to_ladder_block(F.K),
        A=quadrature_to_ladder_block(F.A),
        N=quadrature_to_ladder_block(F.N),
    )
end




"""
    B_from_X_Y_quadrature(Omega, X::AbstractMatrix{<:Real},
        Y::AbstractMatrix{<:Real}; hbar = 1, atol = 0, rtol = ...)

Return the real `2n x 4n` block `B`, from an environment of `2n` modes in
block operator order to the `n` modes of the system, which realizes the
completely positive trace preserving (CPTP) map of the quadrature
transformation `X` and noise `Y` with the environment in the vacuum:
`X*Omega*transpose(X) + B*ΩE*transpose(B) = Omega`, with `ΩE` the block
ordered symplectic form of the environment, and `Y = (hbar/2)*B*transpose(B)`.
`Omega` is the symplectic form of the system and `hbar` the scale of the
covariances of [`is_cptp`](@ref).

`B` is a square root of `(2/hbar)*K`, with
`K = Y + im*(hbar/2)*(Omega - X*Omega*transpose(X))` the matrix of
[`is_cptp`](@ref), which is Hermitian and positive semi-definite for a
CPTP map, and so `Y`, its real part. `K` is judged as `is_cptp` judges
it, to the same tolerance `tol = max(atol, rtol*s)` of the scale
`s = norm(Y) + (hbar/2)*norm(Omega)*(1 + opnorm(X)^2)`: an eigenvalue down
to `-tol` is rounding and taken as zero, and a map `is_cptp` refuses
throws an `ArgumentError`.
"""
function B_from_X_Y_quadrature(Omega::AbstractMatrix, X::AbstractMatrix{<:Real},
    Y::AbstractMatrix{<:Real}; hbar::Real = 1, atol::Real = 0,
    rtol::Real = approxrtol(promote_type(eltype(X), eltype(Y)), atol))
    K, tol = cptpcondition(Omega, X, Y, hbar, atol, rtol)
    if norm(K - K') > tol
        throw(ArgumentError(lazy"The map is not completely positive and trace preserving: `Y` is not symmetric to the tolerance $(tol)."))
    end

    vals, vecs = eigen(Hermitian((K + K') / 2))

    # an eigenvalue below zero by more than the tolerance is a map which is
    # not CPTP; one within it is rounding
    if !isempty(vals) && vals[1] < -tol
        throw(ArgumentError(lazy"The map is not completely positive and trace preserving: `Y + im*(hbar/2)*(Omega - X*Omega*X')` has the eigenvalue $(vals[1]), below zero by more than the tolerance $(tol)."))
    end
    clamp!(vals, 0, Inf)

    # the noise in the units of a vacuum of I, so that Y = (hbar/2)*B*B'
    F = vecs * Diagonal(sqrt.((2 / hbar) .* vals))

    # this is specific to the block form
    B = [imag(F) real(F)]

    return B
end

"""
    B_from_X_Y_quadrature_block(X::AbstractMatrix{<:Real},
        Y::AbstractMatrix{<:Real}; hbar = 1, atol = 0, rtol = ...)

[`B_from_X_Y_quadrature`](@ref) for the quadrature transformation `X` and
noise `Y` in block operator order: the block `B` of a symplectic matrix
`S = [X B; C D]` with the environment in the vacuum, `Y = (hbar/2)*B*B'`
in the units of `hbar` of [`is_cptp`](@ref).
"""
function B_from_X_Y_quadrature_block(X::AbstractMatrix{<:Real},
    Y::AbstractMatrix{<:Real}; kwargs...)

    n = size(X, 1) ÷ 2
    Omega = symplectic_form_block(n)
    B = B_from_X_Y_quadrature(Omega, X, Y; kwargs...)
    return B
end

"""
    B_from_X_Y_quadrature_pair(X::AbstractMatrix{<:Real},
        Y::AbstractMatrix{<:Real}; hbar = 1, atol = 0, rtol = ...)

[`B_from_X_Y_quadrature_block`](@ref) in pair operator order, the columns
of `B` too, in the units of `hbar` of [`is_cptp`](@ref).
"""
function B_from_X_Y_quadrature_pair(X::AbstractMatrix{<:Real},
    Y::AbstractMatrix{<:Real}; kwargs...)

    n = size(X, 1) ÷ 2
    Omega = symplectic_form_pair(n)
    B = B_from_X_Y_quadrature(Omega, X, Y; kwargs...)

    # permute the columns of B to the pair form
    p = block_to_pair_perm(size(B, 2) ÷ 2)
    return B[:, p]
end

"""
    X_Y_to_symplectic_pair(X::AbstractMatrix{<:Real}, Y::AbstractMatrix{<:Real};
        hbar = 1, atol = 0, rtol = ...)

Return a symplectic matrix `S` of `3n` modes, in pair operator order, whose
restriction to the first `n` modes with an environment of `2n` modes in
the vacuum is the completely positive trace preserving (CPTP) map of the
quadrature transformation `X` and noise `Y` of `n` modes: `X` is the upper
left `2n x 2n` block of `S`, and its upper right block `B` adds the noise
`Y = (hbar/2)*B*transpose(B)` from the vacuum `(hbar/2)*I`, in the
convention of [`is_cptp`](@ref). A map which is not CPTP to the
tolerances `atol` and `rtol` of `is_cptp` throws an `ArgumentError`
([`B_from_X_Y_quadrature`](@ref)).

"""
function X_Y_to_symplectic_pair(X::AbstractMatrix{<:Real}, Y::AbstractMatrix{<:Real};
    hbar::Real = 1, atol::Real = 0,
    rtol::Real = approxrtol(promote_type(eltype(X), eltype(Y)), atol))

    # compute B from X and Y
    B = B_from_X_Y_quadrature_pair(X, Y; hbar, atol, rtol)

    return A_B_to_symplectic_pair(X, B)
end


"""
    X_Y_to_symplectic_block(X::AbstractMatrix{<:Real}, Y::AbstractMatrix{<:Real};
        hbar = 1, atol = 0, rtol = ...)

The block ordered form of [`X_Y_to_symplectic_pair`](@ref) for `X` and `Y`
in block order, in the units of `hbar` of [`is_cptp`](@ref) and to its
tolerances: `S` is in the block order of all `3n` modes, so `X` is the
submatrix of the rows and columns `[1:n; 3n+1:4n]` of `S`.

"""
function X_Y_to_symplectic_block(X::AbstractMatrix{<:Real}, Y::AbstractMatrix{<:Real};
    hbar::Real = 1, atol::Real = 0,
    rtol::Real = approxrtol(promote_type(eltype(X), eltype(Y)), atol))

    # compute B from X and Y
    B = B_from_X_Y_quadrature_block(X, Y; hbar, atol, rtol)

    return pair_to_block(A_B_to_symplectic_pair(block_to_pair(X), block_to_pair(B)))
end

"""
    X_Y_to_bogoliubov_pair(X::AbstractMatrix, Y::AbstractMatrix; atol = 0,
        rtol = ...)

The ladder form of [`X_Y_to_symplectic_pair`](@ref): return a Bogoliubov
matrix of `3n` modes, in pair operator order, which realizes the CPTP map
of the ladder transformation `X` and noise `Y` of `n` modes with an
environment in the vacuum. `Y` is a symmetrized covariance with the vacuum
`I/2` ([`is_cptp_ladder_pair`](@ref)), so the block `B` of the Bogoliubov
matrix from the environment to the system adds `Y = B*B'/2`. A map which
is not CPTP to the tolerances `atol` and `rtol` of `is_cptp_ladder_pair`
throws an `ArgumentError`.

"""
function X_Y_to_bogoliubov_pair(X::AbstractMatrix, Y::AbstractMatrix; atol::Real = 0,
    rtol::Real = approxrtol(promote_type(eltype(X), eltype(Y)), atol))

    # R*M*R', with R unitary, is both the similarity which converts the
    # map X and the congruence which converts the covariance Y; R involves
    # no hbar, so the symmetrized covariance, vacuum I/2, becomes the
    # quadrature covariance at hbar = 1
    X_quadrature = real(ladder_to_quadrature_pair(X))
    Y_quadrature = real(ladder_to_quadrature_pair(Y))

    # compute B from X and Y
    B = B_from_X_Y_quadrature_pair(X_quadrature, Y_quadrature; hbar = 1, atol, rtol)

    return quadrature_to_ladder_pair(A_B_to_symplectic_pair(X_quadrature, B))
end


"""
    X_Y_to_bogoliubov_block(X::AbstractMatrix, Y::AbstractMatrix; atol = 0,
        rtol = ...)

The ladder form of [`X_Y_to_symplectic_block`](@ref), with `X`, `Y` and
the Bogoliubov matrix returned in block operator order, and the symmetrized
covariances and the tolerances of [`X_Y_to_bogoliubov_pair`](@ref).

"""
function X_Y_to_bogoliubov_block(X::AbstractMatrix, Y::AbstractMatrix; atol::Real = 0,
    rtol::Real = approxrtol(promote_type(eltype(X), eltype(Y)), atol))

    # as in X_Y_to_bogoliubov_pair, the conversion takes the covariance Y
    # to the quadrature covariance at hbar = 1
    X_quadrature = real(ladder_to_quadrature_block(X))
    Y_quadrature = real(ladder_to_quadrature_block(Y))

    # compute B from X and Y
    B = B_from_X_Y_quadrature_block(X_quadrature, Y_quadrature; hbar = 1, atol, rtol)

    return quadrature_to_ladder_block(pair_to_block(A_B_to_symplectic_pair(block_to_pair(X_quadrature), block_to_pair(B))))
end

"""
    halmos_dilation(S; atol = 0, rtol = ...)

Return the Halmos dilation of the passive lossy scattering parameter matrix
`S`. This converts a passive lossy scattering parameter matrix into a lossless
(unitary) scattering parameter matrix, `[S sqrt(I - S*S'); sqrt(I - S'*S) -S']`,
whose ports after the first `size(S, 1)` carry the loss of `S`. An `n` by
`m` matrix `S` dilates to an `n + m` by `n + m` one, twice the number of
ports for a square `S`.

`S` is passive when none of its singular values exceeds one. A singular
value above one by no more than `max(atol, rtol)`, as rounding leaves those
of a lossless `S`, is taken as one; beyond that `S` is refused. The default
`rtol` is that of `isapprox` and [`is_unitary`](@ref), the square root of
the machine epsilon of the element type of `S`, and zero when `atol` is
given.

# Examples
```jldoctest
julia> JosephsonCircuits.is_unitary(JosephsonCircuits.halmos_dilation([0.1 0;0 0.1]))
true
```

# References
[1] P. L. Robinson, “Julia operators and Halmos dilations,” Mar. 25, 2018,
    arXiv:1803.09329. doi: 10.48550/arXiv.1803.09329.
[2] B. Sz.-Nagy, C. Foias, H. Bercovici, and L. Kérchy, Harmonic Analysis of
    Operators on Hilbert Space. New York, NY: Springer, 2010.
    doi: 10.1007/978-1-4419-6094-8.
[3] P. R. Halmos, “Normal dilations and extensions of operators,” Summa
    Brasiliensis Mathematicae, vol. II, no. VI, pp. 125–134, Dec. 1950.
[4] J. J. Schäffer, “On Unitary Dilations of Contractions,” Proceedings of the
    American Mathematical Society, vol. 6, no. 2, pp. 322–322, 1955,
    doi: 10.2307/2032368.
[5] B. Szőkefalvi-Nagy, “Sur les contractions de l’espace de Hilbert,”
    ACTA SCIENTIARUM MATHEMATICARUM, vol. 15, pp. 87–92, 1954.
"""
function halmos_dilation(S; atol::Real = 0, rtol::Real = approxrtol(eltype(S), atol))
    n, m = size(S)
    k = min(n, m)

    # the dilation
    # U = [W 0; 0 V]*[Σ sqrt(I-ΣΣ'); sqrt(I-Σ'Σ) -Σ']*[V' 0; 0 W']
    # from the full singular value decomposition S = W Σ V', Σ n by m with
    # the singular values on its diagonal, which equals
    # [S sqrt(I - S S'); sqrt(I - S' S) -S']
    F = svd(S; full = true)

    σmax = isempty(F.S) ? zero(eltype(F.S)) : maximum(F.S)
    if σmax > 1 + max(atol, rtol)
        throw(ArgumentError(lazy"The largest singular value $(σmax) of `S` exceeds one by more than the tolerances, so `S` is not passive."))
    end
    # sqrt(1 - σ^2), zero for a singular value within the tolerances above
    # one, and one for the directions beyond the smaller dimension of S
    c = sqrt.(max.(1 .- F.S .^ 2, 0))
    cn = [c; ones(eltype(c), n - k)]
    cm = [c; ones(eltype(c), m - k)]
    Σ = zeros(eltype(F.S), n, m)
    for i in 1:k
        Σ[i, i] = F.S[i]
    end

    Z = zeros(eltype(F.U), n, m)
    U = [F.U Z; Z' F.V] * [Σ Diagonal(cn); Diagonal(cm) -Σ'] * [F.Vt Z'; Z F.U']
    return U
end

"""
    Ymin_from_X(Omega, X; method = 1, hbar = 1)

Return the noise `Y = (hbar/2)*|im*Δ|`, with `Δ = Omega - X*Omega*transpose(X)`
and `|H|` the absolute value of a Hermitian matrix, which makes the Gaussian
map of the quadrature transformation `X` completely positive and trace
preserving for the symplectic form `Omega`, in the units of `hbar` of
[`is_cptp`](@ref). It is minimal: no noise below it in the positive
semi-definite order makes the map CPTP.

`method` selects the factorization `|im*Δ|` is computed from: `1`, the
default, an eigendecomposition of `im*Δ`; `2`, a Schur decomposition of
`Δ`; `3`, a singular value decomposition of `Δ`.
"""
function Ymin_from_X(Omega, X; method=1, hbar::Real=1)
    checkhbar(hbar)

    Delta = Matrix(Omega .- X * Omega * transpose(X))

    if method == 1
        # method 1, the default: an eigendecomposition, which is stable and
        # the fastest
        F = eigen(Hermitian(im * Delta))
        Ymin = real(F.vectors * Diagonal(abs.(F.values)) * F.vectors')

    elseif method == 2
        # method 2: schur
        F = schur(Delta)
        Ymin = F.Z * Diagonal(abs.(imag(F.values))) * transpose(F.Z)

    elseif method == 3
        # method 3: a singular value decomposition
        F = svd(Delta)
        Ymin = F.V * Diagonal(F.S) * F.Vt

    else
        error(lazy"Unknown method")
    end

    return (hbar / 2) * Ymin
end

"""
    Ymin_from_X_quadrature_pair(X; method = 1, hbar = 1)

[`Ymin_from_X`](@ref) for the quadrature transformation `X` in pair
operator order, in the units of `hbar` of [`is_cptp`](@ref).
"""
function Ymin_from_X_quadrature_pair(X; method=1, hbar::Real=1)
    Omega = symplectic_form_pair(size(X, 1) ÷ 2)
    return Ymin_from_X(Omega, X; method=method, hbar=hbar)
end

"""
    Ymin_from_X_quadrature_block(X; method = 1, hbar = 1)

[`Ymin_from_X`](@ref) for the quadrature transformation `X` in block
operator order, in the units of `hbar` of [`is_cptp`](@ref).
"""
function Ymin_from_X_quadrature_block(X; method=1, hbar::Real=1)
    Omega = symplectic_form_block(size(X, 1) ÷ 2)
    return Ymin_from_X(Omega, X; method=method, hbar=hbar)
end


"""
    A_B_to_symplectic_pair(A::AbstractMatrix, B::AbstractMatrix; atol = 0,
        rtol = ...)

Given `A` (2n×2n) and `B` (2n×4n) such that `A*Ω*A' + B*ΩE*B' = Ω`, return the
symplectic matrix `S = [A B; C D]` with respect to `Ωtot = Ω ⊕ ΩE`, in pair
operator order.

`C` (4n×2n) and `D` (4n×4n) are constructed here; `atol` and `rtol` are the
tolerances of the rank decisions that construction makes.

"""
function A_B_to_symplectic_pair(A::AbstractMatrix, B::AbstractMatrix;
    atol::Real=0, rtol::Real=defaultrtol(A, atol))

    # the number of system modes
    n = size(A, 1) ÷ 2

    # twice as many environment modes as system modes; the symplectic forms
    # are sparse, and so is their direct sum
    Ω = symplectic_form_pair(n)
    ΩE = symplectic_form_pair(2n)

    Ωtot = blockdiag(Ω, ΩE)

    W = hcat(A, B)                      # 2n × 6n

    # the complement rows N with W*Ωtot*N' = 0
    K = nullspace(W * Ωtot; atol=atol, rtol=rtol)  # 6n × k

    if size(K, 2) != 4n
        throw(ArgumentError(lazy"The rows of `[A B]` have rank $(6n - size(K,2)) rather than 2n = $(2n) to the tolerances `atol` and `rtol`, so they are not the first 2n rows of a symplectic matrix."))
    end
    N = K'
    # k by 6n, whose rows span the complement

    # the skew form induced on the complement
    G = N * Ωtot * N'

    if !isapprox(G, -transpose(G))
        error(lazy"G must be skew-symmetric.")
    end

    # the real Schur form G = Q*T*Q', quasi triangular for real G
    F = schur(G)
    Q = F.Z
    T = F.T

    _, d = schurblockscaling(T, 2n; inverse = true)

    Sscale = Diagonal(d)

    # the completed bottom block-row
    Wperp = Sscale * (Q') * N           # 4n × 6n

    C = Wperp[:, 1:2n]
    D = Wperp[:, 2n+1:end]

    return [A B; C D]
end

"""
    wmatrix(ws::AbstractVector, wp::Tuple, modes::AbstractVector{<:Tuple})

Return the `Nmodes` by `Nfreqs` matrix of frequencies for the signal, idlers,
and sidebands given the signal frequencies `ws`, pump frequencies `wp`, and
modes `modes`, each a tuple of as many integers as there are pumps: the
frequency of mode `modes[i]` at signal frequency `ws[j]` is
`ws[j] + sum(modes[i] .* wp)`. The element type is that of `ws` and `wp`
promoted.

# Examples
```jldoctest
julia> JosephsonCircuits.wmatrix(0.1:0.1:1.0,(1.0,),[(1,),(-1,)])
2×10 Matrix{Float64}:
  1.1   1.2   1.3   1.4   1.5   1.6   1.7   1.8   1.9  2.0
 -0.9  -0.8  -0.7  -0.6  -0.5  -0.4  -0.3  -0.2  -0.1  0.0
```
"""
function wmatrix(ws::AbstractVector, wp::NTuple{N,Number},
    modes::AbstractVector{NTuple{N,Int}}) where {N}
    # the output is Nmodes by Nfreqs
    w = zeros(promote_type(eltype(ws), map(typeof, wp)...), length(modes),
        length(ws))

    wmatrix!(w, ws, wp, modes)

    return w
end

"""
    wmatrix!(w::AbstractMatrix, ws::AbstractVector, wp::Tuple,
        modes::AbstractVector{<:Tuple})

In place version of [`wmatrix`](@ref), writing into `w`, of size
`(length(modes), length(ws))`, and returning it.
"""
function wmatrix!(w::AbstractMatrix, ws::AbstractVector, wp::NTuple{N,Number},
    modes::AbstractVector{NTuple{N,Int}}) where {N}
    if size(w) != (length(modes), length(ws))
        throw(DimensionMismatch(lazy"The size $(size(w)) of `w` must be the number of modes by the number of frequencies, $((length(modes), length(ws)))."))
    end
    for (j, wsj) in enumerate(ws)
        for (i, mode) in enumerate(modes)
            w[i, j] = wsj + dot(wp, mode)
        end
    end
    return w
end

"""
    interpolate_scattering(w0::AbstractVector, S::AbstractArray, w::AbstractArray;
        extrap = false, extrap_value = 0.0)

Interpolate the scattering parameters `S`, an array of size
`(nports, nports, length(w0))` tabulated at the strictly increasing
frequencies `w0`, at least three of them, onto the frequencies `w`. `w` may
be a matrix such as the one returned by [`wmatrix`](@ref), in which case
the result has one matrix per entry of it.

Each entry is interpolated as the product of a delay and a slowly varying
complex function: the delay `tau` is the least squares slope of the
unwrapped phase of the samples, and the real and imaginary parts of the
undelayed samples `S*cis(w0*tau)` are interpolated quadratically, then
multiplied by `cis(-w*tau)`. A delay line is then followed exactly however
often its phase winds between the samples, and a transmission zero, where
the magnitude has a kink and the phase a step, is crossed smoothly.

A negative frequency takes the complex conjugate of the value at its
magnitude. A frequency within rounding of an end of `w0`, such as an idler
of `wmatrix` formed from the signal and the pump, takes the value at that
end. With `extrap = true` a frequency whose magnitude is outside `w0`
takes the value `extrap_value` (conjugated at a negative frequency);
otherwise it is a `DomainError`.

# Examples
```jldoctest
w = 0.01:0.01:1.0
S = JosephsonCircuits.ABCD_tline(50,w)
isapprox(S,JosephsonCircuits.interpolate_scattering(w,S,w))

# output
true
```

"""
function interpolate_scattering(w0::AbstractVector, S::AbstractArray,
    w::AbstractArray; extrap::Bool=false, extrap_value=0.0)

    # check that S has only 3 dimensions. the first two are ports and the last
    # is frequencies.
    if ndims(S) != 3
        error("`S` must have 3 dimensions. The first two are ports and the third is frequencies.")
    end

    if length(w0) != size(S,3)
        error("The length of the third dimension of `S` must be equal to the number of frequencies.")
    end

    # the quadratic interpolants need three samples, in order; a repeated
    # frequency would put a zero interval under them
    if length(w0) < 3
        throw(ArgumentError(lazy"At least three frequencies are needed for the quadratic interpolation; got $(length(w0))."))
    end
    if !issorted(w0; lt = <=)
        throw(ArgumentError("The frequencies `w0` must be strictly increasing."))
    end

    # the first two dimensions are those of the scattering matrix; the rest
    # are those of `w`
    sizeout = NTuple{ndims(w) + 2,Int}(i == 1 || i == 2 ? size(S, i) : size(w, i - 2) for i in 1:(ndims(w)+2))

    # complex, whatever the type of the samples, since the interpolant of a
    # real entry carries a delay
    Sout = zeros(complex(float(eltype(S))), sizeout)

    # the band the samples cover, and the rounding admitted at its edges
    # (see edgetolerance); outside it the interpolants are an error unless
    # `extrap` gives the value there
    wmin, wmax = first(w0), last(w0)
    edgetol = edgetolerance(w0)

    # the undelayed samples of an entry
    r = Vector{eltype(Sout)}(undef, length(w0))

    # interpolate each entry over frequency, conjugating at negative
    # frequencies
    for i in 1: size(S, 1)
        for j in 1:size(S, 2)

            s = view(S, i, j, :)
            tau = fitdelay(w0, s)
            r .= s .* cis.(w0 .* tau)
            real_interp = FastInterpolations.quadratic_interp(w0, real.(r))
            imag_interp = FastInterpolations.quadratic_interp(w0, imag.(r))

            for c in CartesianIndices(axes(w))
                wi = abs(w[c])
                Sinterp = if wmin - edgetol <= wi <= wmax + edgetol
                    wi = clamp(wi, wmin, wmax)
                    complex(real_interp(wi), imag_interp(wi)) * cis(-wi * tau)
                elseif extrap
                    extrap_value
                else
                    throw(DomainError(w[c], lazy"The magnitude of the frequency is outside the band [$(wmin), $(wmax)] of the samples; pass `extrap = true` for the value `extrap_value` there."))
                end
                # conjugate at a negative frequency
                if w[c] < 0
                    Sout[i, j, c] = conj(Sinterp)
                else
                    Sout[i, j, c] = Sinterp
                end
            end
        end
    end

    return Sout
end

# the delay of the samples `s` at the frequencies `w`, the least squares
# slope of their unwrapped phase, negated
function fitdelay(w, s)
    phase = unwrap(angle.(s))
    wmean = sum(w) / length(w)
    phasemean = sum(phase) / length(phase)
    num = zero(wmean * phasemean)
    den = zero(wmean * wmean)
    for k in eachindex(w, phase)
        num += (w[k] - wmean) * (phase[k] - phasemean)
        den += (w[k] - wmean)^2
    end
    return -num / den
end
