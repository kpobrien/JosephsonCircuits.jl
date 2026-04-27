
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
    # return Int[0*I(n) I(n);-I(n) 0*I(n)]
end

"""
    direct_sum(A; n::Integer=1)

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
julia> JosephsonCircuits.indefinite_hermitian_form_block(2)
4×4 Diagonal{Int64, Vector{Int64}}:
 1  ⋅   ⋅   ⋅
 ⋅  1   ⋅   ⋅
 ⋅  ⋅  -1   ⋅
 ⋅  ⋅   ⋅  -1
"""
function indefinite_hermitian_form_block(n::Integer)
    d = Vector{Int}(undef, 2 * n)
    for i in 1:n
        d[i] = 1
    end
    for i in n+1:2*n
        d[i] = -1
    end
    # return Int[I(n) 0*I(n);0*I(n) -I(n)]
    return Diagonal(d)
end


function is_positive_semi_definite(M)
   # turn off error checking on cholesky, perform a
   # pivot so it works in the rank deficient case
   # and just check if it works to within numerical
   # error.
   C = cholesky(M,RowMaximum();check=false)
   if isempty(C.p)
      U = C.U
      L = C.L
      A = M
   else
      U = C.U[1:C.rank, :]
      L = C.L[:,1:C.rank]
      A = M[C.p, C.p]
   end
   return isapprox(A,U'*U) && isapprox(A,L*L')
end

function is_positive_definite(M)
    return isposdef(M)
end


"""
    is_unitary(M)

Return `true` if the matrix `M` is unitary, `M ∈ U(n)`, and `false` otherwise.

Tests if `M` satisfies the condition `M*M'==I` where `I` is the identity
matrix.

"""
function is_unitary(M)
    return isapprox(M*adjoint(M),I(size(M,1)))
end

"""
    is_orthogonal(M)

Return `true` if the matrix `M` is orthogonal, `M ∈ O(n)`, and `false`
otherwise.

Tests if `M` satisfies the condition `M*transpose(M)==I` where `I` is the
identity matrix.

"""
function is_orthogonal(M)
    return isapprox(M*transpose(M),I(size(M,1)))
end

"""
    is_symplectic(Ω, S)

Return `true` if the matrix `S` is symplectic, `S ∈ Sp(2n, ℝ)` or
`S ∈ Sp(2n, ℂ)`, and `false` otherwise.

Tests if `S` satisfies the symplectic condition `S*Ω*transpose(S)==Ω` where
`Ω` is a user supplied symplectic form.

See also [`symplectic_form_block`](@ref),
[`symplectic_form_pair`](@ref), [`is_symplectic_block`](@ref), and
[`is_symplectic_pair`](@ref).

"""
function is_symplectic(Ω, S)
    # perhaps i should check if Ω is a symplectic form
    # that would guard against swapping the arguments
    return isapprox(S*Ω*transpose(S), Ω)
end

"""
    is_symplectic_block(S)

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
    is_symplectic_pair(S)

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
    isorthogonal_symplectic_block(M)

Return `true` if the matrix `M` is orthogonal symplectic,
`M ∈ Sp(2n, ℝ) ∩ O(2n) ≅ U(n)`, with block operator order and `false`
otherwise.

"""
function is_orthogonal_symplectic_block(M)
    return is_symplectic_block(M) && is_orthogonal(M)
end

"""
    is_orthogonal_symplectic_pair(M)

Return `true` if the matrix `M` is orthogonal symplectic,
`M ∈ Sp(2n, ℝ) ∩ O(2n) ≅ U(n)`, with pair operator order and `false`
otherwise.

"""
function is_orthogonal_symplectic_pair(M)
    return is_symplectic_pair(M) && is_orthogonal(M)
end

"""
    is_conjugate_symplectic_pair(M)

Return `true` if the matrix `M` is conjugate symplectic,
, with pair operator order and `false`
otherwise.

"""
function is_conjugate_symplectic_pair(M)
    Omega = symplectic_form_pair(size(M, 1) ÷ 2)
    return is_pseudo_unitary(Omega, M)
end

"""
    is_conjugate_symplectic_block(M)

Return `true` if the matrix `M` is conjugate symplectic,
, with block operator order and `false`
otherwise.

"""
function is_conjugate_symplectic_block(M)
    Omega = symplectic_form_block(size(M, 1) ÷ 2)
    return is_pseudo_unitary(Omega, M)
end


"""
    is_pseudo_unitary(Σ,M)

Return `true` if the matrix `M` is pseudo-unitary, `M ∈ U(n, n)`, and `false`
otherwise.

Tests if `M` satisfies the pseudo-unitary condition `M*Σ*M'==Σ` where `Σ` is a
user supplied matrix representing the indefinite Hermitian form.

See also [`indefinite_hermitian_form_pair`](@ref),
[`indefinite_hermitian_form_block`](@ref), [`is_pseudo_unitary_block`](@ref),
and [`is_pseudo_unitary_pair`](@ref).

"""
function is_pseudo_unitary(Sigma, M)
    return isapprox(M*Sigma*M', Sigma)
end

"""
    is_pseudo_unitary_block(M)

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
    is_pseudo_unitary_pair(M)

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
    is_positive_definite_symplectic_block(S)

Return `true` if the matrix `S` is positive definite and symplectic,
`S ∈ Sp(2n, ℝ)` or `S ∈ Sp(2n, ℂ)`, with block operator order and `false`
otherwise.

"""
function is_positive_definite_symplectic_block(S)
    # the tests for ishermitian and issymmetric are exact, but we want
    # to test for equality up to numerical error so use isappox
    return is_symplectic_block(S) && isapprox(S,S') && is_positive_definite(Hermitian(S))
end

"""
    is_positive_definite_symplectic_pair(S)

Return `true` if the matrix `S` is positive definite and symplectic,
`S ∈ Sp(2n, ℝ)` or `S ∈ Sp(2n, ℂ)`, with pair operator order and
`false` otherwise.

"""
function is_positive_definite_symplectic_pair(S)
    # the tests for ishermitian and issymmetric are exact, but we want
    # to test for equality up to numerical error so use isappox
    return is_symplectic_pair(S) && isapprox(S,S') && is_positive_definite(Hermitian(S))
end

"""
    is_bogoliubov_pair(M)

Return `true` if the matrix `M` is Bogoliubov, `M ∈ Sp(2n, ℂ) ∩ U(n, n)`, with
pair operator order and `false` otherwise.

"""
function is_bogoliubov_pair(M)
    return is_pseudo_unitary_pair(M) && is_symplectic_pair(M)
end

"""
    is_bogoliubov_block(M)

Return `true` if the matrix `M` is Bogoliubov, `M ∈ Sp(2n, ℂ) ∩ U(n, n)`, with
block operator order and `false` otherwise.

"""
function is_bogoliubov_block(M)
    return is_pseudo_unitary_block(M) && is_symplectic_block(M)
end

"""
    is_orthogonal_bogoliubov_block(M)

Return `true` if the matrix `M` is orthogonal Bogoliubov,
`M ∈ Sp(2n, ℂ) ∩ U(n, n) ∩ U(2n) ≅ U(n)`, with pair operator order and
`false` otherwise.

"""
function is_orthogonal_bogoliubov_block(M)
    return is_bogoliubov_block(M) && is_unitary(M)
end

"""
    is_orthogonal_bogoliubov_pair(M)

Return `true` if the matrix `M` is orthogonal Bogoliubov,
`M ∈ Sp(2n, ℂ) ∩ U(n, n) ∩ U(2n) ≅ U(n)`, with pair operator order and
`false` otherwise.

"""
function is_orthogonal_bogoliubov_pair(M)
    return is_bogoliubov_pair(M) && is_unitary(M)
end


# should i add is_symplectic_bogoliubov_pair and
# is_symplectic_bogoliubov_block ?
# should i include the bogoliubov form where there are two seperate matrices
# one for the linear processes and the second for parametric processes
# how about a random CPTP map?
# rand_cptp_block(X,Y) and rand_cptp_pair(X,Y) - should these
# satisfy the minimum uncertainty relation?
# does Ymin need to change depending on pair and block? yes
# all of these depend on that for computing delta.

# and is_cptp_block(n) - would have to check if satisfies uncertainty relations
# is_cptp_pair(n)

"""

Eq. 5.37 from Serafini

"""
function is_cptp(Omega, X, Y)
    delta = Matrix(Omega) - X * Omega * X'
    # delta = Matrix(Omega) - X * Omega * transpose(X)
    K = Y + im * delta
    # check if the matrix is approximately Hermitian
    if !isapprox(K, Hermitian(K))
        return false
    else
        return is_positive_semi_definite(Hermitian(K))
    end
end

function is_cptp_block(X, Y)
    n = size(X, 1) ÷ 2
    Omega = symplectic_form_block(n)
    return is_cptp(Omega, X, Y)
end

function is_cptp_pair(X, Y)
    n = size(X, 1) ÷ 2
    Omega = symplectic_form_pair(n)
    return is_cptp(Omega, X, Y)
end

function is_cptp_bogoliubov_pair(X, Y)
    n = size(X, 1) ÷ 2
    # note we need to multiply this by im to get
    # the Robertson-Schrödinger inequality
    Omega = im*indefinite_hermitian_form_pair(n)
    return is_cptp(Omega, X, Y)
end

function is_cptp_bogoliubov_block(X, Y)
    n = size(X, 1) ÷ 2
    # note we need to multiply this by im to get
    # the Robertson-Schrödinger inequality
    Omega = im*indefinite_hermitian_form_block(n)
    return is_cptp(Omega, X, Y)
end

function rand_positive_definite(T, n::Integer)
    A = rand(T,n,n)
    return A*A'
end

function rand_positive_definite(n::Integer)
    return rand_positive_definite(Float64, n)
end

function rand_unitary(T, n::Integer)
    A = rand(T,n,n)
    # make a skew-Hermitian matrix
    H = (A - A')/2
    return exp(H)
end

function rand_unitary(n::Integer)
    return rand_unitary(Complex{Float64},n)
end

function rand_orthogonal(n::Integer)
    return rand_unitary(Float64,n)
end


"""
   rand_positive_semi_definite(T,n,m)


"""
function rand_positive_semi_definite(T,n,m)
   A = rand(T,n+m,n)
   return A*A'
end

function rand_positive_semi_definite(n,m)
   return rand_positive_semi_definite(Float64,n,m)
end


"""
    rand_symplectic_block(T, n::Integer)

Return a random `2n x 2n` symplectic matrix `S`, `S ∈ Sp(2n, ℝ)` or
`S ∈ Sp(2n, ℂ)`, depending on the type `T` with block operator order.

"""
function rand_symplectic_block(T, n::Integer)
    A = randn(T, 2*n, 2*n)
    Omega = symplectic_form_block(n)
    # generate a random symmetric matrix
    M = (transpose(A) + A) / 2
    # matrix exponential seems result in fewer numerical errors than
    # Cayley transform
    return exp(Omega * M)
    # return cayley_transform(Omega,M)
end

"""
    rand_symplectic_block(n::Integer)

Return a random `2n x 2n` symplectic matrix `S`, `S ∈ Sp(2n, ℝ)`, with block
operator order.

"""
function rand_symplectic_block(n::Integer)
    return rand_symplectic_block(Float64, n)
end

function rand_symplectic_pair(T, n::Integer)
    A = randn(T, 2*n, 2*n)
    Omega = symplectic_form_pair(n)
    # generate a random symmetric matrix
    M = (transpose(A) + A) / 2
    # matrix exponential seems result in fewer numerical errors than
    # Cayley transform
    return exp(Omega * M)
    # return cayley_transform(Omega,M)
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
    P = (P-transpose(P))/2

    Q = randn(T, n, n)
    # make Q complex symmetric
    Q = (Q+transpose(Q))/2

    # assemble M
    M = [-Q -P;P -Q]

    # compute the symplectic matrix using the Cayley transform
    Omega = symplectic_form_block(n)
    return cayley_transform(Omega,M)
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

    # M = randn(T,n,n)
    # N = randn(T,n,n)

    # # make a positive definite matrix (symmetric if real or Hermitian if
    # # complex)
    # U = M*M'

    # # make a symmetric matrix
    # V = (N+transpose(N))/2

    # # compute the inverse using a Cholesky decomposition
    # F = cholesky(Hermitian(U))  
    # UinvV = F\V
    # Uinv = F\I

    # # construct a positive definite symplectic matrix
    # return [Uinv -UinvV;-V'*Uinv transpose(U)+V'*UinvV]

    # Uinv = inv(U)
    # return [Uinv -Uinv*V;-V'*Uinv transpose(U)+V'*Uinv*V]

    # the above constructions work and are nice because they are explicit
    # which may be useful for derivations, but they seem to be numerically
    # very brittle. taking the positive definite symplectic term from the
    # polar decomposition of a random symplectic matrix seems to be more
    # robust.
    A = rand_symplectic_block(T,n)
    S,_ = polar(A)
    return S

    # multiplying a symplectic matrix by its adjoint is even easier!
    # A = rand_symplectic_block(T,n)
    # return A*A'

end

"""
    rand_positive_definite_symplectic_block(n::Integer)

Return a random `2n x 2n` positive definite symplectic matrix `S`,
`S ∈ Sp(2n, ℝ)`, with block operator order.

"""
function rand_positive_definite_symplectic_block(n::Integer)
    return rand_positive_definite_symplectic_block(Float64,n)
end

"""
    rand_positive_definite_symplectic_pair(T, n::Integer)

Return a random `2n x 2n` positive definite symplectic matrix `S`,
`S ∈ Sp(2n, ℝ)` or `S ∈ Sp(2n, ℂ)`, depending on the type `T` with pair
operator order.

"""
function rand_positive_definite_symplectic_pair(T, n::Integer)
    A = rand_symplectic_pair(T,n)
    return A*A'
end

"""
    rand_positive_definite_symplectic_pair(n::Integer)

Return a random `2n x 2n` positive definite symplectic matrix `S`,
`S ∈ Sp(2n, ℝ)`, with pair operator order.

"""
function rand_positive_definite_symplectic_pair(n::Integer)
    return rand_positive_definite_symplectic_pair(Float64,n)
end


"""


"""
function rand_conjugate_symplectic_block(T, n::Integer)
    # X = randn(T, 2n, 2n)

    # Omega = symplectic_form_block(n)

    # # project onto sp†(2n)
    # A = 1/2*(X-Omega*X'*Omega')

    # # compute the symplectic matrix using the Cayley transform
    # return cayley_transform(I(2n),A)
    # # return exp(A)

    A = randn(T, 2*n, 2*n)
    Omega = symplectic_form_block(n)
    # generate a random Hermitian matrix
    M = (A + A') / 2
    # matrix exponential seems result in fewer numerical errors than
    # Cayley transform
    return exp(Omega * M)
    # return cayley_transform(Omega,M)

end

function rand_conjugate_symplectic_block(n::Integer)
    return rand_conjugate_symplectic_block(Complex{Float64},n)
end

function rand_conjugate_symplectic_pair(T, n::Integer)
    return block_to_pair(rand_conjugate_symplectic_block(T, n))
end

function rand_conjugate_symplectic_pair(n::Integer)
    return rand_conjugate_symplectic_pair(Complex{Float64},n)
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
    P = (P+transpose(P))/2

    Q = randn(T, n, n)
    # make Q skew-Hermitian
    Q = (Q-Q')/2

    # assemble M
    M = [P Q;transpose(Q) -conj(P)]

    # compute the symplectic matrix using the Cayley transform
    Omega = symplectic_form_block(n)
    # matrix exponential seems result in fewer numerical errors than
    # Cayley transform
    return exp(Omega * M)
    # return cayley_transform(Omega,M)

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
    Q = (Q-Q')/2

    # assemble M
    M = [0*I(n) Q;transpose(Q) 0*I(n)]

    # compute the symplectic matrix using the Cayley transform
    Omega = symplectic_form_block(n)
    return cayley_transform(Omega,M)
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
    A = randn(T, 2*n, 2*n)
    K = indefinite_hermitian_form_block(n)
    # generate a random skew-Hermitian
    M = (A - A')/2
    # return exp(K * M)
    return cayley_transform(K,M)
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
    A = randn(T, 2*n, 2*n)
    K = indefinite_hermitian_form_pair(n)
    # generate a random skew-Hermitian
    M = (A - A')/2
    # return exp(K * M)
    return cayley_transform(K,M)
end

function cayley_transform(Omega,M)
    # add a check to verify the sizes of Omega and M are the same
    n = size(M,1)
    S = (I(n) + Omega*M) * inv(I(n) - Omega*M)
    # S = qr(I(n) - Omega*M)\(I(n) + Omega*M)
    return S
end

function rand_cptp_block(T, nsys::Integer; nenv::Integer=nsys,
    sigma_env=2 * I(2 * nsys))
    # need to start with an pair matrix
    S = rand_symplectic_pair(T, nsys + nenv)
    # now convert each of these blocks to the block form
    A = pair_to_block(S[1:2*nsys, 1:2*nsys])
    B = pair_to_block(S[1:2*nsys, 2*nsys+1:end])

    X = A
    # Y = B * sigma_env * B'
    Y = B*sigma_env*transpose(B)
    return (X=X, Y=Y)
end

function rand_cptp_block(nsys::Integer; nenv::Integer=nsys)
    return rand_cptp_block(Float64, nsys; nenv=nsys)
end

function rand_cptp_pair(T, nsys::Integer; nenv::Integer=nsys,
    sigma_env = 2*I(2*nsys))
    S = rand_symplectic_pair(T, nsys + nenv)
    A = S[1:2*nsys, 1:2*nsys]
    B = S[1:2*nsys, 2*nsys+1:end]

    X = A
    # Y = B * sigma_env * B'
    Y = B*sigma_env*transpose(B)
    return (X=X, Y=Y)
end

function rand_cptp_pair(nsys::Integer; nenv::Integer=nsys)
    return rand_cptp_pair(Float64, nsys; nenv=nsys)
end


function rand_cptp_bogoliubov_pair(T, nsys::Integer; nenv::Integer=nsys,
    sigma_env = 2*I(2*nsys))
    S = rand_bogoliubov_pair(T, nsys + nenv)
    A = S[1:2*nsys, 1:2*nsys]
    B = S[1:2*nsys, 2*nsys+1:end]

    X = A
    Y = B*sigma_env*B'
    return (X=X, Y=Y)
end

function rand_cptp_bogoliubov_pair(nsys::Integer; nenv::Integer=nsys)
    return rand_cptp_bogoliubov_pair(Complex{Float64}, nsys; nenv=nsys)
end

function rand_cptp_bogoliubov_block(T, nsys::Integer; nenv::Integer=nsys,
    sigma_env = 2*I(2*nsys))
    # need to start with an pair matrix
    S = rand_bogoliubov_pair(T, nsys + nenv)
    # now convert each of these blocks to the block form
    A = pair_to_block(S[1:2*nsys, 1:2*nsys])
    B = pair_to_block(S[1:2*nsys, 2*nsys+1:end])

    X = A
    Y = B*sigma_env*B'
    return (X=X, Y=Y)
end

function rand_cptp_bogoliubov_block(nsys::Integer; nenv::Integer=nsys)
    return rand_cptp_bogoliubov_block(Complex{Float64}, nsys; nenv=nsys)
end


# function rand_cptp_block(T,n;method = 1, shift = 100*n*eps(T))
#     X = randn(T,2*n, 2*n)
#     Ymin = Ymin_from_X_block(X;method=method)
#     Y = Ymin+shift*I(2*n)
#     return X,(Y+Y')/2
# end

# function rand_cptp_block(n;method = 1, shift = 100*n*eps(Float64))
#     return rand_cptp_block(Float64,n;method = method, shift = shift)
# end

# function rand_cptp_pair(T,n;method = 1, shift = 100*n*eps(T))
#     X = randn(T,2*n, 2*n)
#     Ymin = Ymin_from_X_pair(X;method=method)
#     Y = Ymin+shift*I(2*n)
#     return X,(Y+Y')/2
# end

# function rand_cptp_pair(n;method = 1, shift = 100*n*eps(Float64))
#     return rand_cptp_pair(Float64,n;method = method, shift = shift)
# end

# function rand_cptp_pair(n;method = 1, shift = 100*eps())
#     X = randn(Float64,2*n, 2*n)
#     Ymin = Ymin_from_X_pair(X;method=method)
#     Y = Ymin+shift*I(2*n)
#     return X,(Y+Y')/2
# end


"""
    block_to_pair_perm(n::Int)

Return a vector `p` which permutes the block operator ordering into the
pair operator ordering.

Return a `2n` length vector `p` which permutes the block operator ordering
`r = [x1,...,xn,p1,...,pn]` to the pair operator ordering
`r[p] = (x_1,p_1,...,x_n,p_n)`.

# Examples
```jldoctest
@variables x1 x2 x3 x4 p1 p2 p3 p4
r = [x1, x2, x3, x4, p1, p2, p3, p4]
p = JosephsonCircuits.block_to_pair_perm(4)
r[p]

# output
8-element Vector{Num}:
 x1
 p1
 x2
 p2
 x3
 p3
 x4
 p4
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
@variables x1 x2 x3 x4 p1 p2 p3 p4
r = [x1, p1, x2, p2, x3, p3, x4, p4]
p = JosephsonCircuits.pair_to_block_perm(4)
r[p]

# output
8-element Vector{Num}:
 x1
 x2
 x3
 x4
 p1
 p2
 p3
 p4
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
    #colptr = Vector{Int}(undef,2*n+1)
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
@variables x1 x2 x3 x4 p1 p2 p3 p4
r = [x1, x2, x3, x4, p1, p2, p3, p4]
JosephsonCircuits.block_to_pair(r)

# output
8-element Vector{Num}:
 x1
 p1
 x2
 p2
 x3
 p3
 x4
 p4
 ```
"""
function block_to_pair(r::AbstractVector)
    p = block_to_pair_perm(length(r) ÷ 2)
    return r[p]
end

"""
    block_to_pair2(r::AbstractVector)

# Examples
```jldoctest
@variables x1 x2 x3 x4 p1 p2 p3 p4
r = [x1, x2, x3, x4, p1, p2, p3, p4]
JosephsonCircuits.block_to_pair2(r)

# output
8-element Vector{Num}:
 x1
 p1
 x2
 p2
 x3
 p3
 x4
 p4
```
"""
function block_to_pair2(r::AbstractVector)
    R = R_block_to_pair(length(r) ÷ 2)
    return R * r
end

"""
    block_to_pair(S::AbstractMatrix)

"""
function block_to_pair(S::AbstractMatrix)
    p1 = block_to_pair_perm(size(S, 1) ÷ 2)
    p2 = block_to_pair_perm(size(S, 2) ÷ 2)
    return S[p1, p2]
end

"""
    block_to_pair2(S::AbstractMatrix)

"""
function block_to_pair2(S::AbstractMatrix)
    R = R_block_to_pair(size(S, 2) ÷ 2)
    return R * S * R'
end

"""
    R_pair_to_block(n::Integer)

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
    #colptr = Vector{Int}(undef,2*n+1)
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
@variables x1 x2 x3 x4 p1 p2 p3 p4
r = [x1, p1, x2, p2, x3, p3, x4, p4]
JosephsonCircuits.pair_to_block(r)

# output
8-element Vector{Num}:
 x1
 x2
 x3
 x4
 p1
 p2
 p3
 p4
```
"""
function pair_to_block(r::AbstractVector)
    p = pair_to_block_perm(length(r) ÷ 2)
    return r[p]
end

"""
    pair_to_block2(r::AbstractVector)

# Examples
```jldoctest
@variables x1 x2 x3 x4 p1 p2 p3 p4
r = [x1, p1, x2, p2, x3, p3, x4, p4]
JosephsonCircuits.pair_to_block2(r)

# output
8-element Vector{Num}:
 x1
 x2
 x3
 x4
 p1
 p2
 p3
 p4
```
"""
function pair_to_block2(r::AbstractVector)
    R = R_pair_to_block(length(r) ÷ 2)
    return R * r
end

"""
    pair_to_block(S::AbstractMatrix)

"""
function pair_to_block(S::AbstractMatrix)
    p1 = pair_to_block_perm(size(S, 1) ÷ 2)
    p2 = pair_to_block_perm(size(S, 2) ÷ 2)
    return S[p1, p2]
end

"""
    pair_to_block2(S::AbstractMatrix)

"""
function pair_to_block2(S::AbstractMatrix)
    R = R_pair_to_block(size(S, 2) ÷ 2)
    return R * S * R'
end

"""
    R_ladder_to_quadrature_pair(n::Integer)

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

# Examples
```jldoctest
@variables a1 a2 a3 a4 adag1 adag2 adag3 adag4
r = [a1, adag1, a2, adag2, a3, adag3, a4, adag4]
JosephsonCircuits.ladder_to_quadrature_pair(r)

# output
8-element Vector{Complex{Num}}:
   0.7071067811865475a1 + 0.7071067811865475adag1
 (-0.7071067811865475a1 + 0.7071067811865475adag1)*im
   0.7071067811865475a2 + 0.7071067811865475adag2
 (-0.7071067811865475a2 + 0.7071067811865475adag2)*im
   0.7071067811865475a3 + 0.7071067811865475adag3
 (-0.7071067811865475a3 + 0.7071067811865475adag3)*im
   0.7071067811865475a4 + 0.7071067811865475adag4
 (-0.7071067811865475a4 + 0.7071067811865475adag4)*im
```
"""
function ladder_to_quadrature_pair(r::AbstractVector)
    R = R_ladder_to_quadrature_pair(length(r) ÷ 2)
    return R * r
end

"""
    ladder_to_quadrature_pair(S::AbstractMatrix)

"""
function ladder_to_quadrature_pair(S::AbstractMatrix)
    R = R_ladder_to_quadrature_pair(size(S, 2) ÷ 2)
    return R * S * R'
end

"""
    R_quadrature_to_ladder_pair(n::Integer)

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

# Examples
```jldoctest
@variables x1 x2 x3 x4 p1 p2 p3 p4
r = [x1, p1, x2, p2, x3, p3, x4, p4]
JosephsonCircuits.quadrature_to_ladder_pair(r)

# output
8-element Vector{Complex{Num}}:
 0.7071067811865475x1 + 0.7071067811865475im*p1
 0.7071067811865475x1 - 0.7071067811865475im*p1
 0.7071067811865475x2 + 0.7071067811865475im*p2
 0.7071067811865475x2 - 0.7071067811865475im*p2
 0.7071067811865475x3 + 0.7071067811865475im*p3
 0.7071067811865475x3 - 0.7071067811865475im*p3
 0.7071067811865475x4 + 0.7071067811865475im*p4
 0.7071067811865475x4 - 0.7071067811865475im*p4
```
"""
function quadrature_to_ladder_pair(r::AbstractVector)
    R = R_quadrature_to_ladder_pair(length(r) ÷ 2)
    return R * r
end

"""
    quadrature_to_ladder_pair(S::AbstractMatrix)

"""
function quadrature_to_ladder_pair(S::AbstractMatrix)
    R = R_quadrature_to_ladder_pair(size(S, 2) ÷ 2)
    return R * S * R'
end

"""
    R_ladder_to_quadrature_block(n::Integer)

# Examples
```jldoctest
julia> JosephsonCircuits.R_ladder_to_quadrature_block(1)
2×2 Matrix{ComplexF64}:
 0.707107+0.0im       0.707107+0.0im
      0.0-0.707107im       0.0+0.707107im
```
"""
function R_ladder_to_quadrature_block(n::Integer)
    # I should convert this to a sparse matrix form
#    return Complex{Float64}[I(n) I(n); -im*I(n) im*I(n)] / sqrt(2)
    return [I(n)/sqrt(2) I(n)/sqrt(2); -im*I(n)/sqrt(2) im*I(n)/sqrt(2)]
end

"""
    ladder_to_quadrature_block(r::AbstractVector)

# Examples
```jldoctest
@variables a1 a2 a3 a4 adag1 adag2 adag3 adag4
r = [a1, a2, a3, a4, adag1, adag2, adag3, adag4]
JosephsonCircuits.ladder_to_quadrature_block(r)

# output
8-element Vector{Complex{Num}}:
   0.7071067811865475a1 + 0.7071067811865475adag1
   0.7071067811865475a2 + 0.7071067811865475adag2
   0.7071067811865475a3 + 0.7071067811865475adag3
   0.7071067811865475a4 + 0.7071067811865475adag4
 (-0.7071067811865475a1 + 0.7071067811865475adag1)*im
 (-0.7071067811865475a2 + 0.7071067811865475adag2)*im
 (-0.7071067811865475a3 + 0.7071067811865475adag3)*im
 (-0.7071067811865475a4 + 0.7071067811865475adag4)*im
```

"""
function ladder_to_quadrature_block(r::AbstractVector)
    R = R_ladder_to_quadrature_block(length(r) ÷ 2)
    return R * r
end

"""
    ladder_to_quadrature_block(S::AbstractMatrix)

"""
function ladder_to_quadrature_block(S::AbstractMatrix)
    R = R_ladder_to_quadrature_block(size(S, 2) ÷ 2)
    return R * S * R'
end

"""
    R_quadrature_to_ladder_block(n::Integer)

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
    # I should convert this to a sparse matrix form
    return [I(n)/sqrt(2) im*I(n)/sqrt(2); I(n)/sqrt(2) -im*I(n)/sqrt(2)]
end

"""
    quadrature_to_ladder_block(r::AbstractVector)

# Examples
```jldoctest
@variables x1 x2 x3 x4 p1 p2 p3 p4
r = [x1, x2, x3, x4, p1, p2, p3, p4]
JosephsonCircuits.quadrature_to_ladder_block(r)

# output
8-element Vector{Complex{Num}}:
 0.7071067811865475x1 + 0.7071067811865475im*p1
 0.7071067811865475x2 + 0.7071067811865475im*p2
 0.7071067811865475x3 + 0.7071067811865475im*p3
 0.7071067811865475x4 + 0.7071067811865475im*p4
 0.7071067811865475x1 - 0.7071067811865475im*p1
 0.7071067811865475x2 - 0.7071067811865475im*p2
 0.7071067811865475x3 - 0.7071067811865475im*p3
 0.7071067811865475x4 - 0.7071067811865475im*p4
```
"""
function quadrature_to_ladder_block(r::AbstractVector)
    R = R_quadrature_to_ladder_block(length(r) ÷ 2)
    return R * r
end

"""
    quadrature_to_ladder_block(S::AbstractMatrix)

"""
function quadrature_to_ladder_block(S::AbstractMatrix)
    R = R_quadrature_to_ladder_block(size(S, 2) ÷ 2)
    return R * S * R'
end


"""
   optimum_eigenvalue_angle(values)

I want to shift the midpoint of the arc to -1+0im or to ±π radians.
idea from here: https://github.com/XanaduAI/thewalrus/pull/403

# Examples
```
using Plots
t = range(0, 4π, length = 100)
values = randn(Complex{Float64},10)
values ./= abs.(values)
optimum_angle, optimum_rotation = optimum_eigenvalue_angle(values)
shift = exp(im*optimum_rotation)
plot(cos.(t), sin.(t))
plot!(real.(values),imag.(values);seriestype=:scatter)
plot!([cos(optimum_angle)],[sin(optimum_angle)];seriestype=:scatter)
plot!(real.(shift*values),imag.(shift*values);seriestype=:scatter)
```
"""
function optimum_eigenvalue_angle(values; target_angle=pi)
    # normalize the eigenvalues so we can view them as vectors on the unit
    # circle
    normalized_values = values ./ abs.(values)

    # sort them by angle
    sorted_values = sort(normalized_values, by=angle)

    # the optimum angle to rotate by is in the middle of the widest
    # arc. find the width of the arc using the dot product.
    max_arc_width = acos(real(conj(sorted_values[1]) * sorted_values[end]))
    # find the midpoint of the arc by averaging the two vectors (don't bother
    # normalizing). do both of these operations for the arc formed by the first
    # and last points.
    arc_midpoint = sorted_values[1] * exp(im * max_arc_width / 2)
    # println(max_arc_width," ",arc_midpoint)

    # perform these operations for the rest of the arcs
    # then return the midpoint of the widest arc
    for i in 1:length(sorted_values)-1
        arc_width = acos(real(conj(sorted_values[i]) * sorted_values[i+1]))
        # println(arc_width," ")
        if arc_width > max_arc_width
            max_arc_width = arc_width
            arc_midpoint = sorted_values[i] * exp(im * max_arc_width / 2)
            # println(max_arc_width," ",arc_midpoint)
        end
    end
    return angle(arc_midpoint), angle(exp(im * target_angle) / arc_midpoint)
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

    # matrix must be square

    # P = sqrt(Hermitian(A*A'))
    F = svd(A)
    P = F.U * Diagonal(F.S) * F.U'

    # W = inv(P)*A
    U = F.U * F.Vt

    # probably should define a struct
    # model off skewlinearalgebra.jl
    return (P=P, U=U)
end

"""
    williamson_pair(M::AbstractMatrix{<:Real})

For a symmetric positive semi-definite matrix `M`, return a vector of values `d`
and a real symplectic matrix `S` such that `M = S Diagonal(d) S^T`. `S` is
symplectic with respect to the pair ordered symplectic form `Ω`.

The values `d` are unique but the matrix `S` is not.

At some point evaluate whether the method in this reference
http://arxiv.org/abs/2108.05364v2 is better than the one we are using. I
switch to the Schur decomposition based method from [2] because anything based
on eigedecomposition may have problems with degenerate eigenvalues.

# References
[1] M. Idel, S. Soto Gaona, and M. M. Wolf, “Perturbation bounds for
Williamson’s symplectic normal form,” Linear Algebra and its Applications,
vol. 525, pp. 45–58, Jul. 2017, doi: 10.1016/j.laa.2017.03.013.
[2] M. Houde, W. McCutcheon, and N. Quesada, “Matrix decompositions in Quantum
Optics: Takagi/Autonne, Bloch-Messiah/Euler, Iwasawa, and Williamson,” Can.
J. Phys., vol. 102, no. 10, pp. 497–507, Oct. 2024, doi: 10.1139/cjp-2024-0070.
"""
function williamson_pair(M::AbstractMatrix{<:Real})
    n = size(M, 1) ÷ 2
    Omega = symplectic_form_pair(n)
    d, S = _williamson(Omega,M)
    return d, S
end

function williamson_block(M::AbstractMatrix{<:Real})
    n = size(M, 1) ÷ 2
    # alternatively, we could call williamson_pair(block_to_pair(M)) and 
    # then convert the outputs back to block ordering. I prefer to do it this
    # way because 
    Omega = symplectic_form_block(n)
    d, S = _williamson(Omega,M)
    return pair_to_block(d), S*R_block_to_pair(n)
end

function _williamson(Omega, M::AbstractMatrix{<:Real})
    n = size(M, 1) ÷ 2

    # compute the pivoted cholesky factorization and turn check to false. this
    # will fail if the matrix is not positive semi-definite. It's more robust
    # against small negative eigenvalues due to floating point errors than
    # non-pivoted cholesky.

    # evaluate this in a function because sparse and non-sparse cholesky have
    # different syntax
    L,rankL = choleskyLr(M)

    K = transpose(L)*Omega*L

    # check that C.rank is even
    # i've noticed that sometimes cholesky produces a rank that is larger
    # than the rank i set. 
    if isodd(rankL)
        # println(rankL)
        # rankL = rankL - 1
        error("The rank must be even.")
    end
    r = rankL ÷ 2

    # K is skew-symmetric so the real schur decomposition will return
    # 2x2 blocks of either all zeros or identical values on the off diagonals
    # differing by a sign.
    # use the schur decomposition of A to compute the symplectric normal form.
    F = schur(K)
    T = Matrix(F.T)
    Z = Matrix(F.Z)

    # loop through the blocks
    d = zeros(eltype(T),2*r)
    phi = zeros(eltype(T),2*r)

    # this naturally gives the pair ordering. for block need to permute the
    # order, which we do outside of this function.
    # for i in 1:2:2*n
    for i in 1:2:2*r
        a = (T[i, i+1] - T[i+1, i]) / 2
        d[i] = d[i+1] = abs(a)
        s = inv(sqrt(abs(a)))
        phi[i] = phi[i+1] = s
        phi[i+1] *= sign(a)
    end

    S1 = L*Z*Diagonal(phi)

    # if r == n then just return the values
    # otherwise, pad it out using the symplectic complement
    if r == n
        return (d = d, S = S1)
    else
        S2 = symplectic_complement(Omega,S1)
        return (d = vcat(d,zeros(eltype(d),2*(n-r))), S = hcat(S1,S2))
    end
end

function choleskyLr(M::AbstractArray)
    # compute the pivoted cholesky factorization and turn check to false. this
    # will fail if the matrix is not positive semi-definite. It's more robust
    # against small negative eigenvalues due to floating point errors than
    # non-pivoted cholesky.
    C = cholesky(M,RowMaximum();check=false)

    if isodd(C.rank)
        # println(rankL)
        rankL = C.rank - 1
        # error("The rank must be even.")
    else
        rankL = C.rank
    end

    L = C.L[invperm(C.p),1:rankL]
    return L,rankL
end

function choleskyLr(M::SparseMatrixCSC)
    C = cholesky(M)
    return Matrix(sparse(C.L)[invperm(C.p),:]), size(M,1)
end

function symplectic_complement(Omega,S1)

    n = size(S1,1) ÷ 2
    r = size(S1,2) ÷ 2

    # it looks like the rank deficient case is failing in the block format
    # make sense because i'm using the pair symplectic form for
    # symplectic_complement. i think i also need to 
    # Omega = symplectic_form_pair(n)

    # ker(S1' Ω) is the symplectic complement of range(S1) = range(M).
    # 2n × (2n-r), orthonormal cols
    S2 = nullspace(transpose(S1) * Omega)
    # (2n-r) × (2n-r), skew, full-rank
    K2 = transpose(S2) * Omega * S2
    F2 = schur(K2)
    T2 = Matrix(F2.T); Z2 = Matrix(F2.Z)
    phi2 = zeros(eltype(T2), 2n - 2r)
    for i in 1:2:(2n - 2r)
        a = (T2[i, i+1] - T2[i+1, i]) / 2
        s = inv(sqrt(abs(a)))
        phi2[i]   = s
        phi2[i+1] = s * sign(a)
    end
    # 2n × (2n-r)
    return S2 * Z2 * Diagonal(phi2)
end



"""
    symplectic_normal_form_pair(A::AbstractMatrix{<:Real})

For a real skew-symmetric matrix `A` return `Q` such that `A = Q Ω Q^T` where
`Q` is an invertible matrix and `Ω` is the pair symplectic form. If `A`
is singular, the decomposition works, but `Q` is no longer invertible.

I should test this function thoroughly to see if the eigenvalues always come
in pairs, especially for singular matrices.

## alternatively, we can implement this using skewchol from
## SkewLinearAgebra.jl
## that method is faster 5x faster for 40x40 matrices, but
## from the paper not sure how stable
## https://etna.ricam.oeaw.ac.at/vol.11.2000/pp85-93.dir/pp85-93.pdf
# using Test, LinearAlgebra
# import SkewLinearAlgebra as sk
# A = randn(Float64,40,40);
# Aa = (A-A')/2
# C = sk.skewchol(Aa)
# Omega = jc.symplectic_form_pair(size(A,1)÷2)
# # undo the pivot and transpose to account for definition differences
# # skewchol is defined such that transpose(C.R) * C.J * C.R ≈ A[C.p,C.p]
# # I want C.R*C.J*transpose(C.R) = A
# R = transpose(C.R[:, invperm(C.p)])
# @test isapprox(Aa,R*Omega*R')

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

    # use the schur decomposition of A to compute the symplectric normal form.
    F = schur(A)
    T = Matrix(F.T)
    Z = Matrix(F.Z)

    # loop through the blocks
    d = Vector{eltype(T)}(undef, n)
    for i in 1:2:n
        a = (T[i, i+1] - T[i+1, i]) / 2
        s = sqrt(abs(a))
        d[i] = s
        d[i+1] = sign(a) * s
    end

    Q = Z * Diagonal(d)
    return Q
end

function symplectic_normal_form_block(A::AbstractMatrix{<:Real})
    Q = symplectic_normal_form_pair(block_to_pair(A))
    return pair_to_block(Q)
end

function inv_symplectic_pair(S)
    Omega = symplectic_form_pair(size(S,2)÷2)
    return transpose(Omega)*transpose(S)*Omega
end

function inv_symplectic_block(S)
    Omega = symplectic_form_block(size(S,2)÷2)
    return transpose(Omega)*transpose(S)*Omega
end

function inv_bogoliubov_pair(S)
    return inv_symplectic_pair(S)
end

function inv_bogoliubov_block(S)
    return inv_symplectic_block(S)
end


"""
    autonne_takagi(M::AbstractMatrix)

Return a vector `Λ` and a unitary matrix `W` for a symmetric complex input
matrix `M` such that `M == W*Diagonal(Λ)*transpose(W)` where `M` satisfies
`M = transpose(M)`. Note that if `M` complex this means `M` is not Hermitian.

# References
[1] A. M. Chebotarev and A. E. Teretenkov, “Singular value decomposition for
the Takagi factorization of symmetric matrices,” Applied Mathematics and
Computation, vol. 234, pp. 380–384, May 2014, doi: 10.1016/j.amc.2014.01.170.
[2] M. Houde, W. McCutcheon, and N. Quesada, “Matrix decompositions in Quantum
Optics: Takagi/Autonne, Bloch-Messiah/Euler, Iwasawa, and Williamson,” Can.
J. Phys., vol. 102, no. 10, pp. 497–507, Oct. 2024, doi: 10.1139/cjp-2024-0070.
[3] https://github.com/XanaduAI/thewalrus/pull/403
"""
function autonne_takagi(M::AbstractMatrix)

    # this test is very strict and checks for equality not
    # just approximate quality. unclear if that's a good thing.
    if !issymmetric(M)
        error(lazy"M must be symmetric.")
    end
    # does svd specialize for a symmetric matrix?
    # F = svd(Symmetric(M))
    # unfortunately for Julia 1.10 and 1.11 svd doesn't have a method defined
    # for a symmetric complex matrix.
    # i could use a try catch block or look for the method so i can use
    # svd(Symmetric(M)) on more recent Julia versions.
    F = svd(Matrix(M))
    # this is how Chebotarev and Teretenkov 2014 define Z
    Z = F.U' * transpose(F.Vt)

    E = schur(Z)
    optimum_angle, optimum_rotation = optimum_eigenvalue_angle(E.values)
    shift = exp(im * optimum_rotation)
    invshifto2 = exp(-im * optimum_rotation / 2)
    Zsqrt = invshifto2 * E.vectors * sqrt(shift * E.Schur) * E.vectors'
    W = F.U * Zsqrt
    return F.S, W
end

# function sqrtm(Z,F;atol=eps(eltype(Z)))
#     Zsqrt = zeros(eltype(Z),size(Z))
#     imin = 1
#     imax = 1
#     for i in 2:size(Z,1)
#         # if the current singular value is within machine
#         # precision of the next one. add to the current
#         # block, otherwise take the sqrt of the current block
#         if abs(F.S[i]-F.S[i-1]) > atol
#             # if this element is degenerate with the previous one
#             # increment the max index of the block
#             # Zsqrt[imin:imax,imin:imax] .= sqrt(Z[imin:imax,imin:imax])
#             if imin == imax
#                 Zsqrt[imin:imax,imin:imax] = sqrt(Z[imin:imax,imin:imax])
#             else
#                 E = eigen(Z[imin:imax,imin:imax])
#                 Zsqrt[imin:imax,imin:imax] .= E.vectors*Diagonal(sqrt.(Complex.(E.values)))*E.vectors'
#             end
#             imin = i
#         end
#         imax = i
#     end
#     # do the last square root
#     # Zsqrt[imin:imax,imin:imax] .= sqrt(Z[imin:imax,imin:imax])
#     if imin == imax
#         Zsqrt[imin:imax,imin:imax] = sqrt(Z[imin:imax,imin:imax])
#     else
#         E = eigen(Z[imin:imax,imin:imax])
#         Zsqrt[imin:imax,imin:imax] .= E.vectors*Diagonal(sqrt.(Complex.(E.values)))*E.vectors'
#     end

#     return Zsqrt
# end

"""
    autonne_takagi(M::AbstractMatrix{<:Real})

Return a vector `Λ` and a unitary matrix `W` for input matrix `M` such that
`M == W*Diagonal(Λ)*transpose(W)` where `M` is a symmetric real matrix
`M = transpose(M)`.

"""
function autonne_takagi(M::AbstractMatrix{<:Real})

    # this test is very strict and checks for equality not
    # just approximate quality. unclear if that's a good thing.
    if !issymmetric(M)
        error(lazy"M must be symmetric.")
    end
    F = eigen(Symmetric(M); sortby=abs)

    # return abs.(F.values),F.vectors*Diagonal(sqrt.(sign.(Complex.(F.values))))

    # not sure if this is necessary for this one.
    optimum_angle, optimum_rotation = optimum_eigenvalue_angle(F.values)
    shift = exp(im * optimum_rotation)
    invshifto2 = exp(-im * optimum_rotation / 2)
    return (Λ=abs.(F.values), M=invshifto2 * F.vectors * Diagonal(sqrt.(sign.(shift * F.values))))

end


"""
    bloch_messiah_block(S::AbstractMatrix{<:Real})

Return the Bloch-Messiah (Euler) decomposition  `O`, `D`, `Q` of the
symplectic matrix `S = O*Diagonal(D)*Q` where `O` and `Q` are
orthogonal-symplectic matrices and `Diagonal(D)` is a symplectic-diagonal and
positive definite matrix. This is also called the symplectic singular value
decomposition (SVD). The matrices are symplectric with respect to the block
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

    # check that S is square

    # get the size
    n = size(S, 1) ÷ 2

    # test if S is symplectic
    if !is_symplectic_block(S)
        error(lazy"A must be symplectic.")
    end

    # perform a polar decomposition
    # S = PY
    P, Y = polar(S)

    # partition the symplectic matrix P
    # P = [A B; B^T C]
    A = view(P, 1:n, 1:n)
    B = view(P, 1:n, n+1:2*n)
    Bt = view(P, n+1:2*n, 1:n)
    C = view(P, n+1:2*n, n+1:2*n)

    # compute
    # M = 1/2*(A - C + i(B+B^T)
    M = 1 / 2 * (A .- C .+ im * B .+ im * Bt)

    # perform Takagi-Autonne decomposition
    # M = W*Λ*W^T
    # @show M
    Λ, W = autonne_takagi(Symmetric(M))
    # Λ, W = autonne_takagi(M)

    # form O = [Re(W) -Im(W);Im(W) Re(W)]
    # and Γ = Λ + sqrt(I + Λ^2), D = Γ direct sum 1/Γ
    Γ = Λ .+ sqrt.(1 .+ Λ .^ 2)
    #
    D = vcat(Γ, inv.(Γ))

    # compute O = [Re(W) -Im(W);Im(W) Re(W)]
    O = [real(W) -imag(W); imag(W) real(W)]

    # compute Q = O^T * Y
    Q = O' * Y

    return (O=O, D=D, Q=Q)

end

function bloch_messiah_pair(S::AbstractMatrix{<:Real})
    F = bloch_messiah_block(pair_to_block(S))
    return (O=block_to_pair(F.O), D=block_to_pair(F.D),
        Q=block_to_pair(F.Q))
end

"""
    pre_iwasawa_block(S::AbstractMatrix)

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

    # partition the symplectic matrix S
    # S = [A B; C D]
    A = view(S, 1:n, 1:n)
    B = view(S, 1:n, n+1:2*n)
    C = view(S, n+1:2*n, 1:n)
    D = view(S, n+1:2*n, n+1:2*n)

    #
    A0 = sqrt(A * transpose(A) + B * transpose(B))
    invA0 = inv(A0)
    C0 = (C * transpose(A) + D * transpose(B)) * invA0
    X = invA0 * A
    Y = invA0 * B

    E = [I(n) 0*I(n); C0*invA0 I(n)]
    D = [A0 0*I(n); 0*I(n) invA0]
    F = [X Y; -Y X]

    return (E=E, D=D, F=F)
end

function pre_iwasawa_pair(S::AbstractMatrix)
    F = pre_iwasawa_block(pair_to_block(S))
    return (E=block_to_pair(F.E), D=block_to_pair(F.D),
        F=block_to_pair(F.F))
end

"""
    iwasawa_block(S::AbstractMatrix)

Return the Iwasawa (KAN) decomposition `K`, `A`, `N` of the symplectic matrix
`S`. `K` is a unitary symplectic matrix (maximal compact), `A` is a diagonal
symplectic matrix (Abelian), and `N` is a upper triangular symplectric matrix
(nilpotent). The symplectric matrices are symplectic with respect to the block
symplectric form `Ω`. This decomposition is unique.

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

    # keep D as a vector
    # D = diag(R11) .^ 2
    # switched to abs2 for complex matrices
    D = abs2.(diag(R11))
    Dsqrt = sqrt.(D)
    Dinvsqrt = inv.(Dsqrt)

    # define K, A, N
    A = Diagonal(vcat(Dsqrt, Dinvsqrt))

    K11 = Q11 * H * Diagonal(Dinvsqrt)
    K12 = -Q21 * H * Diagonal(Dinvsqrt)
    # K = [K11 K12; -K12 K11]
    # added conjugates for complex matrices
    K = [K11 conj(K12); -K12 conj(K11)]


    # changed from transpose to adjoint for complex matrices
    # N = inv(A)*transpose(K)*[S12;S22]
    # N = inv(A) * K' * [S12; S22]
    # switch to ldiv
    N = A\(K' * [S12; S22])

    N12 = view(N, 1:n, 1:n)
    N22 = view(N, n+1:2*n, 1:n)

    N = [U N12; 0*I(n) N22]

    return (K=K, A=A, N=N)
end

function iwasawa_pair(S::AbstractMatrix)
    F = iwasawa_block(pair_to_block(S))
    return (K=block_to_pair(F.K), A=block_to_pair(F.A),
        N=block_to_pair(F.N))
end

function iwasawa_bogoliubov_pair(S::AbstractMatrix)
    F = iwasawa_block(ladder_to_quadrature_block(pair_to_block(S)))
    return (
        K=quadrature_to_ladder_pair(block_to_pair(F.K)),
        A=quadrature_to_ladder_pair(block_to_pair(F.A)),
        N=quadrature_to_ladder_pair(block_to_pair(F.N)),
    )
end

function iwasawa_bogoliubov_block(S::AbstractMatrix)
    F = iwasawa_block(ladder_to_quadrature_block(S))
    return (
        K=quadrature_to_ladder_block(F.K),
        A=quadrature_to_ladder_block(F.A),
        N=quadrature_to_ladder_block(F.N),
    )
end


# function B_from_X_Y_pair0(X::AbstractMatrix{<:Real},
#     Y::AbstractMatrix{<:Real})

#     # check that X and Y are the same size

#     # check that the size is even

#     n = size(X,1) ÷ 2
#     Omega = symplectic_form_pair(n)

#     # test if Y is symmetric
#     if !issymmetric(Y)
#         error("`Y` must be positive definite and thus symmetric.")
#     end

#     # test if Y is positive definite. Not all real symmetric matrices are
#     # positive definite. The matrix must be symmetric and the eigenvalues must
#     # all be positive.
#     if !isposdef(Y)
#         error("`Y` must be positive definite.")
#     end

#     vals, vecs = eigen(Symmetric(Y))
#     # since Y is positive definite, all the eigenvalues must be positive
#     # unless there is numerical error. just take the absolute value
#     # before taking the square root to guard against this.
#     Ysqrt = vecs*Diagonal(sqrt.(abs.(vals)))*vecs'
#     Yminussqrt = vecs*Diagonal(inv.(sqrt.(abs.(vals))))*vecs'

#     Delta = Omega - X*Omega*X'
#     # K = Ysqrt*Omega*Ysqrt
#     K = Yminussqrt*Delta*Yminussqrt

#     O = symplectic_normal_form_pair(K)

#     B = Ysqrt*O

#     return B
# end



function B_from_X_Y(Omega::AbstractMatrix, X::AbstractMatrix{<:Real},
    Y::AbstractMatrix{<:Real})

    # check that X and Y are the same size

    # check that the size is even

    # # test if Y is symmetric
    # if !issymmetric(Y)
    #     error("`Y` must be positive definite and thus symmetric.")
    # end

    # test if Y is positive definite. Not all real symmetric matrices are
    # positive definite. The matrix must be symmetric and the eigenvalues must
    # all be positive.
    if !is_positive_semi_definite(Y)
        error(lazy"`Y` must be positive semi-definite.")
    end

    # vals, vecs = eigen(Symmetric(Y))
    # # since Y is positive definite, all the eigenvalues must be positive
    # # unless there is numerical error. just take the absolute value
    # # before taking the square root to guard against this.
    # Ysqrt = vecs*Diagonal(sqrt.(abs.(vals)))*vecs'
    # Yminussqrt = vecs*Diagonal(1.0./sqrt.(abs.(vals)))*vecs'

    Delta = Matrix(Omega - X * Omega * X')

    Gamma = Hermitian(Y + im * Delta)

    # compute the eigenvalues and eigenvectors
    vals, vecs = eigen(Gamma)

    # println(vals)

    # restrict the eigenvalues to greater than zero
    clamp!(vals, 0, Inf)

    F = vecs * Diagonal(sqrt.(vals))

    # this is specific to the block form
    B = [imag(F) real(F)]

    return B
end

"""
    B_from_X_Y(Omega::AbstractMatrix, X::AbstractMatrix{<:Real},
        Y::AbstractMatrix{<:Real})

Return the `B` part of a symplectic matrix `S=[A B;C D]` from the completely
positive trace preserving (CPTP) map `X`, `Y` assuming the environment is in a
vacuum state.

"""
function B_from_X_Y_block(X::AbstractMatrix{<:Real},
    Y::AbstractMatrix{<:Real})

    n = size(X, 1) ÷ 2
    Omega = symplectic_form_block(n)
    B = B_from_X_Y(Omega, X, Y)
    return B
end

function B_from_X_Y_pair(X::AbstractMatrix{<:Real},
    Y::AbstractMatrix{<:Real})

    n = size(X, 1) ÷ 2
    Omega = symplectic_form_pair(n)
    B = B_from_X_Y(Omega, X, Y)

    # permute the columns of B to the pair form
    p = block_to_pair_perm(size(B, 2) ÷ 2)
    return B[:, p]
end

# """

#     X_Y_to_sympletic(X::AbstractMatrix{<:Real}, Y::AbstractMatrix{<:Real})

# Return the symplectic matrix `S` from the completely positive trace preserving
# (CPTP) map `X`, `Y` assumming a vacuum environment.

# """
# function X_Y_to_sympletic(X::AbstractMatrix{<:Real}, Y::AbstractMatrix{<:Real})

# end

"""

    halmos_dilation(S)

Return the Halmos dilation of the passive potentially lossy scattering
parameter matrix `S`.


X, Y = CP_attenuator(0.4,1.0)
U = jc.halmos_dilation(X)

does this assume a block form?
This is only for passive systems so does the operator
not matter? it would make sense if this is just the a's
but it should also work for a's and adag's if it's just
a beamsplitter interaction

I should add the SVD formula

U = [S sqrt(I(size(S,1)) - S*S');sqrt(I(size(S,1)) - S*S') -S']

We can compute this using the SVD


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
function halmos_dilation(S)
    n = size(S, 1)

    # perform the dilation
    # W, V, sigma
    # U = [W 0;0 V]*[sigma sqrt(I-sigma^2);sqrt(I-sigma^2) -sigma]*[V' 0;0 W']
    F = svd(S)
    U = [F.U 0*I(n); 0*I(n) F.V] * [Diagonal(F.S) Diagonal(sqrt.(1.0 .- F.S .^ 2)); Diagonal(sqrt.(1.0 .- F.S .^ 2)) -Diagonal(F.S)] * [F.Vt 0*I(n); 0*I(n) F.U']
    return U
end

function Ymin_from_X(Omega, X; method=1)

    Delta = Matrix(Omega .- X * Omega * transpose(X))
    # Delta = 0.5 * (Delta - transpose(Delta))

    if method == 1
        # method 1: eigen
        # this method should be stable and is fastest. let's go with this as the
        # default
        F = eigen(Hermitian(im * Delta))
        Ymin = real(F.vectors * Diagonal(abs.(F.values)) * F.vectors')

    elseif method == 2
        # method 2: schur
        F = schur(Delta)
        Ymin = F.Z * Diagonal(abs.(imag(F.values))) * transpose(F.Z)

    elseif method == 3
        # method 3: svd
        # are there any issues with this method?
        F = svd(Delta)
        Ymin = F.V * Diagonal(F.S) * F.Vt

    else
        error(lazy"Unknown method")
    end

    # return (Ymin+Ymin')/2
    return Ymin
    # return Symmetric(Ymin)
end

function Ymin_from_X_pair(X; method=1)
    Omega = symplectic_form_pair(size(X, 1) ÷ 2)
    return Ymin_from_X(Omega, X; method=method)
end

function Ymin_from_X_block(X; method=1)
    Omega = symplectic_form_block(size(X, 1) ÷ 2)
    return Ymin_from_X(Omega, X; method=method)
end


# block diagonal of two matrices
function blockdiag(A::AbstractMatrix, B::AbstractMatrix)
    T = promote_type(eltype(A), eltype(B))
    m, n = size(A)
    p, q = size(B)
    C = zeros(T, m + p, n + q)
    C[1:m, 1:n] .= A
    C[m+1:end, n+1:end] .= B
    return C
end

"""
    A_B_to_symplectic(A, B; rtol=1e-12, atol=0.0, tol_block=1e-10, check=true)

Given A (2n×2n) and B (2n×4n) such that A*Ω*A' + B*ΩE*B' = Ω,
construct C (4n×2n), D (4n×4n) so that S = [A B; C D] is symplectic
with respect to Ωtot = Ω ⊕ ΩE.

Returns (C, D).

"""
function A_B_to_symplectic_pair(A::AbstractMatrix, B::AbstractMatrix;
    atol::Real=0,
    rtol::Real=(min(size(A, 1), size(A, 2)) * eps(real(float(oneunit(eltype(A)))))) * iszero(atol))

    type_out = promote_type(eltype(A), eltype(B))

    # the number of system modes
    n = size(A, 1) ÷ 2

    # the number of environment modes
    # m = 2*n

    # do i always need twice the environment modes?
    # does this symplectic matrix have a vacuum environment?

    # Ω  = omega(n, T)
    # ΩE = omega(2n, T)
    Ω = symplectic_form_pair(n)
    ΩE = symplectic_form_pair(2n)

    Ωtot = blockdiag(Ω, ΩE)

    W = hcat(A, B)                      # 2n × 6n

    # find complement rows N such that W*Ωtot*N' = 0
    # @show W*Ωtot
    K = nullspace(W * Ωtot; atol=atol, rtol=rtol)  # 6n × k

    if size(K, 2) != 4n
        @warn lazy"Expected nullspace dimension 4n=$(4n), got $(size(K,2)). Try adjusting rtol/atol."
    end
    N = K'
    # k × 6n  (rows span complement)

    # Induced skew form on the complement
    G = N * Ωtot * N'
    # enforce skew-symmetry numerically
    # G = 0.5 * (G - G')

    # check if the matrix A is skew-symmetric
    if !isapprox(G, -transpose(G))
        error(lazy"G must be skew-symmetric.")
    end

    # Real Schur: G = Q*Tschur*Q'
    F = schur(G)
    # for real G, returns real quasi-triangular Schur form
    Q = F.Z
    T = F.T

    # loop through the blocks
    d = Vector{eltype(T)}(undef, 4 * n)
    for i in 1:2:4*n
        a = (T[i, i+1] - T[i+1, i]) / 2
        s = inv(sqrt(abs(a)))
        d[i] = s
        d[i+1] = sign(a) * s
    end

    # should i detect zeros and set them to zero? doesn't work
    # d = inv.(sqrt.(abs.(imag(F.values)))).*sign.(imag(F.values))
    Sscale = Diagonal(d)

    # the completed bottom block-row
    Wperp = Sscale * (Q') * N           # 4n × 6n

    C = Wperp[:, 1:2n]
    D = Wperp[:, 2n+1:end]

    return [A B; C D]
end
