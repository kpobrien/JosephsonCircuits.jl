

"""
    Factorization(outofplace,inplace,kwargs)

A structure to hold the factorizations and their keyword arguments.

```
"""
struct Factorization
    factorize
    factorize!
    kwargs
end

function KLUfactorization(;kwargs...)
    return Factorization(KLU.klu,KLU.klu!,kwargs)
end

function LUfactorization(;kwargs...)
    return Factorization(lu,lu!,kwargs)
end

function QRfactorization(;kwargs...)
    return Factorization(qr,nothing,kwargs)
end


"""
    FactorizationCache(factorization)

A cache for the factorization object.

# Examples
```jldoctest
julia> JosephsonCircuits.FactorizationCache(JosephsonCircuits.KLU.klu(JosephsonCircuits.sparse([1, 2], [1, 2], [1/2, 1/2], 2, 2)));

```
"""
mutable struct FactorizationCache
    factorization
end

function FactorizationCache()
    return FactorizationCache(nothing)
end

"""
    tryfactorize!(cache::FactorizationCache,
        factorization::Factorization, A::AbstractArray)

Factorize the matrix `A` using the factorization from `factorization` and
store the result in `cache`. Attempt to reuse the symbolic factorization. Redo the
symbolic factorization if we get a SingularException.

"""
function tryfactorize!(cache::FactorizationCache,
    factorization::Factorization, A::AbstractMatrix)

    # if the factorization cache is empty then generate a factorizaion
    if isnothing(cache.factorization)
        cache.factorization = factorization.factorize(A;
            factorization.kwargs...)
    # otherwise, try to update the factorization, falling back to generating
    # a new one if that fails
    elseif isnothing(factorization.factorize!)
        cache.factorization = factorization.factorize(A;
            factorization.kwargs...)
    else
        try
            # update the factorization. the sparsity structure does 
            # not change so we can reuse the factorization object.
            factorization.factorize!(cache.factorization, A;
                factorization.kwargs...)
        catch e
            if isa(e, SingularException)
                # reusing the symbolic factorization can sometimes
                # lead to numerical problems. if the first linear
                # solve fails try factoring and solving again
                cache.factorization = factorization.factorize(A;
                    factorization.kwargs...)
            else
                throw(e)
            end
        end
    end
    return cache
end


"""
    trysolve!(x,factorization,b)

First try to solve a linear system using ldiv! then if it errors, use \\. The
motivation for this function is some factorizations such as `qr` with sparse
matrices don't support ldiv!. 
"""
function trysolve!(x,factorization,b)
    try
        ldiv!(x,factorization,b)
    catch
        x .= factorization \ b
    end
    return x
end


"""
    linesearch(f, fp, dfdalpha, alphamin)

Quadratic linesearch based on Nocedal and Wright, chapter 3 section 5. `f` is
the value at the first point alpha=0.0, `fp` is the value at the second point,
alpha=1.0, `dfdalpha` is the derivative at the first point, and `alphamin` is
the minimum value of `dfdalpha` below which we will take a full step. The
linesearch will return the fitted minimum of the function with respect to
alpha as (alpha at which minimum occurs, minimum value of function).

"""
function linesearch(f, fp, dfdalpha, alphamin)

    # coefficients of the quadratic equation a*alpha^2+b*alpha+c to interpolate
    # f vs alpha
    a = -dfdalpha + fp - f
    # b = dfdalpha
    # c = f
    alpha1 = -dfdalpha/(2*a)
    f1fit = -dfdalpha*dfdalpha/(4*a) + f

    if isnan(f) || isnan(fp)
        error("NaN in nonlinear solver.")
    end

    if f1fit > fp
        return 1.0, fp
    # if the fitted alpha overshoots the size of the interval (from 0 to 1),
    # then set alpha to 1 and make a full length step.
    elseif alpha1 > 1.0 || alpha1 <= 0.0
        return 1.0, fp
    # if we aren't making sufficient progress, take a step
    # switch to using Armijo rule
    elseif alpha1 <= alphamin
        return 1.0, fp
    # if a is zero, alpha1 will be NaN
    # take a full step
    elseif abs2(a) == 0 
        return 1.0, fp
    else
        return alpha1, f1fit
    end
end

"""
    nlsolve!(fj!::Function, F::AbstractVector{T}, J::AbstractArray{T},
        x::Vector{T}; iterations=1000, ftol=1e-8, switchofflinesearchtol = 1e-5,
        alphamin = 1e-4,factorization = KLUfactorization())

A simple nonlinear solver using Newton's method with linesearch based on
Nocedal and Wright, chapter 3 section 5.

This solver attempts to find x such that f(x) == 0, where f is a
nonlinear function with Jacobian J.

A few points to note:
(1) It uses KLU factorization, so only works on sparse matrices.
(2) The Jacobian J cannot change sparsity structure.
(3) This function attempts to reuse the symbolic factorization which can
    sometimes result in a SingularException, which we catch, then create a
    new factorization object.

# Arguments
- `fj!`: a function to compute a vector-valued objective function and
its Jacobian.
- `F`: matrix for holding intermediate results. Initial values may be
  overwritten and can be bogus values.
- `J`: sparse matrix with with the desired sparsity structure of the
  Jacobian. Initial values may be overwritten and can be bogus values,
  as long as the sparsity structure is correct.
- `x`: initial guess for x.

# Examples
```jldoctest
function fj!(F, J, x)
    if !isnothing(F)
        F[1] = (x[1]+3)*(x[2]^3-7)+18
        F[2] = sin(x[2]*exp(x[1])-1)
    end
    if !isnothing(J)
        J[1, 1] = x[2]^3-7
        J[1, 2] = 3*x[2]^2*(x[1]+3)
        u = exp(x[1])*cos(x[2]*exp(x[1])-1)
        J[2, 1] = x[2]*u
        J[2, 2] = u
    end
    return nothing
end
x = [ 0.1, 1.2]
F = [0.0, 0.0]
J = JosephsonCircuits.sparse([1, 1, 2, 2],[1, 2, 1, 2],[1.3, 0.5, 0.1, 1.2])
JosephsonCircuits.nlsolve!(fj!, F, J, x)
isapprox([0.0,1.0],x)

# output
true
```
"""
function nlsolve!(fj!::Function, F::AbstractVector{T}, J::AbstractArray{T},
    x::Vector{T}; iterations=1000, ftol=1e-8, switchofflinesearchtol = 1e-5,
    alphamin = 1e-4,factorization = KLUfactorization()) where T

    if size(J,1) != size(J,2)
        throw(DimensionMismatch("The Jacobian `J` matrix must be square."))
    end

    if size(J,2) != length(x)
        throw(DimensionMismatch("Second axis of Jacobian `J` must have the same length as the input `x`."))
    end

    if size(J,1) != length(F)
        throw(DimensionMismatch("First axis of the Jacobian `J` must have the same length as the residual `F`."))
    end

    cache = FactorizationCache()
    tryfactorize!(cache,factorization,J)

    deltax = copy(x)

    # Nsamples = 100
    # samples = Float64[]
    fmin = Float64[]
    fvals = Float64[]
    fpvals = Float64[]
    dfdalphavals = Float64[]
    alphas = Float64[]
    normF = Float64[]
    alpha1 = 0.0

    # perform Newton's method with linesearch based on Nocedal and Wright
    # chapter 3 section 5.
    for n in 1:iterations

        if alpha1 == 1.0
            # if alpha was 1, we don't need to update the function 
            # because we have already calculated that in the last
            # loop. just update the jacobian. since we set alpha1=0
            # before the loop, this will never be called on the first
            # iteration.
            fj!(nothing, J, x)
        else
            # update the residual function and the Jacobian
            fj!(F, J, x)
        end

        push!(normF, norm(F))

        # factor the Jacobian
        tryfactorize!(cache,factorization,J)

        # solve the linear system
        trysolve!(deltax, cache.factorization, F)

        # multiply deltax by -1
        rmul!(deltax, -1)

        # calculate the objective function and the derivative of the objective
        # with respect to the scalar variable alpha which parameterizes the
        # path between the old x and the new x. 
        # Note: the dot product takes the complex conjugate of the first vector
        f = real(0.5*dot(F, F))
        dfdalpha = real(dot(F, J, deltax))

        # evaluate the function at the trial point
        fj!(F, nothing, x+deltax)

        fp = real(0.5*dot(F,F))

        # calculate the step size based on the last point, the trial point, and
        # derivative at the first point.
        alpha1, f1fit = linesearch(f,fp,dfdalpha,alphamin)

        # switch to newton once the norm is small enough
        normx = norm(x)
        if normx > 0 && sqrt(fp)/normx <= switchofflinesearchtol && sqrt(f)/normx <= switchofflinesearchtol
            alpha1 = 1.0
            # println("norm(F)/norm(phi): ",sqrt(fp)/norm(x))
        end

        # update x
        x .+= deltax*alpha1
        # push!(alphas,alpha1)

        if norm(F,Inf) <= ftol || ( norm(x) > 0 && norm(F)/norm(x) < ftol)
            # terminate iterations if infinity norm or relative norm are less
            # than ftol. check that norm(x) is greater than zero to avoid
            # divide by zero errors. 
            # println("converged to: infinity norm of : ",norm(F,Inf)," after ",n," iterations")
            # println("norm(F)/norm(phi): ",norm(F)/norm(x))
            break
        end

        if n == iterations
            # Log to warning system
            warning_log("Solver did not converge after maximum iterations of $n")
            warning_log("norm(F)/norm(x): $(norm(F)/norm(x))")
            warning_log("Infinity norm: $(norm(F,Inf))")

            @warn string("Solver did not converge after maximum iterations of ", n,".")
            println("norm(F)/norm(x): ", norm(F)/norm(x))
            println("Infinity norm: ", norm(F,Inf))
            # error(" ")
            # @show alphas
        end
    end
    return nothing
end
