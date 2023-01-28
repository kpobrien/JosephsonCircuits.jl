
"""

    nlsolve!(fj!, F, J, x; iterations=1000, ftol=1e-8,
        switchofflinesearchtol = 1e-5)

A simple nonlinear solver for sparse matrices using Newton's method with
linesearch based on Nocedal and Wright, chapter 3 section 5. A few points to
note:
(1) It uses KLU factorization, so only works on sparse matrices.
(2) The Jacobian J cannot change sparsity structure.
(3) This function attempts to reuse the symbolic factorization which can
    sometimes result in a SingularException, which we catch, then create a
    new factorization object.

# Examples
```jldoctest
function fj!(F, J, x)
    if !(F == nothing)
        F[1] = (x[1]+3)*(x[2]^3-7)+18
        F[2] = sin(x[2]*exp(x[1])-1)
    end
    if !(J == nothing)
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
function nlsolve!(fj!, F, J, x; iterations=1000, ftol=1e-8,
    switchofflinesearchtol = 1e-5)

    factorization = KLU.klu(J)

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

        # # update the residual function and the Jacobian
        # fj!(F, J, x)

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

        # solve the linear system
        try
            # update the factorization. the sparsity structure does 
            # not change so we can reuse the factorization object.
            KLU.klu!(factorization, J)
        catch e
            if isa(e, SingularException)
                # reusing the symbolic factorization can sometimes
                # lead to numerical problems. if the first linear
                # solve fails try factoring and solving again
                factorization = KLU.klu(J)
            else
                throw(e)
            end
        end

        # solve the linear system
        ldiv!(deltax, factorization, F)

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

        # coefficients of the quadratic equation a*alpha^2+b*alpha+c to interpolate
        # f vs alpha
        a = -dfdalpha + fp - f
        b = dfdalpha
        c = f
        alpha1 = -b/(2*a)
        f1fit = -b*b/(4*a) + c

        if f1fit > fp
            f1fit = fp
            alpha1 = 1.0
        end

        # if the fitted alpha overshoots the size of the interval (from 0 to 1),
        # then set alpha to 1 and make a full length step. 
        if alpha1 > 1.0 || alpha1 <= 0
            alpha1 = 1.0
            f1fit = fp
        end

        # if a is zero, alpha1 will be NaN
        if abs2(a) == 0 
            alpha1 = 1.0
        end

        # # switch to newton once the norm is small enough
        # if fp <= switchofflinesearchtol && f <= switchofflinesearchtol && f1fit <= switchofflinesearchtol
        #     alpha1 = 1.0
        # end

        # switch to newton once the norm is small enough
        normx = norm(x)
        if normx > 0 && sqrt(fp)/normx <= switchofflinesearchtol && sqrt(f)/normx <= switchofflinesearchtol
            alpha1 = 1.0
            # println("norm(F)/norm(phi): ",sqrt(fp)/norm(x))
        end

        # update x
        x .+= deltax*alpha1

        if norm(F,Inf) <= ftol || ( norm(x) > 0 && norm(F)/norm(x) < ftol)
            # terminate iterations if infinity norm or relative norm are less
            # than ftol. check that norm(x) is greater than zero to avoid
            # divide by zero errors. 
            # println("converged to: infinity norm of : ",norm(F,Inf)," after ",n," iterations")
            # println("norm(F)/norm(phi): ",norm(F)/norm(x))
            break
        end

        if n == iterations
            @warn string("Solver did not converge after maximum iterations of ", n,".")
            println("norm(F)/norm(x): ", norm(F)/norm(x))
            println("Infinity norm: ", norm(F,Inf))
        end
    end

    return nothing
end