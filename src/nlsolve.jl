

function nlsolve!(fj!, F, J, x; iterations=1000, ftol=1e-8,
    switchofflinesearchtol = 1e-5)

    factorization = KLU.klu(J)

    deltax = copy(x)

    Nsamples = 100
    samples = Float64[]
    fmin = Float64[]
    fvals = Float64[]
    fpvals = Float64[]
    dfdalphavals = Float64[]
    alphas = Float64[]
    normF = Float64[]

    # perform Newton's method with linesearch based on Nocedal and Wright
    # chapter 3 section 5.
    for n = 1:iterations

        # update the residual function and the Jacobian
        fj!(F, J, x)

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

        # evaluate the residual function but not the Jacobian
        # calcfj3!(F, nothing, x, wmodesm, wmodes2m, Rbnm, Rbnmt, invLnm,
        #     Cnm, Gnm, bnm, Ljb, Ljbm, Nmodes,
        #     Nbranches, Lmean, AoLjbmvector, AoLjbm,
        #     AoLjnmindexmap, invLnmindexmap, Gnmindexmap, Cnmindexmap,
        #     AoLjbmindices, conjindicessorted,
        #     freqindexmap, conjsourceindices, conjtargetindices, phimatrix,
        #     AoLjnm, xbAoLjnm, AoLjbmRbnm, xbAoLjbmRbnm,
        #     phimatrixtd, irfftplan, rfftplan,
        # )
        fj!(F, nothing, x)

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
            alpha1 = 1
        end
        # if the fitted alpha overshoots the size of the interval (from 0 to 1),
        # then set alpha to 1 and make a full length step. 
        if alpha1 > 1 || alpha1 <= 0
            alpha1 = 1
            f1fit = fp
        end

        # switch to newton once the norm is small enough
        # switchofflinesearchtol = 1e-5
        if fp <= switchofflinesearchtol && f <= switchofflinesearchtol && f1fit <= switchofflinesearchtol
            alpha1 = 1
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