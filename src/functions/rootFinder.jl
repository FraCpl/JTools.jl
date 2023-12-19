function rootFinder(f::Function, x0::Vector; tol=1e-9, maxIter=500, δxMax=[Inf for _ in 1:length(x0)], verbose=true)

    # Compute Jacobian and initialize variables
    x = copy(x0)
    lx = length(x)
    Jx = zeros(lx, lx)
    fx = zeros(lx)
    δx = zeros(lx)

    # Start iterations
    for iter in 1:maxIter

        # Check function value
        fx[:] = f(x)
        maxf = maximum(abs.(fx))
        if verbose
            @show iter, maxf
        end
        if maxf ≤ tol
            break
        end

        # Compute and trim correction, if needed
        ForwardDiff.jacobian!(Jx, f, x)
        δx[:] = Jx\fx
        for i = 1:lx
            if abs(δx[i]) > δxMax[i]
                δx[i] = sign(δx[i])*δxMax[i]
            end
        end

        # Apply correction
        x[:] -= δx
    end

    # Return solution
    if verbose
        @show "Solution", x, f(x)
    end
    return x
end
