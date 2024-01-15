function rootFinder(
        f::Function,
        x0::Vector;
        tol=1e-9,
        maxIter=500,
        dxMax=[Inf for _ in 1:length(x0)],
        verbose=true,
        useFiniteDiff=false,
        )

    # Compute Jacobian and initialize variables
    x = copy(x0)
    lx = length(x)
    Jx = zeros(lx, lx)
    fx = zeros(lx)
    dx = zeros(lx)

    # Start iterations
    for iter in 1:maxIter

        # Check function value
        fx .= f(x)
        maxf = maximum(abs.(fx))
        if verbose
            @show iter, maxf
        end
        if maxf ≤ tol
            break
        end

        # Compute and trim correction, if needed
        if useFiniteDiff
            Jx .= FiniteDiff.finite_difference_jacobian(f, x)
        else
            Jx .= ForwardDiff.jacobian(f, x)
        end
        dx .= Jx\fx
        for i in eachindex(dx)
            if abs(dx[i]) > dxMax[i]
                dx[i] = sign(dx[i])*dxMax[i]
            end
        end

        # Apply correction
        x .-= dx
    end

    # Return solution
    if verbose
        @show "Solution", x, f(x)
    end
    return x
end
