"""
    x = rootFinder(y, f!, x0)

Finds the solution ``x`` to the nonlinear multivariate root finding problem y = f(x) = 0,
given as input an initial guess ``x0``.

In-place function formulation
f!(y, x)

Author: F. Capolupo\\
European Space Agency, 2024
"""
function rootFinder!(
    y::Vector{T},
    f!::Function,
    x0::Vector{T};
    tol=1e-9,
    maxIter=500,
    dxMax=[Inf for _ in 1:length(x0)],
    verbose=true,
    method=:Broyden,     # :NewtonRaphson, :Broyden, :ModifiedBroyden
    derivatives=:FiniteDiff,   # :ForwardDiff or :FiniteDiff (only for NewtonRaphson)
) where {T}

    # Initialize variables
    x = copy(x0)
    nx = length(x)
    nf = length(y)
    J = zeros(nf, nx)
    H = zeros(nx, nf)
    dx = zeros(nx)
    yOld = copy(y)

    # Start iterations
    for iter in 1:maxIter

        # Check function value
        f!(y, x)
        maxf = maximum(abs, y)
        if verbose
            @show iter, maxf
        end
        if maxf ≤ tol
            break
        end

        # Compute and trim correction, if needed
        if method == :NewtonRaphson
            # Newton-Raphson method
            if derivatives == :FiniteDiff
                # Numerical finite differences
                FiniteDiff.finite_difference_jacobian!(J, f!, x)
            else
                # Autodiff with ForwardDiff
                ForwardDiff.jacobian!(J, f!, x)
            end
            dx .= -J\y

        elseif method == :Broyden
            # Broyden method
            # https://en.wikipedia.org/wiki/Broyden%27s_method
            if iter == 1
                FiniteDiff.finite_difference_jacobian!(J, f!, x)
            else
                J .+= 1/(dx'*dx)*(y - yOld - J*dx)*dx'
            end
            dx .= -J\y
            yOld .= y

        else
            # Modified Broyden method
            # https://documentation.help/GMAT/DifferentialCorrector.html
            # https://en.wikipedia.org/wiki/Broyden%27s_method
            if iter == 1
                FiniteDiff.finite_difference_jacobian!(J, f!, x)
                H .= pinv(J)  # Pinv also works for non square matrices
            else
                H .+= (dx - H*(y - yOld))*(dx'*H/(dx'*H*(y - yOld)))
            end
            dx .= -H*y
            yOld .= y
        end

        for i in eachindex(dx)
            if abs(dx[i]) > dxMax[i]
                dx[i] = sign(dx[i])*dxMax[i]
            end
            x[i] += dx[i]   # Apply correction
        end
    end

    # Return solution
    if verbose
        println("Solution: x = $x, f(x) = $(f(x))")
    end
    return x
end
