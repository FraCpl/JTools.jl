"""
    x = rootFinder(f, x0)

Finds the solution ``x`` to the nonlinear multivariate root finding problem f(x) = 0,
given as input an initial guess ``x0``.

Author: F. Capolupo\\
European Space Agency, 2024
"""
function rootFinder(
    f::Function,
    x0::Vector;
    tol=1e-9,
    maxIter=500,
    dxMax=[Inf for _ in 1:length(x0)],
    verbose=true,
    method=:Broyden,     # :NewtonRaphson, :Broyden, :ModifiedBroyden
    derivatives=:FiniteDiff,   # :ForwardDiff or :FiniteDiff (only for NewtonRaphson)
)

    # Initialize variables
    x = copy(x0)
    y = f(x)
    nx = length(x)
    nf = length(y)
    J = zeros(nf, nx)
    H = zeros(nx, nf)
    dx = zeros(nx)
    yOld = copy(y)

    # Start iterations
    for iter in 1:maxIter

        # Check function value
        if iter > 1
            y .= f(x)
        end
        maxf = maximum(abs.(y))
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
                J .= FiniteDiff.finite_difference_jacobian(f, x)
            else
                # Autodiff with ForwardDiff
                J .= ForwardDiff.jacobian(f, x)
            end
            dx .= -J\y

        elseif method == :Broyden
            # Broyden method
            # https://en.wikipedia.org/wiki/Broyden%27s_method
            if iter == 1
                J .= FiniteDiff.finite_difference_jacobian(f, x)
            else
                J .+= 1/(dx'*dx)*(y - yOld - J*dx)*dx'
            end
            dx .= -J\y
            yOld .= copy(y)

        else
            # Modified Broyden method
            # https://documentation.help/GMAT/DifferentialCorrector.html
            # https://en.wikipedia.org/wiki/Broyden%27s_method
            if iter == 1
                H .= pinv(FiniteDiff.finite_difference_jacobian(f, x))  # Pinv also works for non square matrices
            else
                H .+= (dx - H*(y - yOld))*(dx'*H/(dx'*H*(y - yOld)))
            end
            dx .= -H*y
            yOld .= copy(y)
        end

        for i in eachindex(dx)
            if abs(dx[i]) > dxMax[i]
                dx[i] = sign(dx[i])*dxMax[i]
            end
        end

        # Apply correction
        x .+= dx
    end

    # Return solution
    if verbose
        println("Solution: x = $x, f(x) = $(f(x))")
    end
    return x
end
