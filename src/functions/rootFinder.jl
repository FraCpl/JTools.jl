"""
    x = rootFinder(f, x0)

Finds the solution ``x`` to the nonlinear multivariate root finding problem f(x) = 0,
given as input an initial guess ``x0``.

Author: F. Capolupo\\
European Space Agency, 2024
"""
function rootFinder(
    f::F,
    x0::Vector{T};
    tol=1e-9,
    maxIter=500,
    dxMax=fill(Inf, length(x0)),
    verbose=true,
    method=:Broyden,     # :NewtonRaphson, :Broyden, :ModifiedBroyden
    derivatives=:FiniteDiff,   # :ForwardDiff or :FiniteDiff
) where {T, F}

    # Initialize variables
    x = copy(x0)
    y = f(x)
    nx = length(x)
    nf = length(y)
    J = zeros(nf, nx)
    H = zeros(nx, nf)
    dx = zeros(nx)
    yOld = copy(y)

    jac = derivatives == :FiniteDiff ?
        FiniteDiff.finite_difference_jacobian :
        ForwardDiff.jacobian

    # Start iterations
    for iter in 1:maxIter

        # Check function value
        if iter > 1
            y .= f(x)
        end
        maxf = maximum(abs, y)
        verbose && (@show iter, maxf, maximum(abs, dx))
        maxf ≤ tol && break

        # Compute and trim correction, if needed
        if method == :NewtonRaphson
            # Newton-Raphson method
            J .= jac(f, x)
            dx .= -J \ y

        elseif method == :Broyden
            # Broyden method
            # https://en.wikipedia.org/wiki/Broyden%27s_method
            if iter == 1
                J .= jac(f, x)
            else
                J .+= 1 / dot(dx, dx) * (y - yOld - J * dx) * dx'
            end
            dx .= -J \ y
            yOld .= y

        elseif method == :ModifiedBroyden
            # Modified Broyden method
            # https://documentation.help/GMAT/DifferentialCorrector.html
            # https://en.wikipedia.org/wiki/Broyden%27s_method
            if iter == 1
                H .= pinv(jac(f, x))  # Pinv also works for non square matrices
            else
                H .+= (dx - H * (y - yOld)) * (dx' * H / dot(dx, H*(y - yOld)))
            end
            dx .= -H * y
            yOld .= y
        else
            error("Unrecognized method $method")
        end

        kMax = 1.0
        for i in eachindex(dx)
            kMax = max(kMax, abs(dx[i]) / dxMax[i])
        end

        for i in eachindex(dx)
            dxApplied = dx[i] / kMax#clamp(dx[i], -dxMax[i], dxMax[i])
            x[i] += dxApplied   # Apply correction
            dx[i] = dxApplied   # TODO: Check
        end
    end

    # Return solution
    if verbose
        println("Solution: x = $x, f(x) = $(f(x))")
    end
    return x
end
