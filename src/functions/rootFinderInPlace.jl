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
    f!::F,
    x0::Vector{T};
    tol=1e-9,
    maxIter=500,
    dxMax=fill(Inf, length(x0)),
    verbose=true,
    method=:Broyden,     # :NewtonRaphson, :Broyden, :ModifiedBroyden
    derivatives=:FiniteDiff,   # :ForwardDiff or :FiniteDiff
) where {T, F}

    # Initialize variables
    y .= zero(eltype(y))
    x = copy(x0)
    nx = length(x)
    nf = length(y)
    J = zeros(nf, nx)
    H = zeros(nx, nf)
    dx = zeros(nx)
    yOld = copy(y)
    yJac = copy(y)

    jac!(J, f!, y, x) = derivatives == :FiniteDiff ?
        FiniteDiff.finite_difference_jacobian!(J, f!, x) :
        ForwardDiff.jacobian!(J, f!, y, x)

    # Start iterations
    @inbounds for iter in 1:maxIter

        # Check function value
        f!(y, x)
        maxf = maximum(abs, y)
        verbose && (@show iter, maxf, maximum(abs, dx))
        maxf ≤ tol && break

        # Compute and trim correction, if needed
        if method == :NewtonRaphson
            # Newton-Raphson method
            jac!(J, f!, yJac, x)
            dx .= -J \ y

        elseif method == :Broyden
            # Broyden method
            # https://en.wikipedia.org/wiki/Broyden%27s_method
            if iter == 1
                jac!(J, f!, yJac, x)
            else
                J .+= 1 / dot(dx, dx) * (y - yOld - J * dx) * dx'
            end
            dx .= -J\y
            yOld .= y

        elseif method == :ModifiedBroyden
            # Modified Broyden method
            # https://documentation.help/GMAT/DifferentialCorrector.html
            # https://en.wikipedia.org/wiki/Broyden%27s_method
            if iter == 1
                jac!(J, f!, yJac, x)
                H .= pinv(J)  # Pinv also works for non square matrices
            else
                H .+= (dx - H * (y - yOld)) * (dx' * H / dot(dx, H * (y - yOld)))
            end
            dx .= -H*y
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
        f!(y, x)
        println("Solution: x = $x, f(x) = $y")
    end
    return x
end
