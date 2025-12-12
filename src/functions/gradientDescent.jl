function gradientDescent(f, x0; maxIter=300, dxMax=Inf*ones(length(x0)), tol=1e-6, η=1.0, derivatives=:FiniteDiff, method=:grad, verbose=true)
    x = copy(x0)
    nx = length(x)
    δx = zeros(nx)
    ∇f = zeros(nx)
    H = zeros(nx, nx)
    mx = zeros(nx)
    vx = zeros(nx)

    for iter in 1:maxIter
        if derivatives == :ForwardDiff
            ∇f .= ForwardDiff.gradient(f, x)
        else
            ∇f .= FiniteDiff.finite_difference_gradient(f, x)
        end

        if method == :rmsprop
            mx .= 0.9*mx + 0.1*∇f .^ 2
            δx .= -η ./ (sqrt.(mx) .+ 1e-8) .* ∇f

        elseif method == :momentum
            mx .= 0.3*mx + η*∇f
            δx .= -mx

        elseif method == :adam
            # Adam gradient update
            # [1] Ruder, An overview of Gradient Descent Optimization Algorithms
            # [2] https://arxiv.org/pdf/1412.6980.pdf
            β1 = 0.9
            β2 = 0.999
            b1 = 1 - β1^iter
            b2 = 1 - β2^iter
            mx .= β1*mx + (1 - β1)*∇f
            vx .= β2*vx + (1 - β2)*∇f .^ 2
            δx .= -η ./ (sqrt.(vx ./ b2) .+ 1e-8) .* mx ./ b1

        elseif method == :gradls
            # Gradient descent with line search
            ηOpt = goldenSectionSearch(n -> f(x - n*∇f), 0.0, 10η; verbose=false)
            δx .= -ηOpt*∇f

        elseif method == :newton    # {CHECK: IS THIS WORKING?}
            if derivatives == :ForwardDiff
                H .= ForwardDiff.hessian(f, x)
            else
                H .= FiniteDiff.finite_difference_hessian(f, x)
            end
            δx .= -H\∇f

        else
            # Gradient method
            δx .= -η*∇f
        end

        for i in eachindex(δx)
            if abs(δx[i]) > dxMax[i]
                δx[i] = sign(δx[i])*dxMax[i]
            end
        end

        x .+= δx
        err = norm(δx)
        if verbose
            @show iter, f(x)
        end
        if err ≤ tol || any(isnan.(x))
            break
        end
    end

    return x
end
