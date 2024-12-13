"""
    x, P = lsq(r, x0)

Nonlinear Least Squares
This function computes the parameters 'x' that minimize the weighted norm of a nonlinear
multivariate function r(x) = [r₁(x); r₂(x); ...; rₙ(x)].

minₓ J = 1/2 rᵀ(x) W r(x)

Alternatively, the residuals function can be provided as a vector of functions, and the
equivalent problem to be minimized becomes

minₓ J = 1/2 ∑ᵢ rᵢᵀ(x) Wᵢ rᵢ(x)

where each rᵢ can be **the same** multivariate function evaluated at the i-th point.
Caution: rᵢ shall return a vector, so, if the output is a scalar, rᵢ shall return [rᵢ].

NB: usually in estimation problems r(x) = h(x) - y, where y are stacked noisy measurements
and h(x) is the nonlinear measurement model.

Inputs:
r         Residuals function handle, r(x) = [r₁; r₂; ...], or r = [r₁(x), r₂(x), ...]
x0        Initial guess for x

Outputs:
x         Solution to the least-squares problem
P         Covariance matrix of the solution (only meaningful when W is provided)

Author: F. Capolupo
European Space Agency, 2020
"""
function lsq(
        r,                          # Residual function r(x)
        x0;                         # Initial guess
        lb=-Inf,                    # Lower bound on x [can be a vector]
        ub=+Inf,                    # Upper bound on x [can be a vector]
        algorithm=:lm,              # Algorithm choice, can be :lm (Levemberg-Marquardt), or :grad (classic gradient descent)
        maxIter=1000,               # Maximum number of iterations
        dxMax=+Inf,                 # Maximum correction step amplitude
        tol=1e-9,                   # Tolerance on |dx| to declare convergence
        tolRes=1e-9,                # Tolerance on the resiudals, J = r'*W*r
        λ=1e-3,                     # Levemberg-Marquardt parameter
        maxIterStuck=10,            # Number of iterations to declare the algorithm stuck
        relTolStuck=0.1/100,        # Relative tolerance to declare the algorithm stuck
        verbose=true,               # Show progress
        userJacobian=false,         # input function does provide jacobian, i.e., r, H = r(x), where H = ∂r/∂x
        W=(r) -> ones(length(r)),   # Residuals weighting factor - shall return a vector
    )

    if algorithm != :lm; λ = 0.0; end

    x = copy(x0)
    δx = similar(x)
    isFun = isfun(r)
    nx = length(x0)
    apply(x, δx) = min.(max.(x + δx, lb), ub)
    iter = 0; Jold = 0.0; iStuck = 0
    while iter < maxIter

        # Evaluate residual, jacobian, and additional auxiliary matrices
        J, HᵀWH, HᵀWy = lsq_getJacobians(r, x, W, userJacobian, isFun, nx)

        # Check convergence criterion
        if J ≤ tolRes
            if verbose; println("Solution found: residual below tolerance"); end
            break
        end

        # Check if it is stuck
        if iter > 0
            if (Jold - J)/J ≥ relTolStuck
                iStuck = 0.0
            else
                iStuck += 1
                if iStuck > maxIterStuck
                    if verbose; println("Solver stalled"); end
                    break
                end
            end
        end
        Jold = J

        # Compute correction: Levemberg-Marquardt method
        @label iterLS
        iter += 1
        #δx .= -(HᵀWH + λ*diagm(diag(HᵀWH)))\HᵀWy
        δx .= -(HᵀWH + λ*I)\HᵀWy

        # Saturate correction
        dxm = maximum(abs.(δx))
        if dxm > dxMax
            δx .*= dxMax/dxm    #avoid overshoots
        end

        # Check convergence criterion
        err = norm(δx)
        if verbose; @printf("Iteration: %s    δx: %.3e    λ: %.3e    res: %.3e\n", rpad("$iter", 5, " "), err, λ, J); end

        if err ≤ tol
            if verbose; println("Solution found: correction norm below tolerance"); end
            x .= apply(x, δx)   # Correct the estimate before exiting the loop
            break
        end

        # Adapt LM parameter
        if λ > 0.0
            Jtest = lsq_getResiduals(r, apply(x, δx), W, userJacobian, isFun)
            if Jtest < J
                λ /= 10.0
            else
                # In this case the correction term is not accepted, and a new correction
                # term with increased λ is computed, until Jtest decreases.
                λ *= 10.0
                @goto iterLS
            end
        end

        # Correct Estimate
        x .= apply(x, δx)     # Saturate bounds
    end

    if iter == maxIter && verbose
        println("Max number of iterations reached")
    end

    ~, HᵀWH, ~ = lsq_getJacobians(r, x, W, userJacobian, isFun, nx)
    return x, inv(HᵀWH)
end

isfun(f) = !isempty(methods(f))

function lsq_getResiduals(r, x, W, userJacobian, isFun)
    if isFun
        # r(x) is a function
        res = userJacobian ? r(x)[1] : r(x)
        Wk = W(res)
        J = res'*(Wk.*res)
        return J/2
    end

    # r(x) is a vector of functions
    J = 0.0
    for k in eachindex(r)
        res = userJacobian ? r[k](x)[1] : r[k](x)
        Wk = W(res)
        J += res'*(Wk.*res)
    end
    return J/2
end

function lsq_getJacobians(r, x, W, userJacobian, isFun, nx)
    if isFun
        # r(x) is a function
        res, H = userJacobian ? r(x) : r(x), ForwardDiff.jacobian(r, x)
        Wk = W(res)
        J = res'*(Wk.*res)
        HᵀW = H'*diagm(Wk)
        HᵀWH = HᵀW*H
        HᵀWy = HᵀW*res
        return J/2, HᵀWH, HᵀWy
    end

    # r(x) is a vector of functions
    res, H = userJacobian ? r[1](x) : r[1](x), ForwardDiff.jacobian(r[1], x)
    J = 0.0; HᵀWH = zeros(nx, nx)
    HᵀWy = zeros(nx, length(res))
    for k in eachindex(r)
        if k > 1
            res, H = userJacobian ? r[k](x) : r[k](x), ForwardDiff.jacobian(r[k], x)
        end
        Wk = W(res)
        J += res'*(Wk.*res)
        HᵀW = H'*diagm(Wk)
        HᵀWH += HᵀW*H
        HᵀWy += HᵀW*res
    end
    return J/2, HᵀWH, HᵀWy
end

# The following function provides a robust and generic lsq weighting function to be used
# when the lsq data is corrupted by outliers.
# c must be > 0.0
#
# References:
# [1] Barron, General and Adaptive Robust Loss Function
# [2] Chebrolu, Labe, Vysotska, Behley, Stachniss, Adaptive Robust Kernels for Non-Linear
#     Least Squares problems
# [3] Zhang, Parameter Estimation Techniques: A Tutorial with Application to Conic Fitting
lsqWeight(x, α, c) =  ρ′.(x, α, c)./x

function ρ′(x, α, c)
    if abs(α - 2.0) < 1e-8; return x/c^2;        end
    if abs(α) < 1e-8; return 2x/(x^2 + 2c^2);    end
    if α == -Inf; return x/c^2*exp(-0.5(x/c)^2); end
    return x/c^2*((x/c)^2/abs(α - 2) + 1)^(α/2 - 1)
end

# ~~~ LS TEST CODE ~~~ #

# using JTools
# using GLMakie

# function main()
#     t = collect(range(0, 3, 100))
#     xTrue = [1.3; 0.3]
#     xEst0 = [4.0; 1.0]

#     fun(t, x) = 1/x[2]^3*exp(-t*x[1])
#     y = fun.(t, Ref(xTrue)) + 0.5*randn(size(t))
#     W(x) = lsqWeight(x, 0.0, 5.0)

#     # Add outliers
#     y[50] += 50
#     y[70] -= 30

#     f(x) = fun.(t, Ref(x)) - y
#     @time xEst1, ~ = lsq(f, xEst0; W=W)
#     res = [x -> [fun(t[k], x) - y[k]] for k in eachindex(t)]
#     @time xEst2, ~ = lsq(res, xEst0; W=W)

#     @show xEst1
#     @show xEst2

#     fig, ax = scatter(t, y); display(fig)
#     lines!(ax, t, fun.(t, Ref(xTrue)); linewidth=3, color=:black)
#     lines!(ax, t, fun.(t, Ref(xEst1)); linewidth=3, color=:green)
#     lines!(ax, t, fun.(t, Ref(xEst2)); linewidth=3, color=:red, linestyle=:dash)
# end
# main()


# # fun(x) = [10.0*(x[2] - x[1]^2); 1.0 - x[1]]
# # x0 = [-2.0; 2.0]
# # xEst = lsq(fun, x0)[1]
