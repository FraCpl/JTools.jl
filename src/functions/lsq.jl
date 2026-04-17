"""
    x, P = lsq(r, x0)

Nonlinear Least Squares
This function computes the parameters 'x' that minimize the weighted norm of a nonlinear
multivariate function r(x) = [r₁(x); r₂(x); ...; rₙ(x)].

minₓ J = 1/2 rᵀ(x) W r(x)

Alternatively, the residuals function can be provided as a vector of functions, and the
equivalent problem to be minimized becomes

minₓ J = 1/2 ∑ᵢ rᵢᵀ(x) Wᵢ rᵢ(x)

where each rᵢ can be a different multivariate function of `x`.
Caution: rᵢ shall return a vector, so, if the output is a scalar, rᵢ shall return [rᵢ].

NB: usually in estimation problems r(x) = h(x) - y, where y are stacked noisy measurements
and h(x) is the nonlinear measurement model.

Inputs:
r         Residuals function handle, r(x) = [r₁; r₂; ...], or r = (r₁(x), r₂(x), ...)
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
    lb=(-Inf),                  # Lower bound on x [can be a vector]
    ub=(+Inf),                  # Upper bound on x [can be a vector]
    maxIter=1000,               # Maximum number of iterations
    dxMax=(+Inf),               # Maximum correction step amplitude
    tol=1e-9,                   # Tolerance on |dx| to declare convergence
    tolRes=1e-9,                # Tolerance on the resiudals, J = r'*W*r
    λ=1e-3,                     # Levemberg-Marquardt parameter (set to 0 to use classic gradient descent)
    ν=2.0,                      # Levemberg-Marquardt update parameter
    maxIterStuck=10,            # Number of iterations to declare the algorithm stuck
    relTolStuck=0.1/100,        # Relative tolerance to declare the algorithm stuck
    verbose=true,               # Show progress
    userJacobian=false,         # input function does provide jacobian, i.e., r, H = r(x), where H = ∂r/∂x
    W=r->ones(length(r)),       # Residuals weighting factor - shall return a vector
)

    x = copy(x0)
    xtest = zero(x0)
    δx = zero(x)
    rfun = isfun(r) ? (r, ) : r   # make r(x) a tuple of functions
    nx = length(x0)
    HᵀWH = zeros(nx, nx)
    H_lm = zeros(nx, nx)
    HᵀWy = zeros(nx)
    iter = 0; iStuck = -1
    Jold = 0.0; J = 0.0
    repeatIter = false
    while iter < maxIter
        if !repeatIter
            # Evaluate residual, jacobian, and additional auxiliary matrices
            J = lsq_getJacobians!(HᵀWH, HᵀWy, rfun, x, W, userJacobian)

            # Check convergence criterion
            lsq_checkConvergence(J, tolRes, verbose) && break

            # Check if it is stuck
            iStuck, stuckFlag = lsq_checkStall(J, Jold, relTolStuck, iStuck, maxIterStuck, verbose)
            stuckFlag && break
            Jold = J
        end

        # Compute correction: Levemberg-Marquardt method
        iter += 1
        repeatIter = false
        H_lm .= HᵀWH
        @inbounds for i in 1:nx
            H_lm[i, i] += λ * HᵀWH[i, i]
        end
        δx .= -H_lm \ HᵀWy

        # Saturate correction
        dxm = maximum(abs, δx)
        if dxm > dxMax
            δx .*= dxMax ./ dxm    #avoid overshoots
        end

        # Check convergence criterion
        verbose && @printf("Iteration: %s    δx: %.3e    res: %.3e    λ: %.3e\n", rpad("$iter", 5, " "), dxm, J, λ)
        if dxm ≤ tol
            verbose && println("Solution found: correction below tolerance")
            apply!(x, δx, lb, ub)   # Correct the estimate before exiting the loop
            break
        end

        # Adapt LM parameter
        if λ > 0
            xtest .= x
            apply!(xtest, δx, lb, ub)
            Jtest = lsq_getResiduals(rfun, xtest, W, userJacobian)

            repeatIter = J < Jtest
            if repeatIter
                # In this case the correction term is not accepted, and a new correction
                # term with increased λ is computed, until Jtest decreases.
                λ *= ν
                continue
            end

            λ /= ν
            x .= xtest  # Accept correction
        else
            # Apply correction
            apply!(x, δx, lb, ub)
        end
    end

    verbose && iter == maxIter && println("Max number of iterations reached")

    lsq_getJacobians!(HᵀWH, HᵀWy, rfun, x, W, userJacobian)      # To update HᵀWH
    res = [(userJacobian ? rk(x)[1] : rk(x)) for rk in rfun]       # Get residuals [overhead: this is already computed in lsq_getJacobians]
    return x, (inv(HᵀWH), res)
end

isfun(f) = isa(f, Function)

function lsq_getResiduals(r, x, W, userJacobian)
    # r(x) is a vector of functions
    J = 0.0
    for rk in r
        res = userJacobian ? rk(x)[1] : rk(x)
        Wk = W(res)
        J += dot(res, (Wk .* res))
    end
    return J / 2
end

function lsq_getJacobians!(HᵀWH, HᵀWy, r, x, W, userJacobian)
    # r(x) is a vector of functions
    J = 0.0
    fill!(HᵀWH, 0.0)
    fill!(HᵀWy, 0.0)
    for rk in r
        if userJacobian
            res, H = rk(x)
        else
            res = rk(x)
            # H = ForwardDiff.jacobian(rk, x)
            H = FiniteDiff.finite_difference_jacobian(rk, x)
        end
        # res, H = userJacobian ? rk(x) : rk(x), ForwardDiff.jacobian(rk, x)
        Wk = W(res)
        J += dot(res, (Wk .* res))
        HᵀW = (Wk .* H)'# H'*diagm(Wk)
        HᵀWH .+= HᵀW*H
        HᵀWy .+= HᵀW*res
    end
    return J / 2
end

function lsq_checkConvergence(J, tolRes, verbose)
    if J ≤ tolRes
        if verbose
            println("Solution found: residual below tolerance")
        end
        return true
    end
    return false
end

function lsq_checkStall(J, Jold, relTolStuck, iStuck, maxIterStuck, verbose)
    flag = false
    iStuck += 1
    if (Jold - J) / J ≥ relTolStuck
        iStuck = -1
    else
        iStuck += 1
        if iStuck > maxIterStuck
            verbose && println("Solver stalled")
            flag = true
        end
    end
    return iStuck, flag
end

# The following function provides a robust and generic lsq weighting function to be used
# when the lsq data is corrupted by outliers.
# c must be > 0.0
#
# References:
# [1] Barron, General and Adaptive Robust Loss Function
#     https://openaccess.thecvf.com/content_CVPR_2019/papers/Barron_A_General_and_Adaptive_Robust_Loss_Function_CVPR_2019_paper.pdf
# [2] Chebrolu, Labe, Vysotska, Behley, Stachniss, Adaptive Robust Kernels for Non-Linear
#     Least Squares problems, https://arxiv.org/pdf/2004.14938
# [3] Zhang, Parameter Estimation Techniques: A Tutorial with Application to Conic Fitting
lsqWeight(x, α, c) = ρ′.(x, α, c) ./ x

function ρ′(x, α, c)
    if abs(α - 2.0) < 1e-8
        return x/c^2
    end
    if abs(α) < 1e-8
        return 2x/(x^2 + 2c^2)
    end
    if α == -Inf
        return x/c^2*exp(-0.5(x/c)^2)
    end
    return x/c^2*((x/c)^2/abs(α - 2) + 1)^(α/2 - 1)
end

function apply!(x::Vector{T}, dx, lb::Vector{T}, ub::Vector{T}) where {T}
    @inbounds for i in eachindex(x)
        x[i] = clamp(x[i] + dx[i], lb[i], ub[i])
    end
    return x
end

function apply!(x::Vector{T}, dx, lb::T, ub::T) where {T}
    @inbounds for i in eachindex(x)
        x[i] = clamp(x[i] + dx[i], lb, ub)
    end
    return x
end
