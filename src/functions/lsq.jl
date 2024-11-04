"""
    x, P = lsq(f, x0)

Nonlinear Least Squares
This function computes the parameters 'x' that minimize the weighted norm of a nonlinear
multivariate function f(x) = [f₁(x); f₂(x); ...; fₙ(x)].

minₓ: J = fᵀ(x) W f(x)

NB: usually in estimation problems f(x) = h(x) - y, where y are stacked noisy measurements
and h(x) is the nonlinear measurement model.

Inputs:
f         Residuals function handle
x0        Initial guess for x

Outputs:
x         Solution to the least-squares problem
P         Covariance matrix of the solution (only meaningful when W is provided)

Author: F. Capolupo
European Space Agency, 2020
"""
function lsq(
        f,                          # Residual function f(x)
        x0;                         # Initial guess
        lb=-Inf,                    # Lower bound on x [can be a vector]
        ub=+Inf,                    # Upper bound on x [can be a vector]
        algorithm=:lm,              # Algorithm choice, can be :lm (Levemberg-Marquardt), or :grad (classic gradient descent)
        maxIter=1000,               # Maximum number of iterations
        dxMax=+Inf,                 # Maximum correction step amplitude
        tol=1e-9,                   # Tolerance on |dx| to declare convergence
        tolRes=1e-9,                # Tolerance on the resiudals, J = f'*W*f
        λ=1e-3,                     # Levemberg-Marquardt parameter
        maxIterStuck=10,            # Number of iterations to declare the algorithm stuck
        relTolStuck=0.1/100,        # Relative tolerance to declare the algorithm stuck
        verbose=true,               # Show progress
        userJacobian=false,         # input function does provide jacobian, i.e., r, J = f(x)
        W=(x, y) -> 1.0,            # Residuals weighting matrix (i.e., inv(R))
    )

    if algorithm != :lm; λ = 0.0; end
    iStuck = 0
    x = copy(x0)
    δx = similar(x)

    fun(x) = userJacobian ? f(x) : f(x), ForwardDiff.jacobian(f, x)
    res(x) = userJacobian ? f(x)[1] : f(x)      # This saves some computation when the user does not provide
    apply(x, δx) = min.(max.(x + δx, lb), ub)
    iter = 0; Vold = 0.0
    while iter < maxIter

        # Evaluate residual and jacobian
        y, J = fun(x)
        Wk = W(x, y)
        V = y'*Wk*y

        # Check convergence criterion
        if V ≤ tolRes
            if verbose; println("Solution found: residual below tolerance"); end
            break
        end

        # Check if it is stuck
        if iter > 0
            if (Vold - V)/V < relTolStuck
                iStuck += 1
                if iStuck > maxIterStuck
                    if verbose; println("Solver stalled"); end
                    break
                end
            else
                iStuck = 0
            end
        end
        Vold = V

        # Auxiliary matrices
        Jt = J'*Wk
        H = Jt*J

        # Compute correction: Levemberg-Marquardt method
        @label iterLS
        iter += 1
        #δx .= -(H + λ*diagm(diag(H)))\(Jt*y)
        δx .= -(H + λ*I)\(Jt*y)

        # Saturate correction
        dxm = maximum(abs.(δx))
        if dxm > dxMax
            δx .*= dxMax/dxm    #avoid overshoots
        end

        # Check convergence criterion
        err = norm(δx)
        if verbose; @printf("Iteration: %d    δx: %.3e    λ: %.3e    res: %.3e\n", iter, err, λ, V); end

        if err ≤ tol
            if verbose; println("Solution found: correction norm below tolerance"); end
            x .= apply(x, δx)   # Correct the estimate before exiting the loop
            break
        end

        # Adapt LM parameter
        if λ > 0.0
            ytest = res(apply(x, δx))
            Vtest = ytest'*Wk*ytest
            if Vtest < V
                λ /= 10.0
            else
                # In this case the correction term is not accepted, and a new correction
                # term with increased λ is computed, until Vtest decreases.
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

    y, J = fun(x)
    return x, inv(J'*W(x, y)*J)
end


# ~~~ LS TEST CODE ~~~ #
#=
using JTools
using GLMakie

d = collect(range(0, 3, 100))
xTrue = [1.3; 0.3]
xEst0 = [4.0; 1.0]
fun(x) = 1/x[2]^3*exp.(-d*x[1])
y = fun(xTrue) + 0.5*randn(size(d))
f(x) = fun(x) - y

xEst, P = lsq(f, xEst0)

@show xEst

fig, ax = scatter(d, y)
lines!(ax, d, fun(xEst))
display(fig)
=#
