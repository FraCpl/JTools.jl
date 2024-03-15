"""
    x, P = lsq(f, x0)

Nonlinear Least Squares
This function computes the parameters 'x' that minimize the weighted norm of a nonlinear
multivariate function f(x) = [f₁(x); f₂(x); ...; fₙ(x)].

minₓ: J = fᵀ(x) W f(x)

NB: usually in estimation problems f(x) = y - h(x), where y are noisy measurements and h(x)
is the nonlinear measurement model.

Inputs:
f         Residuals function handle
x0        Initial guess for x

Outputs:
x         Solution to the least-squares problem
P         Covariance matrix of the solution (only meaningful when W is provided)

Author: F. Capolupo
European Space Agency, 2020
"""
function lsq(f, x0;
        W=1.0,                  # Residuals weighting matrix (i.e., inv(R))
        lb=-Inf*abs.(x0),       # Lower bound on x
        ub=Inf*abs.(x0),        # Upper bound on x
        algorithm=:lm,          # Algorithm choice, can be :lm (Levemberg-Marquardt), or :grad (classic gradient descent)
        maxIter=1000,           # Maximum number of iterations
        dxMax=1e3,              # Maximum correction step amplitude
        tol=1e-9,               # Tolerance on |dx| to declare convergence
        λ=1e-3,                 # Levemberg-Marquardt parameter
        maxIterStuck=10,        # Number of iterations to declare the algorithm stuck
        relTolStuck=0.1/100,    # Relative tolerance to declare the algorithm stuck
        verbose=true,           # Show progress
    )

    iStuck = 0
    x = copy(x0)
    δx = similar(x)
    Vold = 0.0; normyold = 0.0
    for iter in 1:maxIter

        # Evaluate residual and jacobian
        y = f(x)
        J = ForwardDiff.jacobian(f, x)

        # Auxiliary matrices
        Jt = J'*W
        H = Jt*J

        # Compute LS correction
        if algorithm == :lm
            # Levemberg-Marquardt method
            δx .= -(H + λ*diagm(diag(H)))\(Jt*y)

            # Adapt parameter
            V = y'*y
            if iter > 1
                if V < Vold
                    # If cost is decreasing, decrease damping parameter (--> Newton)
                    λ /=10.0
                else
                    # If cost is incrasing, increase damping parameter (--> Gradient)
                    λ *= 10.0
                end
            end
            Vold = V
        else
            # Classic method
            δx .= -H\Jt*y
        end

        # Correct Estimate
        dxm = maximum(abs.(δx))
        if dxm > dxMax
            δx .*= dxMax/dxm    #avoid overshoots
        end
        x .= min.(max.(x + δx, lb), ub)     # Saturate bounds

        # Check for convergence
        err = norm(δx)
        if verbose
            println("Iteration: $iter   δx: $err")
        end
        if err < tol
            if verbose
                println("Solution found")
            end
            break
        end

        # Check if it is stuck
        normy = norm(y)
        if (normyold - normy)/normy < relTolStuck
            iStuck += 1
            if iStuck > maxIterStuck
                if verbose
                    println("Solver stalled")
                end
                break
            end
        else
            iStuck = 0
        end
        normyold = normy
    end

    J = ForwardDiff.jacobian(f, x)
    return x, inv(J'*W*J)
end


# ~~~ LS TEST CODE ~~~ #
#=
using JTools
using GLMakie

d = collect(range(0, 3, 100))
xTrue = [1.3]
xEst0 = [4.0]
y = exp.(-xTrue[1]*d) + 0.05*randn(size(d))
f(x) = exp.(-d*x[1]) - y

x, P = lsq(f, xEst0)#; algorithm=:grad)

fig = Figure()
ax = Axis(fig[1, 1])
scatter!(ax, d, y)
lines!(ax, d, exp.(-d*x[1]))
display(fig)
=#
