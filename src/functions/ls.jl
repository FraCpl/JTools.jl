# Nonlinear least squares
function ls(f, x0; lb=-Inf*abs.(x0), ub=Inf*abs.(x0), algorithm=:lm, maxIter=1000, dxMax=1e3, tol=1e-9, λ=1e-3, W=1.0, maxIterStuck=10, relTolStuck=0.1/100)

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
        Hx = Jt*J

        # Compute LS correction
        if algorithm == :lm
            # Levemberg-Marquardt method
            δx .= -(Hx + λ*diagm(diag(Hx)))\(Jt*y)

            # Adapt parameter
            V = y'*y
            if iter > 1
                if V < Vold
                    # If cost is decreasing, decrease damping parameter (--> Newton)
                    λ /=10
                else
                    # If cost is incrasing, increase damping parameter (--> Gradient)
                    λ *= 10
                end
            end
            Vold = V
        else
            # Classic method
            δx .= -Hx\Jt*y
        end

        # Correct Estimate
        dxm = maximum(abs.(δx))
        if dxm > dxMax
            #avoid overshoots
            δx .*= dxMax/dxm
        end

        Main.dbg = x, δx, lb, ub
        #x .= min.(max.(x + δx, lb), ub)     # Saturate bounds
        x .+= δx
        err = norm(δx)
        normy = norm(y)
        println("Iteration: $iter   δx: $err")

        if err < tol
            break
        end

        # Check if iStuck
        if (normyold - normy)/normy < relTolStuck
            iStuck += 1
        else
            iStuck = 0
        end
        normyold = normy
        if iStuck > maxIterStuck
            break
        end
    end

    J = ForwardDiff.jacobian(f, x)
    return x, inv(J'*W*J)
end
