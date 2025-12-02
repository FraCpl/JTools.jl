# Find minimum of unimodal scalar function f(x) of scalar variable x, with x ∈ [lb, ub]
# caution: f(a) must be != than f(b)
function bisection(f::Function, lb, ub; tol = 1e-6, maxIter = 1000, verbose = true)

    # Initialize variables
    a = lb
    b = ub
    fa = f(a)
    if sign(f(b)) == sign(fa)
        if verbose
            println("Ill-posed problem: sign(f(lb)) == sign(f(ub))!")
        end
        return NaN
    end

    # Start iterations
    for iter = 1:maxIter
        c = (a + b)/2
        fc = f(c)
        if abs(fc) < tol
            return c
        end
        if sign(fc) == sign(fa)
            a = c
            fa = fc
        else
            b = c
        end
        if verbose
            @show iter, fc
        end
    end
end
