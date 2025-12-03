# Find minimum of unimodal scalar function f(x) of scalar variable x, with x ∈ [lb, ub]
function goldenSectionSearch(f::Function, lb, ub; tol=1e-6, maxIter=1000, verbose=true)

    # Initialize variables
    invϕ = 2/(√5 + 1)
    a = lb
    b = ub
    c = b - (b - a)*invϕ
    d = a + (b - a)*invϕ

    # Start iterations
    for iter in 1:maxIter
        if abs(c - d) ≤ tol
            break
        end
        fc = f(c)
        fd = f(d)
        if fc < fd
            b = d
        else
            a = c
        end

        c = b - (b - a)*invϕ
        d = a + (b - a)*invϕ
        if verbose
            @show iter, min(fc, fd)
        end
    end

    # Return solution
    a = (a + b)/2
    if verbose
        @show "Solution", a, f(a)
    end
    return a
end
