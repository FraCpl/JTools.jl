# Find minimum of unimodal scalar function f(x) of scalar variable x, with x ∈ [lb, ub]
function goldenSectionSearch(f, lb, ub; tol=1e-6, maxIter=1000, verbose=false)
    invϕ = 2 / (sqrt(5) + 1)
    a, b = lb, ub

    c = b - invϕ * (b - a)
    d = a + invϕ * (b - a)

    fc, fd = f(c), f(d)
    for iter in 1:maxIter
        if (b - a) ≤ tol
            break
        end

        if fc < fd
            b, d = d, c
            c = b - invϕ * (b - a)
            fd = fc
            fc = f(c)
        else
            a, c = c, d
            d = a + invϕ * (b - a)
            fc = fd
            fd = f(d)
        end

        if verbose
            @show iter, min(fc, fd)
        end
    end

    return (a + b) / 2
end
