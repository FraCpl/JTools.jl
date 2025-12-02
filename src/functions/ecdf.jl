# Compute Empirical Cumulative Distribution Function F(y) of input data 'y'.
#
# Author: F. Capolupo
# European Space Agency, 2021
function ecdf(y::Vector)

    # Check inputs
    n = length(y)
    if n < 2 | any(isnan.(y) .| isinf.(y))
        return NaN, NaN
    end

    # Return Empirical F(x)
    return sort(y), (1:n) ./ n         # x, F(x)
end

# Compute F[xi]
function interpxEcdf(x, F, xi)
    if xi > x[end]
        ;
        return 1.0;
    end
    if xi < x[1]
        ;
        return 0.0;
    end
    return interpfEcdf(F, x, xi)
end

# Compute x[Fi]
function interpfEcdf(x, F, Fi)
    if Fi ≥ 1.0
        ;
        return x[end];
    end
    id1 = findfirst(F .> Fi)
    if id1 == 1
        ;
        return NaN;
    end        # TODO: CHECK!
    F0 = F[id1-1]
    F1 = F[id1]
    x0 = x[id1-1]
    x1 = x[id1]
    return x0 + (x1 - x0)/(F1 - F0)*(Fi - F0)
end
