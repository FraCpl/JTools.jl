# Compute Empirical Cumulative Distribution Function F(y) of input data 'y'.
#
# Author: F. Capolupo
# European Space Agency, 2021
function ecdf(y::Vector)

    # Check inputs
    x = copy(y)
    filter!(x -> ~(isnan.(x) .| isinf.(x)), x)   # Remove NaN and Inf entries

    if isempty(x)
        return NaN, NaN
    end

    # Compute Empirical F(x)
    sort!(x)
    n = length(x)
    F = (1:n)./n

    return x, F
end

# Compute F[xi]
function interpxEcdf(x, F, xi)
    if xi > x[end]; return 1.0; end
    if xi < x[1]; return 0.0; end
    return interpfEcdf(F, x, xi)
end

# Compute x[Fi]
function interpfEcdf(x, F, Fi)
    id1 = findfirst(F .> Fi)
    F0 = F[id1 - 1]
    F1 = F[id1]
    x0 = x[id1 - 1]
    x1 = x[id1]
    return x0 + (x1 - x0)/(F1 - F0)*(Fi - F0)
end
