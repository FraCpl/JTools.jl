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
