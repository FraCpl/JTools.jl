@inline crossmat(v) = [0.0 -v[3]  v[2]; v[3] 0.0 -v[1]; -v[2] v[1] 0.0];
@inline logrange(a, b, length=100) = exp10.(range(a; stop=b, length=length))

@inline function unwrap!(x, period=2π)
	y = convert(eltype(x), period)
	v = first(x)
	@inbounds for k = eachindex(x)
		x[k] = v = v + rem(x[k] - v,  y, RoundNearest)
	end
end

# Robust modulo operation (Matlab-like)
# modup=false (default): modd(x, x) = 0
# modup=true: modd(x, x) = x
@inline function modd(x, y; modup=false)
    r = rem(x, y, RoundNearest)
    # if isapprox(rem(x, y, RoundNearest), 0.0; atol=1.49e-8)     # 1.49e-8 = sqrt(eps())
    if abs(r) < 1.49e-8
        return modup ? y : 0.0
    end
    return mod(x, y)
end

# This function returns true if x is an integer multiple of y
@inline isMultiple(x, y) = modd(x, y) == 0.0

@inline signum(x) = x ≥ 0.0 ? 1.0 : -1.0

# Returns interpolating coefficients C, so that
# y ≈ C[1] + C[2]*x + C[3]*x^2 + ...
@inline function polyfit(x, y, n)
    A = ones(length(x), n + 1)
    @inbounds for i in 1:n
        A[:, i+1] = x.^i
    end
    return A\y
end

# Evaluate polynomial from coefficients C, so that
# y = C[1] + C[2]*x + C[3]*x^2 + ...
@inline function polyval(C, x)
    y = 0.0
    for i in eachindex(C)
        y += C[i]*x^(i - 1)
    end
    return y
end

@inline function interp1(x, y, xi)
    if xi ≥ x[end]
        return y[end]
    elseif xi < x[1]
        return y[1]
    end
    # id0 = sum(xi .≥ x)
    id0 = searchsortedlast(x, xi)
    dx = (xi - x[id0])/(x[id0 + 1] - x[id0])
    return (1 - dx).*y[id0] + dx.*y[id0 + 1]
end
