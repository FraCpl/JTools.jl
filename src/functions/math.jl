crossmat(v) = [0.0 -v[3]  v[2]; v[3] 0.0 -v[1]; -v[2] v[1] 0.0];
logrange(a, b, length=100) = exp10.(range(a; stop=b, length=length))
mag2db(x) = 20log10(x)
db2mag(x) = 10^(x/20)

function unwrap!(x, period=2π)
	y = convert(eltype(x), period)
	v = first(x)
	@inbounds for k = eachindex(x)
		x[k] = v = v + rem(x[k] - v,  y, RoundNearest)
	end
end

# Robust modulo operation (Matlab-like)
# modup=false (default): modd(x, x) = 0
# modup=true: modd(x, x) = x
function modd(x, y; modup=false)
    if isapprox(rem(x, y, RoundNearest), 0.0; atol=1.49e-8)     # 1.49e-8 = sqrt(eps())
        return modup ? y : 0.0
    end
    return mod(x, y)
end

# This function returns true if x is an integer multiple of y
isMultiple(x, y) = modd(x, y) == 0.0

signum(x) = x ≥ 0.0 ? 1.0 : -1.0
