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
