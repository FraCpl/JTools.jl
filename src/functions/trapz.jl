function trapz(t, x)
    r = 0.0.*x[1]
    for k in 2:lastindex(t)
        r += 0.5*(x[k] + x[k-1])*(t[k] - t[k-1])
    end
    return r
end

function cumtrapz(t, x)
    r = copy(x)
    r[1] .= 0.0
    for k in 2:lastindex(t)
        r[k] = r[k-1] + 0.5*(x[k] + x[k-1])*(t[k] - t[k-1])
    end
    return r
end
