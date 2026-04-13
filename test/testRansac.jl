using LinearAlgebra
using JTools
using BenchmarkTools

function distFun!(dist, x, y, id)        # non allocating....
    # Fit line parameters with data[id]
    @assert length(id) == 2
    i1, i2 = id
    x1, x2 = x[i1], x[i2]
    y1, y2 = y[i1], y[i2]
    invd = 1 / (x1 - x2)

    a = invd * (y1 - y2)
    b = invd * (x1 * y2 - x2 * y1)

    # Compute distance for all data points
    dd = 1 / sqrt(a *a + 1)
    @inbounds for i in eachindex(dist)
        dist[i] = abs(a * x[i] - y[i] + b) * dd # orthogonal distance to line
    end
    return dist
end

function testRansac()
    a = randn()
    b = randn()

    x = collect(range(0, 30, 10))
    y = a*x .+ b
    y[3] *= 100
    y[4] += 200
    dist = zeros(length(x))

    idtest = [2; 3]
    @btime distFun!($dist, $x, $y, $idtest)
    f(id) = distFun!(dist, x, y, id)

    Nmin = 2; N = length(x)

    @show ransac(f, N, Nmin; verbose=true)

    @btime ransac($f, $N, $Nmin; verbose=false)
end

testRansac()
