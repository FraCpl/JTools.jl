using LinearAlgebra
using JTools
using BenchmarkTools

a = randn()
b = randn()

x = collect(range(0, 30, 10))
y = a*x .+ b
y[3] *=100
on = ones(2)
A = [x ones(length(x))]

function distFun(id)
    a, b = A[id, :]\y[id]
    return abs.(a*x - y .+ b)./sqrt(a*a + 1)  # orthogonal distance to line
end

@btime ransac(distFun, length(x), 2; verbose=false)
