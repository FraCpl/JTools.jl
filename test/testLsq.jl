using JTools
using GLMakie

d = collect(range(0, 6, 300))
xTrue = [1.3; 0.4]
xEst0 = ones(length(xTrue))

fun(x) = 1/x[2]^3*exp.(-d*x[1])
y = fun(xTrue) + 0.5*randn(size(d))
f(x) = fun(x) - y

xEst, P = lsq(f, xEst0)
@show xEst

fig, ax = scatter(d, y)
lines!(ax, d, fun(xEst); color=:red, linewidth=3)
display(fig)
