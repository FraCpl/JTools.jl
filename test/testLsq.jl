using JTools
using GLMakie

function main()
    t = collect(range(0, 3, 100))
    xTrue = [1.3; 0.3]
    xEst0 = [4.0; 1.0]

    fun(t, x) = 1/x[2]^3*exp(-t*x[1])
    y = fun.(t, Ref(xTrue)) + 0.8*randn(size(t))
    W(x) = lsqWeight(x, 0.0, 5.0)

    # Add outliers
    y[50] += 50
    y[70] -= 30

    f(x) = fun.(t, Ref(x)) - y
    @time xEst1, ~ = lsq(f, xEst0; W=W, verbose=false)
    res = [x -> [fun(t[k], x) - y[k]] for k in eachindex(t)]
    @time xEst2, ~ = lsq(res, xEst0; W=W, verbose=false)

    @show xEst1
    @show xEst2

    fig, ax = scatter(t, y); display(fig)
    lines!(ax, t, fun.(t, Ref(xTrue)); linewidth=3, color=:black)
    lines!(ax, t, fun.(t, Ref(xEst1)); linewidth=3, color=:green)
    lines!(ax, t, fun.(t, Ref(xEst2)); linewidth=3, color=:red, linestyle=:dash)
end
main()


# fun(x) = [10.0*(x[2] - x[1]^2); 1.0 - x[1]]
# x0 = [-2.0; 2.0]
# xEst = lsq(fun, x0)[1]
