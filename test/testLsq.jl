using JTools
using GLMakie

function main()
    t = collect(range(-1.1, 2.5, 100))
    xTrue = [1.3; 0.3]
    xEst0 = [4.0; 1.0]

    fun(t, x) = x[2] + x[1] * t * t #1 / (x[2]^3 + 0.01) * exp(-t * x[1])
    y = fun.(t, Ref(xTrue)) + 0.8*randn(size(t))
    W(x) = lsqWeight(x, 0.0, 5.0)

    # Add outliers
    y[50] += 50
    y[70] -= 30

    f(x) = fun.(t, Ref(x)) - y
    res = [x -> [fun(t[k], x) - y[k]] for k in eachindex(t)]
    resh = [x -> (fun(t[k], x) - y[k], [t[k] * t[k] 1]) for k in eachindex(t)]  # With jacobian

    @time xEst1, _ = lsq(f, xEst0; W=W, verbose=true)
    @time xEst2, _ = lsq(res, xEst0; W=W, verbose=false)
    @time xEst3, _ = lsq(f, xEst0; W=W, λ=0.0, verbose=true)
    @time xEst4, _ = lsq(resh, xEst0; W=W, verbose=true, userJacobian=true)

    @show xEst1
    @show xEst2
    @show xEst3
    @show xEst4

    fig, ax = scatter(t, y); display(fig)
    lines!(ax, t, fun.(t, Ref(xTrue)); linewidth=3, color=:black)
    lines!(ax, t, fun.(t, Ref(xEst1)); linewidth=3, color=:green)
    lines!(ax, t, fun.(t, Ref(xEst2)); linewidth=3, color=:red, linestyle=:dash)
    lines!(ax, t, fun.(t, Ref(xEst3)); linewidth=3, color=:orange, linestyle=:dot)
    lines!(ax, t, fun.(t, Ref(xEst4)); linewidth=3, color=:cyan, linestyle=:dash)
end
main()

# fun(x) = [10.0*(x[2] - x[1]^2); 1.0 - x[1]]
# x0 = [-2.0; 2.0]
# xEst = lsq(fun, x0)[1]
