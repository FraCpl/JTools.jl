# f!(dx, t, x)
@inline function ode38(f!, x0, tspan; nSteps=100, K1=similar(x0), K2=similar(x0), K3=similar(x0), K4=similar(x0), K5=similar(x0))
    # 3/8 Runge-Kutta Method
    # http://www.mymathlib.com/diffeq/runge-kutta/runge_kutta_3_8.html
    x = copy(x0)
    return ode38!(x, f!, x0, tspan, K1, K2, K3, K4, K5; nSteps=nSteps)
end

function ode38!(x, f!, x0, tspan, K1, K2, K3, K4, K5; nSteps=100)
    # 3/8 Runge-Kutta Method
    # http://www.mymathlib.com/diffeq/runge-kutta/runge_kutta_3_8.html
    x .= x0
    t = tspan[1]
    h = (tspan[2] - tspan[1])/nSteps
    for _ in 1:nSteps
        f!(K1, t, x)

        @. K5 = x + h*K1/3
        f!(K2, t + 1/3*h, K5)

        @. K5 = x + h*(-K1/3 + K2)
        f!(K3, t + 2/3*h, K5)

        @. K5 = x + h*(K1 - K2 + K3)
        f!(K4, t + h, K5)

        t += h
        @. x += h*(K1 + 3K2 + 3K3 + K4)/8
    end
    return x
end
