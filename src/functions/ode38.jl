# f!(dx, t, x)
function ode38(f!, x0, tspan; nSteps = 100)
    # 3/8 Runge-Kutta Method
    # http://www.mymathlib.com/diffeq/runge-kutta/runge_kutta_3_8.html
    x = copy(x0)
    t = tspan[1]
    h = (tspan[2] - tspan[1])/nSteps
    K1 = similar(x);
    K2 = similar(x)
    K3 = similar(x);
    K4 = similar(x)
    for _ = 1:nSteps
        f!(K1, t, x)
        f!(K2, t + 1/3*h, x + h*K1/3)
        f!(K3, t + 2/3*h, x + h*(-K1/3 + K2))
        f!(K4, t + h, x + h*(K1 - K2 + K3))
        t += h
        x += h*(K1 + 3K2 + 3K3 + K4)/8
    end
    return x
end
