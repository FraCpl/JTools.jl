# f!(dx, t, x)
@inline function ode38(f!, x0, tspan; nSteps=100, K1=similar(x0), K2=similar(x0), K3=similar(x0), K4=similar(x0), Ktmp=similar(x0))
    # 3/8 Runge-Kutta Method
    # http://www.mymathlib.com/diffeq/runge-kutta/runge_kutta_3_8.html
    x = copy(x0)
    return ode38!(x, f!, x0, tspan, K1, K2, K3, K4, Ktmp; nSteps=nSteps)
end

function ode38!(x, f!, x0, tspan, K1, K2, K3, K4, Ktmp; nSteps=100)
    # 3/8 Runge-Kutta Method
    # http://www.mymathlib.com/diffeq/runge-kutta/runge_kutta_3_8.html
    x .= x0
    t = tspan[1]
    h = (tspan[2] - tspan[1]) / nSteps
    for _ in 1:nSteps
        f!(K1, t, x)

        @. Ktmp = x + h / 3 * K1
        f!(K2, t + h / 3, Ktmp)

        @. Ktmp = x - h / 3 * K1 + h * K2
        f!(K3, t + 2 / 3 * h, Ktmp)

        @. Ktmp = x + h * (K1 - K2 + K3)
        f!(K4, t + h, Ktmp)

        t += h
        @. x += h / 8 * (K1 + 3 * K2 + 3 * K3 + K4)
    end
    return x
end

# f!(dx, t, x, p)
@inline function ode38(f!, x0, tspan, p; nSteps=100, K1=similar(x0), K2=similar(x0), K3=similar(x0), K4=similar(x0), Ktmp=similar(x0))
    # 3/8 Runge-Kutta Method
    # http://www.mymathlib.com/diffeq/runge-kutta/runge_kutta_3_8.html
    x = copy(x0)
    return ode38!(x, f!, x0, tspan, p, K1, K2, K3, K4, Ktmp; nSteps=nSteps)
end

@inline function ode38!(x, f!, x0, tspan, p, K1, K2, K3, K4, Ktmp; nSteps=100)
    # 3/8 Runge-Kutta Method
    # http://www.mymathlib.com/diffeq/runge-kutta/runge_kutta_3_8.html
    x .= x0
    t = tspan[1]
    h = (tspan[2] - tspan[1]) / nSteps
    for _ in 1:nSteps
        f!(K1, t, x, p)

        @. Ktmp = x + h / 3 * K1
        f!(K2, t + h / 3, Ktmp, p)

        @. Ktmp = x - h / 3 * K1 + h * K2
        f!(K3, t + 2 / 3 * h, Ktmp, p)

        @. Ktmp = x + h * (K1 - K2 + K3)
        f!(K4, t + h, Ktmp, p)

        t += h
        @. x += h / 8 * (K1 + 3 * K2 + 3 * K3 + K4)
    end
    return x
end
