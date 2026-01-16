using JTools

xTrue = [1.0; 3.0]
f(x) = (x[1] + 2x[2] - 7)^2 + (2x[1] + x[2] - 5)^2

# xTrue = [3.0; 0.5]
# function f(x)
#     (1.5 - x[1] + x[1]*x[2])^2 + (2.25 - x[1] + x[1]*x[2]^2)^2 + (2.625 - x[1] + x[1]*x[2]^3)^2
# end

# f(x) = ()

x0 = [0.0; 0.0]
@show x = gradientDescent(f, x0; optimizer=Grad(η=0.1))
@show x = gradientDescent(f, x0; optimizer=RMSProp(η=0.1))
@show x = gradientDescent(f, x0; optimizer=Adam(η=0.1))
@show x = gradientDescent(f, x0; optimizer=Momentum(η=0.1))
# @show x = gradientDescent(f, x0; η=0.1, method=:gradls, verbose=false, maxIter=1000)
# @show x = gradientDescent(f, x0; η=0.1, method=:newton, verbose=false, maxIter=1000)
