using JTools

f(x) = [3x[1]^3+2-sin(x[1]); exp(x[2])-10x[2]]
x0 = [-1.0; 0.0]

x1 = rootFinder(f, x0; method = :NewtonRaphson, derivatives = :ForwardDiff)
x2 = rootFinder(f, x0; method = :NewtonRaphson, derivatives = :FiniteDiff)
x3 = rootFinder(f, x0; method = :Broyden)
@time x4 = rootFinder(f, x0; method = :ModifiedBroyden)
