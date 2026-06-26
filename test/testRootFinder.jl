using JTools, Test

function testRootFinder()
    function f!(y, x)
        y[1] = 3 * x[1]^3 + 2 - sin(x[1])
        y[2] = exp(x[2]) - 10 * x[2]
        return y
    end

    f(x) = f!(zero(x), x)

    x0 = [-1.0; 0.0]
    y = zeros(2)

    @testset "rootFinder.jl" begin
        for d in (:ForwardDiff, :FiniteDiff), m in (:NewtonRaphson, :Broyden, :ModifiedBroyden)
            x1 = rootFinder(f, x0; method=m, derivatives=d)
            x2 = rootFinder!(y, f!, x0; method=m, derivatives=d)
            x3 = rootFinder(f, x0; dxMax=1e-1*ones(2), method=m, derivatives=d)
            x4 = rootFinder!(y, f!, x0; dxMax=1e-1*ones(2), method=m, derivatives=d)

            @test norm(f(x1)) < 1e-10
            @test norm(f(x2)) < 1e-10
            @test norm(f(x3)) < 1e-9
            @test norm(f(x4)) < 1e-9
        end
    end
end

testRootFinder()
