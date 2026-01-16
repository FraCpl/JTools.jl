abstract type AbstractLocalOptimizer end

# GRADIENT OPTIMIZER ----------------------------------- #
struct Grad <: AbstractLocalOptimizer
    η::Float64
end

Grad(; η=1.0) = Grad(η)

initOptimizer!(opt::Grad, nx) = nothing

function stepOptimizer!(opt::Grad, δx, ∇f)
    @inbounds for i in eachindex(δx)
        δx[i] = -opt.η*∇f[i]
    end
    return nothing
end

# RMSPROP OPTIMIZER ----------------------------------- #
mutable struct RMSProp <: AbstractLocalOptimizer
    η::Float64
    ρ::Float64
    mx::Vector{Float64}
end

RMSProp(; η=1.0, ρ=0.9) = RMSProp(η, ρ, Float64[])

function initOptimizer!(opt::RMSProp, nx)
    opt.mx = zeros(nx)
    return nothing
end

function stepOptimizer!(opt::RMSProp, δx, ∇f)
    mx = opt.mx
    @inbounds for i in eachindex(δx)
        mx[i] = opt.ρ * mx[i] + (1.0 - opt.ρ) * ∇f[i] * ∇f[i]
        δx[i] = -opt.η / (sqrt(mx[i]) + 1e-8) * ∇f[i]
    end
    return nothing
end

# MOMENTUM OPTIMIZER ----------------------------------#
# https://www.cs.toronto.edu/%7Ehinton/absps/momentum.pdf
mutable struct Momentum <: AbstractLocalOptimizer
    η::Float64
    μ::Float64
    mx::Vector{Float64}
end

Momentum(; η=1.0, μ=0.3) = Momentum(η, μ, Float64[])

function initOptimizer!(opt::Momentum, nx)
    opt.mx = zeros(nx)
    return nothing
end

function stepOptimizer!(opt::Momentum, δx, ∇f)
    mx = opt.mx
    @inbounds for i in eachindex(δx)
        mx[i] = opt.μ * mx[i] - opt.η * ∇f[i]
        δx[i] = mx[i]
    end
    return nothing
end

# Adam OPTIMIZER ----------------------------------#
mutable struct Adam <: AbstractLocalOptimizer
    η::Float64
    β1::Float64
    β2::Float64
    b1::Float64
    b2::Float64
    mx::Vector{Float64}
    vx::Vector{Float64}
end

Adam(; η=1.0, β1=0.9, β2=0.999) = Adam(η, β1, β2, 1.0, 1.0, Float64[], Float64[])

function initOptimizer!(opt::Adam, nx)
    opt.mx = zeros(nx)
    opt.vx = zeros(nx)
    opt.b1 = 1.0
    opt.b2 = 1.0
    return nothing
end

function stepOptimizer!(opt::Adam, δx, ∇f)
    # [1] Ruder, An overview of Gradient Descent Optimization Algorithms
    # [2] https://arxiv.org/pdf/1412.6980.pdf
    mx = opt.mx
    vx = opt.vx
    opt.b1 = opt.b1*opt.β1
    opt.b2 = opt.b2*opt.β2
    @inbounds for i in eachindex(δx)
        mx[i] = opt.β1 * mx[i] + (1 - opt.β1) * ∇f[i]
        vx[i] = opt.β2 * vx[i] + (1 - opt.β2) * ∇f[i] * ∇f[i]
        δx[i] = -opt.η / (sqrt(vx[i] / (1 - opt.b2)) + 1e-8) * mx[i] / (1 - opt.b1)
    end
    return nothing
end


# gradient!(∇f, f, x) with y = f(x)
#
# examples:
# gradient! = FiniteDiff.finite_difference_gradient!
# gradient! = ForwardDiff.gradient!
function gradientDescent(f, x0; optimizer::T=Grad(), gradient!::F=FiniteDiff.finite_difference_gradient!, maxIter=300, dxMax=Inf*ones(length(x0)), tol=1e-6, verbose=false) where {T<:AbstractLocalOptimizer, F}
    # Initialize problem and optimizer
    x = copy(x0)
    nx = length(x)
    δx = zeros(nx)
    ∇f = zeros(nx)
    # H = zeros(nx, nx)
    # mx = zeros(nx)
    # vx = zeros(nx)

    initOptimizer!(optimizer, nx)
    tolSq = tol*tol

    for iter in 1:maxIter
        gradient!(∇f, f, x)                 # Compute gradient
        stepOptimizer!(optimizer, δx, ∇f)   # Call optimizer step

        # elseif method == :gradls
        #     # Gradient descent with line search
        #     ηOpt = goldenSectionSearch(n -> f(x - n*∇f), 0.0, 10η; verbose=false)
        #     @inbounds for i in eachindex(δx)
        #         δx[i] = -ηOpt*∇f[i]
        #     end

        # elseif method == :newton    # {CHECK: IS THIS WORKING?}
        #     if derivatives == :ForwardDiff
        #         H .= ForwardDiff.hessian(f, x)
        #     else
        #         H .= FiniteDiff.finite_difference_hessian(f, x)
        #     end
        #     δx .= -H\∇f


        # Correct solution (check bounds)
        err = 0.0
        for i in eachindex(δx)
            dx = clamp(δx[i], -dxMax[i], dxMax[i])
            x[i] += dx
            err += dx*dx
        end

        if verbose
            @show iter, f(x)
        end
        if isnan(err) || err ≤ tolSq
            break
        end
    end

    # Return solution
    return x
end
